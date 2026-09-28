//! BAM access through noodles: record accessors, indexed region queries, and writing.
//!
//! Target IDs are exposed as `i32` with `-1` for "no reference", and positions as
//! 0-based `i64` with `-1` for "no position", matching SAM/BAM on-disk conventions.

use crate::types::GenomicRegion;
use noodles::bam;
use noodles::core::{Position, Region};
use noodles::csi::{self, BinningIndex};
use noodles::sam::{
    self,
    alignment::{
        io::Write as _,
        record::{cigar::op::Kind, data::field::Tag},
        record_buf::data::field::Value,
    },
};
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

/// A decoded alignment record.
pub type Record = sam::alignment::RecordBuf;
/// A SAM/BAM header.
pub type Header = sam::Header;
/// A decoded CIGAR.
pub type Cigar = sam::alignment::record_buf::Cigar;

/// Quality value reported for every base when a record has no qualities (SAM `*`).
const MISSING_QUALITY: u8 = 0xff;

/// Coordinate and flag accessors in BAM's signed conventions.
pub trait RecordExt {
    /// Reference sequence ID, or `-1` when unplaced.
    fn tid(&self) -> i32;
    /// 0-based alignment start, or `-1` when absent.
    fn pos(&self) -> i64;
    /// Mate reference sequence ID, or `-1` when unplaced.
    fn mtid(&self) -> i32;
    /// 0-based mate alignment start, or `-1` when absent.
    fn mpos(&self) -> i64;
    /// Mapping quality, `255` when unavailable.
    fn mapq(&self) -> u8;
    /// Raw SAM flag bits.
    fn flag_bits(&self) -> u16;
    /// Exclusive 0-based reference end implied by the CIGAR.
    fn end_pos(&self) -> i64;
    /// Base qualities, with a `0xff` per base when the record stores none.
    fn qual(&self) -> Cow<'_, [u8]>;
}

fn signed_id(id: Option<usize>) -> i32 {
    id.map_or(-1, |id| id as i32)
}

fn zero_based(position: Option<Position>) -> i64 {
    position.map_or(-1, |p| usize::from(p) as i64 - 1)
}

impl RecordExt for Record {
    fn tid(&self) -> i32 {
        signed_id(self.reference_sequence_id())
    }

    fn pos(&self) -> i64 {
        zero_based(self.alignment_start())
    }

    fn mtid(&self) -> i32 {
        signed_id(self.mate_reference_sequence_id())
    }

    fn mpos(&self) -> i64 {
        zero_based(self.mate_alignment_start())
    }

    fn mapq(&self) -> u8 {
        self.mapping_quality().map_or(255, u8::from)
    }

    fn flag_bits(&self) -> u16 {
        self.flags().bits()
    }

    fn end_pos(&self) -> i64 {
        crate::utils::calculate_end_pos(self.pos(), self.cigar())
    }

    fn qual(&self) -> Cow<'_, [u8]> {
        let qual = self.quality_scores().as_ref();
        let seq_len = self.sequence().len();
        if qual.is_empty() && seq_len > 0 {
            Cow::Owned(vec![MISSING_QUALITY; seq_len])
        } else {
            Cow::Borrowed(qual)
        }
    }
}

/// Whether a CIGAR op kind consumes reference bases.
pub fn consumes_reference(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip
    )
}

/// Replace a record's base qualities.
pub fn set_quality_scores(record: &mut Record, qual: Vec<u8>) {
    *record.quality_scores_mut() = qual.into();
}

/// The `MC` (mate CIGAR) tag as a string, when present.
pub fn mate_cigar(record: &Record) -> Option<&str> {
    match record.data().get(&Tag::MATE_CIGAR) {
        Some(Value::String(mc)) => std::str::from_utf8(mc).ok(),
        _ => None,
    }
}

/// Query name bytes (empty when absent).
pub fn qname(record: &Record) -> &[u8] {
    record.name().map_or(&[], |name| name.as_ref())
}

/// `(name, length)` of every reference sequence, in target-ID order.
pub fn reference_sequences(header: &Header) -> Vec<(String, usize)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, rs)| (name.to_string(), rs.length().get()))
        .collect()
}

/// Target ID to reference name mapping.
pub fn tid_to_name(header: &Header) -> HashMap<i32, String> {
    reference_sequences(header)
        .into_iter()
        .enumerate()
        .map(|(tid, (name, _))| (tid as i32, name))
        .collect()
}

/// Parse one SAM record line against `header`.
pub fn record_from_sam(header: &Header, line: &[u8]) -> io::Result<Record> {
    let mut reader = sam::io::Reader::new(line);
    let mut record = sam::Record::default();
    reader.read_record(&mut record)?;
    Record::try_from_alignment_record(header, &record)
}

/// A sequential BAM reader that decodes records.
pub struct BamReader {
    reader: bam::io::Reader<noodles::bgzf::io::Reader<File>>,
    header: Header,
}

impl BamReader {
    /// Open `path` and read its header.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let mut reader = bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        Ok(Self { reader, header })
    }

    /// The BAM header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Decoded records in file order.
    pub fn records(&mut self) -> impl Iterator<Item = io::Result<Record>> + '_ {
        self.reader.record_bufs(&self.header)
    }
}

/// A BAM file with its header and index loaded once, queried by region from any thread.
pub struct IndexedBam {
    path: PathBuf,
    header: Header,
    index: Box<dyn BinningIndex + Send + Sync>,
}

impl IndexedBam {
    /// Open `path`, reading its header and its `.bai` or `.csi` index.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let header = BamReader::open(&path)?.header;
        let index = read_index(&path)?;
        Ok(Self {
            path,
            header,
            index,
        })
    }

    /// The BAM header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Call `on_record` for every record overlapping `region` (0-based, half-open).
    pub fn for_each_in_region(
        &self,
        region: &GenomicRegion,
        mut on_record: impl FnMut(Record),
    ) -> io::Result<()> {
        if region.start >= region.end {
            return Ok(());
        }
        let name = usize::try_from(region.tid)
            .ok()
            .and_then(|tid| self.header.reference_sequences().get_index(tid))
            .map(|(name, _)| name.to_string())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown target ID {}", region.tid),
                )
            })?;
        let to_position = |p: i64| {
            usize::try_from(p)
                .ok()
                .and_then(Position::new)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid region"))
        };
        let query_region = Region::new(
            name,
            to_position(region.start + 1)?..=to_position(region.end)?,
        );

        let mut reader = File::open(&self.path).map(bam::io::Reader::new)?;
        let query = reader.query(&self.header, &self.index, &query_region)?;
        for result in query.records() {
            let record = result?;
            on_record(Record::try_from_alignment_record(&self.header, &record)?);
        }
        Ok(())
    }
}

fn read_index(bam_path: &Path) -> io::Result<Box<dyn BinningIndex + Send + Sync>> {
    let with_suffix = |suffix: &str| {
        let mut name = bam_path.as_os_str().to_owned();
        name.push(suffix);
        PathBuf::from(name)
    };
    let bai = with_suffix(".bai");
    if bai.exists() {
        return Ok(Box::new(bam::bai::fs::read(bai)?));
    }
    let csi = with_suffix(".csi");
    if csi.exists() {
        return Ok(Box::new(csi::fs::read(csi)?));
    }
    let sibling_bai = bam_path.with_extension("bai");
    if sibling_bai.exists() {
        return Ok(Box::new(bam::bai::fs::read(sibling_bai)?));
    }
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!("no .bai or .csi index found for {}", bam_path.display()),
    ))
}

/// A BAM writer; call [`BamWriter::finish`] to flush the final BGZF block and EOF marker.
pub struct BamWriter {
    writer: bam::io::Writer<noodles::bgzf::io::Writer<Box<dyn Write>>>,
    header: Header,
}

impl BamWriter {
    /// Create a writer to `path`, or to stdout when `path` is `None`, and write `header`.
    pub fn create(path: Option<&Path>, header: Header) -> io::Result<Self> {
        let sink: Box<dyn Write> = match path {
            Some(path) => Box::new(File::create(path)?),
            None => Box::new(io::stdout().lock()),
        };
        let mut writer = bam::io::Writer::new(sink);
        writer.write_header(&header)?;
        Ok(Self { writer, header })
    }

    /// Write one record.
    pub fn write(&mut self, record: &Record) -> io::Result<()> {
        self.writer.write_alignment_record(&self.header, record)
    }

    /// Flush and finish the BGZF stream.
    pub fn finish(mut self) -> io::Result<()> {
        self.writer.try_finish()
    }
}
