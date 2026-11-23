//! Optical-duplicate comparisons.
//!
//! An optical duplicate and its original are two clusters grown from the same library
//! molecule, so a base where they disagree points at the sequencer (cluster amplification,
//! imaging, base calling) rather than the library. This module pairs each optical duplicate
//! with its original and counts every compared base by read cycle.

use crate::processing::{read_axis_position, record_read_num};
use crate::types::{OVERLAP_QUALITY_BIN_WIDTH, ProcessingConfig};
use crate::utils::{base_to_char, complement};
use rust_htslib::bam::Record;
use rust_htslib::bam::record::{Aux, CigarString, CigarStringView};
use std::collections::HashMap;
use std::io::Write;

/// How a record was marked by duplicate marking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DuplicateKind {
    /// Not flagged as a duplicate (0x400 unset).
    NotDuplicate,
    /// Flagged, with `DT:Z:SQ` (Picard) or `dt:Z:SQ` (samtools markdup -d).
    Optical,
    /// Flagged, with `DT:Z:LB` / `dt:Z:LB`: a PCR/library duplicate.
    Library,
    /// Flagged with no (recognised) duplicate-type tag.
    Untagged,
}

/// Classify a record by its duplicate flag and `DT`/`dt` tag.
pub fn duplicate_kind(record: &Record) -> DuplicateKind {
    if !record.is_duplicate() {
        return DuplicateKind::NotDuplicate;
    }
    for tag in [b"DT", b"dt"] {
        if let Ok(Aux::String(value)) = record.aux(tag) {
            return match value {
                "SQ" => DuplicateKind::Optical,
                "LB" => DuplicateKind::Library,
                _ => DuplicateKind::Untagged,
            };
        }
    }
    DuplicateKind::Untagged
}

/// The coordinates duplicate markers group on: each mate's unclipped 5' end and strand.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct DuplicateSetKey {
    pub tid: i32,
    pub five_prime: i64,
    pub is_reverse: bool,
    pub read_num: u8,
    pub mate_tid: i32,
    pub mate_five_prime: i64,
    pub mate_is_reverse: bool,
}

fn unclipped_five_prime(cigar: &CigarStringView, is_reverse: bool) -> i64 {
    if is_reverse {
        cigar.end_pos() + cigar.trailing_softclips() + cigar.trailing_hardclips()
    } else {
        cigar.pos() - cigar.leading_softclips() - cigar.leading_hardclips()
    }
}

/// Duplicate-set key for a mapped read with a mapped mate. The mate's 5' end comes from
/// the `MC` tag; without it the mate's aligned start is used, which can split a set when
/// mates are clipped differently (pairs are then missed, never mismatched).
pub fn duplicate_set_key(record: &Record) -> Option<DuplicateSetKey> {
    if !record.is_paired() || record.is_unmapped() || record.is_mate_unmapped() {
        return None;
    }
    let mate_is_reverse = record.is_mate_reverse();
    let mate_five_prime = match record.aux(b"MC") {
        Ok(Aux::String(mc)) => {
            let view = CigarString::try_from(mc).ok()?.into_view(record.mpos());
            unclipped_five_prime(&view, mate_is_reverse)
        }
        _ => record.mpos(),
    };
    Some(DuplicateSetKey {
        tid: record.tid(),
        five_prime: unclipped_five_prime(&record.cigar(), record.is_reverse()),
        is_reverse: record.is_reverse(),
        read_num: record_read_num(record),
        mate_tid: record.mtid(),
        mate_five_prime,
        mate_is_reverse,
    })
}

/// Key for one compared base between an optical duplicate and its original.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct OpticalDuplicateKey {
    /// Read number by SAM flag (1 or 2).
    pub read_num: u8,
    /// Duplicate's position on the tasmanian-mismatch read-mode axis (1-based, from 5').
    pub position: usize,
    /// Original's base, in read orientation (`N` for anything other than A/C/G/T).
    pub original_base: char,
    /// Duplicate's base, in read orientation.
    pub duplicate_base: char,
    /// Lower bound of the bin holding min(original quality, duplicate quality).
    pub min_quality_bin: u8,
}

/// Duplicate sets seen and why some were not compared.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct OpticalDuplicateStats {
    /// Duplicate-flagged reads carrying an optical (`SQ`) tag.
    pub optical_reads: usize,
    /// Sets with at least one optical duplicate whose duplicates were compared.
    pub sets_compared: usize,
    /// Optical duplicates compared with their original.
    pub duplicates_compared: usize,
    /// Sets skipped because they also hold library or untagged duplicates.
    pub sets_skipped_mixed: usize,
    /// Sets skipped because no unflagged original passed the filters in this region.
    pub sets_skipped_no_original: usize,
    /// Sets skipped because more than one unflagged read shares the coordinates.
    pub sets_skipped_multiple_originals: usize,
}

impl std::ops::AddAssign for OpticalDuplicateStats {
    fn add_assign(&mut self, o: Self) {
        self.optical_reads += o.optical_reads;
        self.sets_compared += o.sets_compared;
        self.duplicates_compared += o.duplicates_compared;
        self.sets_skipped_mixed += o.sets_skipped_mixed;
        self.sets_skipped_no_original += o.sets_skipped_no_original;
        self.sets_skipped_multiple_originals += o.sets_skipped_multiple_originals;
    }
}

fn passes_filters(record: &Record, config: &ProcessingConfig) -> bool {
    !(record.is_unmapped()
        || record.is_secondary()
        || record.is_supplementary()
        || record.is_quality_check_failed()
        || record.mapq() < config.min_map_quality)
}

/// Compare an optical duplicate with its original cycle by cycle from the 5' end.
///
/// Both reads come from the same molecule and start at the same unclipped 5' end, so cycle
/// `i` of each read covers the same template base whatever the aligner did with clipping or
/// indels. Bases are written in read orientation; only bases where both reads pass
/// `min_base_quality` are counted, concordant ones included.
pub fn compare_optical_duplicate(
    original: &Record,
    duplicate: &Record,
    config: &ProcessingConfig,
    counts: &mut HashMap<OpticalDuplicateKey, usize>,
) {
    let is_reverse = duplicate.is_reverse();
    let (orig_seq, orig_qual) = (original.seq(), original.qual());
    let (dup_seq, dup_qual) = (duplicate.seq(), duplicate.qual());
    let (orig_len, dup_len) = (orig_seq.len(), dup_seq.len());
    let read_num = record_read_num(duplicate);
    let stored = |cycle: usize, len: usize| if is_reverse { len - 1 - cycle } else { cycle };
    let oriented = |byte: u8| {
        let c = base_to_char(byte).unwrap_or('N');
        if is_reverse { complement(c) } else { c }
    };

    for cycle in 0..orig_len.min(dup_len) {
        let (oi, di) = (stored(cycle, orig_len), stored(cycle, dup_len));
        let (oq, dq) = (orig_qual[oi], dup_qual[di]);
        if oq < config.min_base_quality || dq < config.min_base_quality {
            continue;
        }
        let min_qual = oq.min(dq);
        let key = OpticalDuplicateKey {
            read_num,
            position: read_axis_position(di, dup_len, is_reverse, config.mode_len),
            original_base: oriented(orig_seq[oi]),
            duplicate_base: oriented(dup_seq[di]),
            min_quality_bin: min_qual - min_qual % OVERLAP_QUALITY_BIN_WIDTH,
        };
        *counts.entry(key).or_insert(0) += 1;
    }
}

#[derive(Default)]
struct DuplicateSet<'a> {
    originals: Vec<&'a Record>,
    optical: Vec<&'a Record>,
    other_duplicates: usize,
}

/// Group records into duplicate sets and compare every optical duplicate with its original.
///
/// A set is used only when it has exactly one unflagged read (the original) and all its
/// duplicates are optical. A set that also holds library or untagged duplicates is skipped,
/// because the optical duplicate may be an optical copy of a PCR duplicate rather than of the
/// original, which would let library errors into the comparison. `-f`/`-F`/`-G` are not
/// applied here (they would usually drop the duplicates); mapping quality, unmapped,
/// secondary, supplementary and QC-fail reads are.
pub fn collect_optical_duplicate_comparisons<'a>(
    records: impl IntoIterator<Item = &'a Record>,
    config: &ProcessingConfig,
    counts: &mut HashMap<OpticalDuplicateKey, usize>,
    stats: &mut OpticalDuplicateStats,
) {
    let mut sets: HashMap<DuplicateSetKey, DuplicateSet<'a>> = HashMap::new();
    let mut originals: Vec<(DuplicateSetKey, &'a Record)> = Vec::new();

    for record in records {
        if !passes_filters(record, config) {
            continue;
        }
        let kind = duplicate_kind(record);
        let Some(key) = duplicate_set_key(record) else {
            continue;
        };
        match kind {
            DuplicateKind::NotDuplicate => originals.push((key, record)),
            DuplicateKind::Optical => {
                stats.optical_reads += 1;
                sets.entry(key).or_default().optical.push(record);
            }
            DuplicateKind::Library | DuplicateKind::Untagged => {
                sets.entry(key).or_default().other_duplicates += 1;
            }
        }
    }
    // Only originals whose coordinates hold a duplicate are needed.
    for (key, record) in originals {
        if let Some(set) = sets.get_mut(&key) {
            set.originals.push(record);
        }
    }

    for set in sets.values() {
        if set.optical.is_empty() {
            continue;
        }
        if set.other_duplicates > 0 {
            stats.sets_skipped_mixed += 1;
            continue;
        }
        let original = match set.originals.as_slice() {
            [] => {
                stats.sets_skipped_no_original += 1;
                continue;
            }
            [original] => *original,
            _ => {
                stats.sets_skipped_multiple_originals += 1;
                continue;
            }
        };
        stats.sets_compared += 1;
        for duplicate in &set.optical {
            compare_optical_duplicate(original, duplicate, config, counts);
            stats.duplicates_compared += 1;
        }
    }
}

/// Write optical-duplicate comparisons as TSV. Rows with equal bases are the concordant
/// comparisons, so `count` summed over rows with different bases, divided by `count` summed
/// over all rows, is the discordance rate.
pub fn write_optical_duplicates_to_writer<W: Write>(
    counts: &HashMap<OpticalDuplicateKey, usize>,
    writer: &mut W,
) -> std::io::Result<()> {
    let mut rows: Vec<_> = counts.iter().collect();
    rows.sort_unstable_by_key(|&(key, _)| key);
    writeln!(
        writer,
        "read_num\tposition\toriginal_base\tduplicate_base\tmin_quality_bin\tcount"
    )?;
    for (key, count) in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            key.read_num,
            key.position,
            key.original_base,
            key.duplicate_base,
            key.min_quality_bin,
            count
        )?;
    }
    Ok(())
}

/// Write optical-duplicate comparisons to a TSV file.
pub fn write_optical_duplicates_tsv(
    counts: &HashMap<OpticalDuplicateKey, usize>,
    output_path: &str,
) -> std::io::Result<()> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(output_path)?);
    write_optical_duplicates_to_writer(counts, &mut writer)?;
    writer.flush()
}
