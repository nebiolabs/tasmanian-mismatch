use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
use std::fs;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

pub fn unique_temp_dir(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("{}_{}", prefix, nanos));
    fs::create_dir_all(&dir).expect("failed to create temp dir");
    dir
}

pub fn repo_log_path(prefix: &str) -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let log_dir = manifest_dir.join("test_output/integration_logs");
    fs::create_dir_all(&log_dir).expect("failed to create integration log dir");
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    log_dir.join(format!("{}_{}.log", prefix, nanos))
}

pub fn log_line(log_path: &std::path::Path, message: &str) {
    let mut log = OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)
        .expect("failed to open integration log");
    writeln!(log, "{}", message).expect("failed to write integration log");
}

pub fn log_command(log_path: &std::path::Path, binary: &str, args: &[&str]) {
    let cmd_str = format!("{} {}", binary, args.join(" "));
    log_line(log_path, &format!("Executing: {}", cmd_str));
}

/// Write a BAM with two 8 bp contigs (chr1, chr2), one mismatched read on each.
#[allow(dead_code)] // not every test crate that includes this module uses it
pub fn write_two_chrom_bam(path: &Path) {
    let mut header = Header::new();
    for name in ["chr1", "chr2"] {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", name);
        sq.push_tag(b"LN", 8);
        header.push_record(&sq);
    }

    let header_view = HeaderView::from_header(&header);
    let mut writer =
        Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");

    // Both references are ACGTACGT. read1 on chr1 has A->T; read2 on chr2 has T->A.
    for sam_line in [
        &b"read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:1"[..],
        &b"read2\t0\tchr2\t1\t60\t8M\t*\t0\t0\tACGTACGA\tIIIIIIII\tNM:i:1"[..],
    ] {
        let record = Record::from_sam(&header_view, sam_line).expect("failed to parse SAM line");
        writer.write(&record).expect("failed to write BAM record");
    }
}
