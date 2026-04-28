use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::index;
use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_temp_dir(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("{}_{}", prefix, nanos));
    fs::create_dir_all(&dir).expect("failed to create temp dir");
    dir
}

fn write_test_bam(path: &Path) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", &"chr1");
    sq.push_tag(b"LN", &8);
    header.push_record(&sq);

    let header_view = HeaderView::from_header(&header);
    let mut writer = Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");

    // Reference is ACGTACGT. This read has one mismatch (A->T).
    let sam_line = b"read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:1";
    let record = Record::from_sam(&header_view, sam_line).expect("failed to parse SAM line");
    writer.write(&record).expect("failed to write BAM record");
}

#[test]
fn integration_diagnostics_fixture_bam_produces_expected_outputs() {
    let temp_dir = unique_temp_dir("diagnostics_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");
    let discounts_tsv = temp_dir.join("discounts.tsv");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
    write_test_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let binary = env!("CARGO_BIN_EXE_tasmanian-diagnostics");
    let status = Command::new(binary)
        .arg("-q")
        .arg("0")
        .arg("--min-map-quality")
        .arg("0")
        .arg("--genomic-threshold")
        .arg("1")
        .arg("--genomic-depth-threshold")
        .arg("1")
        .arg("--variants-output")
        .arg(&variants_tsv)
        .arg("--inconsistencies-output")
        .arg(&inconsistencies_tsv)
        .arg("--discount-output")
        .arg(&discounts_tsv)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .status()
        .expect("failed to execute diagnostics binary");

    assert!(status.success(), "diagnostics command failed");

    let variants = fs::read_to_string(&variants_tsv).expect("failed to read variants output");
    assert!(
        variants.contains("chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth")
    );
    assert!(
        variants.lines().any(|line| line == "chr1\t4\tA\tT\t1\t1"),
        "expected variant row in output, got:\n{}",
        variants
    );

    let inconsistencies =
        fs::read_to_string(&inconsistencies_tsv).expect("failed to read inconsistencies output");
    assert_eq!(
        inconsistencies,
        "read1_position\tread2_position\tdiscordance_type\tcount\n"
    );

    let discounts = fs::read_to_string(&discounts_tsv).expect("failed to read discount output");
    assert!(discounts.contains("mismatch_type\tread_num\tread_position\tdiscount_count"));
    assert!(
        discounts.lines().any(|line| line == "A>T\t1\t4\t1"),
        "expected discount row in output, got:\n{}",
        discounts
    );
}

#[test]
fn integration_diagnostics_can_write_discounts_to_stdout() {
    let temp_dir = unique_temp_dir("diagnostics_stdout_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
    write_test_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let binary = env!("CARGO_BIN_EXE_tasmanian-diagnostics");
    let output = Command::new(binary)
        .arg("-q")
        .arg("0")
        .arg("--min-map-quality")
        .arg("0")
        .arg("--genomic-threshold")
        .arg("1")
        .arg("--genomic-depth-threshold")
        .arg("1")
        .arg("--variants-output")
        .arg(&variants_tsv)
        .arg("--inconsistencies-output")
        .arg(&inconsistencies_tsv)
        .arg("--discount-output")
        .arg("-")
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .output()
        .expect("failed to execute diagnostics binary");

    assert!(output.status.success(), "diagnostics command failed");

    let stdout = String::from_utf8(output.stdout).expect("stdout was not valid utf-8");
    assert!(
        stdout.contains("mismatch_type\tread_num\tread_position\tdiscount_count"),
        "expected discount header in stdout, got:\n{}",
        stdout
    );
    assert!(
        stdout.lines().any(|line| line == "A>T\t1\t4\t1"),
        "expected discount row in stdout, got:\n{}",
        stdout
    );
}