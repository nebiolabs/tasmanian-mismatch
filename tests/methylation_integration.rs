mod test_utils;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::index;
use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
use std::fs;
use std::path::Path;
use std::process::Command;
use test_utils::{log_command, log_line, repo_log_path, unique_temp_dir};

/// Reference: "AAAAGAAA" + "CCCCTCCC" (16 bp, chr1).
///
/// - read_bisulfite: reverse strand, covers [0, 8). Reference has G at 0-based
///   position 4; the read carries A there. In reference-forward orientation
///   this looks like a G>A mismatch, but on the reverse strand it is the
///   bisulfite signature of an unmethylated C (C>T on the original strand),
///   which methylation mode should collapse back to a non-event (G>G) in
///   *reference* orientation.
/// - read_snp: forward strand, covers [8, 16). Reference has T at 0-based
///   position 12; the read carries G there. This is an ordinary mismatch
///   unrelated to bisulfite conversion and must be reported identically
///   with methylation mode on or off.
fn write_test_bam(path: &Path) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 16);
    header.push_record(&sq);

    let header_view = HeaderView::from_header(&header);
    let mut writer =
        Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");

    let bisulfite_line =
        b"read_bisulfite\t16\tchr1\t1\t60\t8M\t*\t0\t0\tAAAAAAAA\tIIIIIIII\tNM:i:1";
    let bisulfite_record =
        Record::from_sam(&header_view, bisulfite_line).expect("failed to parse SAM line");
    writer
        .write(&bisulfite_record)
        .expect("failed to write BAM record");

    let snp_line = b"read_snp\t0\tchr1\t9\t60\t8M\t*\t0\t0\tCCCCGCCC\tIIIIIIII\tNM:i:1";
    let snp_record = Record::from_sam(&header_view, snp_line).expect("failed to parse SAM line");
    writer
        .write(&snp_record)
        .expect("failed to write BAM record");
}

fn run_diagnostics(
    log_path: &Path,
    binary: &str,
    fixture_bam: &Path,
    reference_fa: &Path,
    variants_tsv: &Path,
    inconsistencies_tsv: &Path,
    discounts_tsv: &Path,
    methylation: bool,
) {
    let mut args: Vec<&str> = vec![
        "-q",
        "0",
        "--min-map-quality",
        "0",
        "--genomic-threshold",
        "1",
        "--genomic-depth-threshold",
        "1",
    ];
    if methylation {
        args.push("-m");
    }
    let variants_str = variants_tsv.to_string_lossy().to_string();
    let inconsistencies_str = inconsistencies_tsv.to_string_lossy().to_string();
    let discounts_str = discounts_tsv.to_string_lossy().to_string();
    let bam_str = fixture_bam.to_string_lossy().to_string();
    let reference_str = reference_fa.to_string_lossy().to_string();
    args.extend([
        "--variants-output",
        &variants_str,
        "--inconsistencies-output",
        &inconsistencies_str,
        "--discount-output",
        &discounts_str,
        &bam_str,
        &reference_str,
    ]);

    log_command(log_path, binary, &args);
    let output = Command::new(binary)
        .args(&args)
        .output()
        .expect("failed to execute tasmanian-diagnostics");
    log_line(
        log_path,
        &format!(
            "diagnostics (methylation={}) status: {}",
            methylation, output.status
        ),
    );
    log_line(
        log_path,
        &format!(
            "diagnostics stdout:\n{}",
            String::from_utf8_lossy(&output.stdout)
        ),
    );
    log_line(
        log_path,
        &format!(
            "diagnostics stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        ),
    );
    assert!(
        output.status.success(),
        "tasmanian-diagnostics (methylation={}) failed",
        methylation
    );
}

#[test]
fn integration_methylation_mode_collapses_bisulfite_signature_but_not_real_snp() {
    let temp_dir = unique_temp_dir("methylation_integration");
    let log_path = repo_log_path("methylation_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");

    log_line(&log_path, "Starting methylation-mode integration test");

    fs::write(&reference_fa, ">chr1\nAAAAGAAACCCCTCCC\n").expect("failed to write reference");
    write_test_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let diagnostics_bin = env!("CARGO_BIN_EXE_tasmanian-diagnostics");

    // Baseline run: methylation mode off.
    let off_variants = temp_dir.join("variants_off.tsv");
    run_diagnostics(
        &log_path,
        diagnostics_bin,
        &fixture_bam,
        &reference_fa,
        &off_variants,
        &temp_dir.join("inconsistencies_off.tsv"),
        &temp_dir.join("discounts_off.tsv"),
        false,
    );
    let off_text = fs::read_to_string(&off_variants).expect("failed to read variants_off.tsv");
    log_line(&log_path, &format!("variants_off.tsv:\n{}", off_text));

    assert!(
        off_text
            .lines()
            .any(|line| line == "chr1\t4\tG\tA\t1\t1"),
        "expected uncollapsed G>A row at position 4 without methylation mode, got:\n{}",
        off_text
    );
    assert!(
        off_text
            .lines()
            .any(|line| line == "chr1\t12\tT\tG\t1\t1"),
        "expected T>G SNP row at position 12 without methylation mode, got:\n{}",
        off_text
    );

    // Methylation mode on: the bisulfite-consistent G>A signature must collapse
    // back to a same-base (non-mismatch) row in *reference* orientation, while
    // the unrelated T>G SNP must be reported identically.
    let on_variants = temp_dir.join("variants_on.tsv");
    run_diagnostics(
        &log_path,
        diagnostics_bin,
        &fixture_bam,
        &reference_fa,
        &on_variants,
        &temp_dir.join("inconsistencies_on.tsv"),
        &temp_dir.join("discounts_on.tsv"),
        true,
    );
    let on_text = fs::read_to_string(&on_variants).expect("failed to read variants_on.tsv");
    log_line(&log_path, &format!("variants_on.tsv:\n{}", on_text));

    assert!(
        !on_text.lines().any(|line| line == "chr1\t4\tG\tA\t1\t1"),
        "G>A row at position 4 should have collapsed under methylation mode, got:\n{}",
        on_text
    );
    assert!(
        on_text.lines().any(|line| line == "chr1\t4\tG\tG\t1\t1"),
        "expected collapsed G>G row (reference orientation) at position 4 with methylation mode, got:\n{}",
        on_text
    );
    assert!(
        on_text.lines().any(|line| line == "chr1\t12\tT\tG\t1\t1"),
        "expected T>G SNP row at position 12 to be unaffected by methylation mode, got:\n{}",
        on_text
    );
}
