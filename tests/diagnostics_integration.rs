mod test_utils;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::index;
use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
use std::fs;
use std::path::Path;
use std::process::Command;
use test_utils::{log_command, log_line, repo_log_path, unique_temp_dir};

fn write_test_bam(path: &Path) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 8);
    header.push_record(&sq);

    let header_view = HeaderView::from_header(&header);
    let mut writer =
        Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");

    // Reference is ACGTACGT. This read has one mismatch (A->T).
    let sam_line = b"read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:1";
    let record = Record::from_sam(&header_view, sam_line).expect("failed to parse SAM line");
    writer.write(&record).expect("failed to write BAM record");
}

fn write_paired_test_bam(path: &Path) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 32);
    header.push_record(&sq);

    let header_view = HeaderView::from_header(&header);
    let mut writer =
        Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");

    let read1 = Record::from_sam(
        &header_view,
        b"pair1\t99\tchr1\t1\t60\t8M\t=\t5\t12\tACGTACGT\tIIIIIIII\tMC:Z:8M",
    )
    .expect("failed to parse read1 SAM line");
    let read2 = Record::from_sam(
        &header_view,
        b"pair1\t147\tchr1\t5\t60\t8M\t=\t1\t-12\tACGTACGT\tIIIIIIII\tMC:Z:8M",
    )
    .expect("failed to parse read2 SAM line");

    writer.write(&read1).expect("failed to write read1");
    writer.write(&read2).expect("failed to write read2");
}

#[test]
fn integration_diagnostics_fixture_bam_produces_expected_outputs() {
    let temp_dir = unique_temp_dir("diagnostics_integration");
    let log_path = repo_log_path("diagnostics_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");
    let discounts_tsv = temp_dir.join("discounts.tsv");

    log_line(
        &log_path,
        "Starting diagnostics integration test (file output)",
    );

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
    write_test_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");
    log_line(
        &log_path,
        &format!(
            "Prepared inputs: bam={}, reference={}",
            fixture_bam.display(),
            reference_fa.display()
        ),
    );

    let binary = env!("CARGO_BIN_EXE_tasmanian-diagnostics");
    log_command(
        &log_path,
        binary,
        &[
            "-q",
            "0",
            "--min-map-quality",
            "0",
            "--genomic-threshold",
            "1",
            "--genomic-depth-threshold",
            "1",
            "--variants-output",
            &variants_tsv.to_string_lossy(),
            "--inconsistencies-output",
            &inconsistencies_tsv.to_string_lossy(),
            "--discount-output",
            &discounts_tsv.to_string_lossy(),
            &fixture_bam.to_string_lossy(),
            &reference_fa.to_string_lossy(),
        ],
    );
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
        .arg(&discounts_tsv)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .output()
        .expect("failed to execute diagnostics binary");
    log_line(&log_path, &format!("Command status: {}", output.status));
    log_line(
        &log_path,
        &format!(
            "Command stdout:\n{}",
            String::from_utf8_lossy(&output.stdout)
        ),
    );
    log_line(
        &log_path,
        &format!(
            "Command stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        ),
    );

    assert!(output.status.success(), "diagnostics command failed");

    let variants = fs::read_to_string(&variants_tsv).expect("failed to read variants output");
    log_line(&log_path, &format!("variants.tsv:\n{}", variants));
    assert!(variants.contains("chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth"));
    assert!(
        variants.lines().any(|line| line == "chr1\t4\tA\tT\t1\t1"),
        "expected variant row in output, got:\n{}",
        variants
    );

    let inconsistencies =
        fs::read_to_string(&inconsistencies_tsv).expect("failed to read inconsistencies output");
    log_line(
        &log_path,
        &format!("inconsistencies.tsv:\n{}", inconsistencies),
    );
    assert_eq!(
        inconsistencies,
        "read1_position\tread2_position\tdiscordance_type\tcount\n"
    );

    let discounts = fs::read_to_string(&discounts_tsv).expect("failed to read discount output");
    log_line(&log_path, &format!("discounts.tsv:\n{}", discounts));
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
    let log_path = repo_log_path("diagnostics_stdout_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");

    log_line(
        &log_path,
        "Starting diagnostics integration test (stdout discounts)",
    );

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
    write_test_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let binary = env!("CARGO_BIN_EXE_tasmanian-diagnostics");
    log_command(
        &log_path,
        binary,
        &[
            "-q",
            "0",
            "--min-map-quality",
            "0",
            "--genomic-threshold",
            "1",
            "--genomic-depth-threshold",
            "1",
            "--variants-output",
            &variants_tsv.to_string_lossy(),
            "--inconsistencies-output",
            &inconsistencies_tsv.to_string_lossy(),
            "--discount-output",
            "-",
            &fixture_bam.to_string_lossy(),
            &reference_fa.to_string_lossy(),
        ],
    );
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
    log_line(&log_path, &format!("Command status: {}", output.status));
    log_line(
        &log_path,
        &format!(
            "Command stdout:\n{}",
            String::from_utf8_lossy(&output.stdout)
        ),
    );
    log_line(
        &log_path,
        &format!(
            "Command stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        ),
    );

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

#[test]
fn integration_diagnostics_processes_paired_reads() {
    let temp_dir = unique_temp_dir("diagnostics_paired_integration");
    let log_path = repo_log_path("diagnostics_paired_integration");
    let fixture_bam = temp_dir.join("paired_input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");
    let discounts_tsv = temp_dir.join("discounts.tsv");

    log_line(
        &log_path,
        "Starting diagnostics integration test (paired read branch)",
    );

    fs::write(&reference_fa, ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
        .expect("failed to write reference");
    write_paired_test_bam(&fixture_bam);
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
        .arg(&discounts_tsv)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .output()
        .expect("failed to execute diagnostics binary");

    assert!(output.status.success(), "diagnostics command failed");

    let variants = fs::read_to_string(&variants_tsv).expect("failed to read variants output");
    assert!(variants.contains("chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth"));

    let inconsistencies =
        fs::read_to_string(&inconsistencies_tsv).expect("failed to read inconsistencies output");
    assert!(
        inconsistencies.starts_with("read1_position\tread2_position\tdiscordance_type\tcount\n")
    );

    let discounts = fs::read_to_string(&discounts_tsv).expect("failed to read discounts output");
    assert!(discounts.contains("mismatch_type\tread_num\tread_position\tdiscount_count"));
}

fn write_two_chrom_bam(path: &Path) {
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

#[test]
fn integration_diagnostics_bed_filter_mode_include_and_filter_are_inverses() {
    let temp_dir = unique_temp_dir("diagnostics_bed_include_integration");
    let log_path = repo_log_path("diagnostics_bed_include_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let bed_file = temp_dir.join("regions.bed");

    log_line(
        &log_path,
        "Starting diagnostics --bed-filter-mode include/filter integration test",
    );

    fs::write(&reference_fa, ">chr1\nACGTACGT\n>chr2\nACGTACGT\n")
        .expect("failed to write reference");
    write_two_chrom_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");
    // Covers all of chr1 and none of chr2, so every chr2 chunk has no BED intervals.
    fs::write(&bed_file, "chr1\t0\t8\n").expect("failed to write BED file");

    let binary = env!("CARGO_BIN_EXE_tasmanian-diagnostics");

    let run = |mode: &str| -> String {
        let variants_tsv = temp_dir.join(format!("{mode}_variants.tsv"));
        let inconsistencies_tsv = temp_dir.join(format!("{mode}_inconsistencies.tsv"));
        let discounts_tsv = temp_dir.join(format!("{mode}_discounts.tsv"));
        let args = [
            "-q",
            "0",
            "--min-map-quality",
            "0",
            "--genomic-threshold",
            "1",
            "--genomic-depth-threshold",
            "1",
            "-b",
            &bed_file.to_string_lossy(),
            "--bed-filter-mode",
            mode,
            "--variants-output",
            &variants_tsv.to_string_lossy(),
            "--inconsistencies-output",
            &inconsistencies_tsv.to_string_lossy(),
            "--discount-output",
            &discounts_tsv.to_string_lossy(),
            &fixture_bam.to_string_lossy(),
            &reference_fa.to_string_lossy(),
        ];
        log_command(&log_path, binary, &args);
        let output = Command::new(binary)
            .args(args)
            .output()
            .unwrap_or_else(|_| panic!("failed to execute diagnostics binary in {mode} mode"));
        log_line(&log_path, &format!("[{mode}] status: {}", output.status));
        log_line(
            &log_path,
            &format!(
                "[{mode}] stderr:\n{}",
                String::from_utf8_lossy(&output.stderr)
            ),
        );
        assert!(
            output.status.success(),
            "diagnostics command failed in {mode} mode"
        );
        fs::read_to_string(&variants_tsv).expect("failed to read variants output")
    };

    let include_variants = run("include");
    log_line(&log_path, &format!("include variants:\n{include_variants}"));
    assert!(
        include_variants.lines().any(|l| l.starts_with("chr1\t")),
        "include mode should keep chr1's variant, got:\n{include_variants}"
    );
    assert!(
        !include_variants.lines().any(|l| l.starts_with("chr2\t")),
        "include mode should drop chr2's variant (chr2 isn't in the BED), got:\n{include_variants}"
    );

    let filter_variants = run("filter");
    log_line(&log_path, &format!("filter variants:\n{filter_variants}"));
    assert!(
        !filter_variants.lines().any(|l| l.starts_with("chr1\t")),
        "filter mode should drop chr1's variant (chr1 is in the BED), got:\n{filter_variants}"
    );
    assert!(
        filter_variants.lines().any(|l| l.starts_with("chr2\t")),
        "filter mode should keep chr2's variant, got:\n{filter_variants}"
    );
}
