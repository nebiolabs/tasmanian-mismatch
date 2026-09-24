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

#[test]
fn integration_mismatch_fixture_bam_produces_expected_counts() {
    let temp_dir = unique_temp_dir("mismatch_integration");
    let log_path = repo_log_path("mismatch_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let output_tsv = temp_dir.join("mismatch.tsv");

    log_line(&log_path, "Starting mismatch integration test");

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

    let binary = env!("CARGO_BIN_EXE_tasmanian-mismatch");
    log_command(
        &log_path,
        binary,
        &[
            "-q",
            "0",
            "-m",
            "0",
            "--position-mode",
            "read",
            "-o",
            &output_tsv.to_string_lossy(),
            &fixture_bam.to_string_lossy(),
            &reference_fa.to_string_lossy(),
        ],
    );
    let output = Command::new(binary)
        .arg("-q")
        .arg("0")
        .arg("-m")
        .arg("0")
        .arg("--position-mode")
        .arg("read")
        .arg("-o")
        .arg(&output_tsv)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .output()
        .expect("failed to execute mismatch binary");
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

    assert!(output.status.success(), "mismatch command failed");

    let output = fs::read_to_string(&output_tsv).expect("failed to read mismatch output");
    log_line(&log_path, &format!("Mismatch output file:\n{}", output));
    let mut lines = output.lines();
    assert_eq!(
        lines.next(),
        Some("base_change\tread_num\treference_order\tread_position\tcount")
    );
    assert!(
        output.lines().any(|line| line == "A>T\t1\t1\t5\t1"),
        "expected A>T mismatch row in output, got:\n{}",
        output
    );
}

#[test]
fn integration_min_max_position_excludes_bases_outside_range() {
    let temp_dir = unique_temp_dir("mismatch_position_range_integration");
    let log_path = repo_log_path("mismatch_position_range_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let output_tsv = temp_dir.join("mismatch.tsv");

    log_line(&log_path, "Starting min/max position integration test");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
    write_test_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let binary = env!("CARGO_BIN_EXE_tasmanian-mismatch");
    let args = [
        "-q",
        "0",
        "-m",
        "0",
        "--position-mode",
        "read",
        "--min-read-position",
        "6",
        "--max-read-position",
        "8",
        "-o",
        &output_tsv.to_string_lossy(),
        &fixture_bam.to_string_lossy(),
        &reference_fa.to_string_lossy(),
    ];
    log_command(&log_path, binary, &args);
    let output = Command::new(binary)
        .args(args)
        .output()
        .expect("failed to execute mismatch binary");
    log_line(&log_path, &format!("Command status: {}", output.status));
    log_line(
        &log_path,
        &format!(
            "Command stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        ),
    );

    assert!(output.status.success(), "mismatch command failed");

    let output = fs::read_to_string(&output_tsv).expect("failed to read mismatch output");
    log_line(&log_path, &format!("Mismatch output file:\n{}", output));

    // The A>T mismatch at read_position 5 falls outside [6, 8] and must be excluded.
    assert!(
        !output.contains("A>T"),
        "A>T mismatch at position 5 should be filtered out, got:\n{}",
        output
    );
    // Every remaining row's read_position (4th column) must fall within [6, 8].
    for line in output.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        let read_position: usize = fields[3].parse().expect("read_position should be numeric");
        assert!(
            (6..=8).contains(&read_position),
            "row outside requested range in output:\n{}",
            output
        );
    }
}
