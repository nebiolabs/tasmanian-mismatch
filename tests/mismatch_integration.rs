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

fn write_two_chrom_bam(path: &Path) {
    let mut header = Header::new();
    let mut sq1 = HeaderRecord::new(b"SQ");
    sq1.push_tag(b"SN", "chr1");
    sq1.push_tag(b"LN", 8);
    header.push_record(&sq1);
    let mut sq2 = HeaderRecord::new(b"SQ");
    sq2.push_tag(b"SN", "chr2");
    sq2.push_tag(b"LN", 8);
    header.push_record(&sq2);

    let header_view = HeaderView::from_header(&header);
    let mut writer =
        Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");

    // chr1 reference ACGTACGT; read has one mismatch (A->T) at position 5.
    let read1 = Record::from_sam(
        &header_view,
        b"read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:1",
    )
    .expect("failed to parse SAM line");
    writer.write(&read1).expect("failed to write BAM record");

    // chr2 reference ACGTACGT; read has one mismatch (T->A) at position 8.
    let read2 = Record::from_sam(
        &header_view,
        b"read2\t0\tchr2\t1\t60\t8M\t*\t0\t0\tACGTACGA\tIIIIIIII\tNM:i:1",
    )
    .expect("failed to parse SAM line");
    writer.write(&read2).expect("failed to write BAM record");
}

#[test]
fn integration_bed_filter_mode_include_keeps_only_overlapping_reads() {
    let temp_dir = unique_temp_dir("mismatch_bed_include_integration");
    let log_path = repo_log_path("mismatch_bed_include_integration");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let bed_file = temp_dir.join("regions.bed");

    log_line(
        &log_path,
        "Starting --bed-filter-mode include integration test",
    );

    fs::write(&reference_fa, ">chr1\nACGTACGT\n>chr2\nACGTACGT\n")
        .expect("failed to write reference");
    write_two_chrom_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");
    // Covers all of chr1 (read1's mismatch A>T), none of chr2 (read2's mismatch T>A).
    fs::write(&bed_file, "chr1\t0\t8\n").expect("failed to write BED file");

    let binary = env!("CARGO_BIN_EXE_tasmanian-mismatch");

    let run = |mode: &str, output_name: &str| -> String {
        let output_tsv = temp_dir.join(output_name);
        let args = [
            "-q",
            "0",
            "-m",
            "0",
            "--position-mode",
            "read",
            "-b",
            &bed_file.to_string_lossy(),
            "--bed-filter-mode",
            mode,
            "-o",
            &output_tsv.to_string_lossy(),
            &fixture_bam.to_string_lossy(),
            &reference_fa.to_string_lossy(),
        ];
        log_command(&log_path, binary, &args);
        let output = Command::new(binary)
            .args(args)
            .output()
            .unwrap_or_else(|_| panic!("failed to execute mismatch binary in {mode} mode"));
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
            "mismatch command failed in {mode} mode"
        );
        fs::read_to_string(&output_tsv).expect("failed to read mismatch output")
    };

    let include_output = run("include", "include.tsv");
    log_line(&log_path, &format!("include output:\n{include_output}"));
    assert!(
        include_output.contains("A>T"),
        "include mode should keep chr1's A>T mismatch, got:\n{include_output}"
    );
    assert!(
        !include_output.contains("T>A"),
        "include mode should drop chr2's T>A mismatch (chr2 isn't in the BED), got:\n{include_output}"
    );

    // filter mode is the exact logical inverse of include on the same BAM+BED.
    let filter_output = run("filter", "filter.tsv");
    log_line(&log_path, &format!("filter output:\n{filter_output}"));
    assert!(
        !filter_output.contains("A>T"),
        "filter mode should drop chr1's A>T mismatch (chr1 is in the BED), got:\n{filter_output}"
    );
    assert!(
        filter_output.contains("T>A"),
        "filter mode should keep chr2's T>A mismatch, got:\n{filter_output}"
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

fn run_mismatch(temp_dir: &Path, log_path: &Path, output_name: &str, extra: &[&str]) -> String {
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let output_tsv = temp_dir.join(output_name);
    let binary = env!("CARGO_BIN_EXE_tasmanian-mismatch");
    let mut args = vec![
        "-q".to_string(),
        "0".to_string(),
        "-m".to_string(),
        "0".to_string(),
        "--position-mode".to_string(),
        "read".to_string(),
    ];
    args.extend(extra.iter().map(|a| a.to_string()));
    args.extend([
        "-o".to_string(),
        output_tsv.to_string_lossy().to_string(),
        fixture_bam.to_string_lossy().to_string(),
        reference_fa.to_string_lossy().to_string(),
    ]);
    let arg_refs: Vec<&str> = args.iter().map(String::as_str).collect();
    log_command(log_path, binary, &arg_refs);
    let output = Command::new(binary)
        .args(&args)
        .output()
        .expect("failed to execute mismatch binary");
    log_line(
        log_path,
        &format!("stderr:\n{}", String::from_utf8_lossy(&output.stderr)),
    );
    assert!(
        output.status.success(),
        "mismatch command failed: {extra:?}"
    );
    fs::read_to_string(&output_tsv).expect("failed to read mismatch output")
}

fn prepare_two_chrom_inputs(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
    let temp_dir = unique_temp_dir(name);
    let log_path = repo_log_path(name);
    fs::write(
        temp_dir.join("reference.fa"),
        ">chr1\nACGTACGT\n>chr2\nACGTACGT\n",
    )
    .expect("failed to write reference");
    let fixture_bam = temp_dir.join("input.bam");
    write_two_chrom_bam(&fixture_bam);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");
    (temp_dir, log_path)
}

#[test]
fn integration_window_size_reports_each_window_separately() {
    let (temp_dir, log_path) = prepare_two_chrom_inputs("mismatch_window_integration");
    // 4 bp windows: each 8 bp read starts at 0, so it lands in the first window of its contig.
    let output = run_mismatch(&temp_dir, &log_path, "windows.tsv", &["--window-size", "4"]);
    log_line(&log_path, &format!("window output:\n{output}"));

    let mut lines = output.lines();
    assert_eq!(
        lines.next(),
        Some(
            "chrom\tstart\tend\tread_num\tbase_change\tcount\tref_total\trate\tall_bases\trate_all"
        )
    );
    // Each read covers two A and two T reference bases.
    assert!(
        output
            .lines()
            .any(|l| l == "chr1\t0\t4\t1\tA>T\t1\t2\t0.500000\t8\t0.125000"),
        "expected chr1's A>T in its first window, got:\n{output}"
    );
    assert!(
        output
            .lines()
            .any(|l| l == "chr2\t0\t4\t1\tT>A\t1\t2\t0.500000\t8\t0.125000"),
        "expected chr2's T>A in its first window, got:\n{output}"
    );
    // Windows are written in contig order, and windows without reads are omitted.
    let mut windows: Vec<String> = output
        .lines()
        .skip(1)
        .map(|l| l.splitn(4, '\t').take(3).collect::<Vec<_>>().join("\t"))
        .collect();
    windows.dedup();
    assert_eq!(windows, ["chr1\t0\t4", "chr2\t0\t4"]);
}

#[test]
fn integration_window_output_does_not_depend_on_region_size() {
    let (temp_dir, log_path) = prepare_two_chrom_inputs("mismatch_window_region_integration");
    let default_chunks = run_mismatch(&temp_dir, &log_path, "default.tsv", &["--window-size", "4"]);
    // A 3 bp chunk rounds up to one 4 bp window, so each window is its own chunk.
    let small_chunks = run_mismatch(
        &temp_dir,
        &log_path,
        "small.tsv",
        &["--window-size", "4", "-r", "3"],
    );
    assert_eq!(default_chunks, small_chunks);
}

#[test]
fn integration_bootstrap_adds_reproducible_interval_columns() {
    let (temp_dir, log_path) = prepare_two_chrom_inputs("mismatch_bootstrap_integration");
    let plain = run_mismatch(&temp_dir, &log_path, "plain.tsv", &[]);
    let boot_args = ["-r", "4", "--bootstrap", "50"];
    let boot = run_mismatch(&temp_dir, &log_path, "boot.tsv", &boot_args);
    log_line(&log_path, &format!("bootstrap output:\n{boot}"));

    assert_eq!(
        boot.lines().next(),
        Some(
            "base_change\tread_num\treference_order\tread_position\tcount\tfrequency\tci_low\tci_high"
        )
    );
    // The leading columns are the plain table's, row for row.
    let leading = |table: &str| -> Vec<String> {
        table
            .lines()
            .skip(1)
            .map(|l| l.splitn(6, '\t').take(5).collect::<Vec<_>>().join("\t"))
            .collect()
    };
    assert_eq!(leading(&boot), leading(&plain));

    for line in boot.lines().skip(1) {
        let fields: Vec<f64> = line
            .split('\t')
            .skip(5)
            .map(|v| v.parse().expect("frequency/CI columns should be numeric"))
            .collect();
        let [freq, low, high] = fields[..] else {
            panic!("expected frequency, ci_low, ci_high in: {line}");
        };
        assert!(
            (0.0..=1.0).contains(&low) && low <= high && high <= 1.0,
            "bad interval in: {line}"
        );
        assert!((0.0..=1.0).contains(&freq), "bad frequency in: {line}");
    }

    let again = run_mismatch(&temp_dir, &log_path, "boot_again.tsv", &boot_args);
    assert_eq!(boot, again, "same seed should reproduce the same intervals");
}
