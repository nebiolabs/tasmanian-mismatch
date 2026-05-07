mod test_utils;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::index;
use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
use std::fs;
use std::io::Read;
use std::path::Path;
use std::process::{Command, Stdio};
use test_utils::{log_command, log_line, repo_log_path, unique_temp_dir};

fn write_test_bam(path: &Path) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", &"chr1");
    sq.push_tag(b"LN", &8);
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
fn integration_all_binaries_with_options_and_three_way_pipe() {
    let temp_dir = unique_temp_dir("all_binaries_and_pipes");
    let log_path = repo_log_path("all_binaries_and_pipes");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let bed_file = temp_dir.join("regions.bed");
    let discount_file = temp_dir.join("discounts_in.tsv");
    let matrix_file = temp_dir.join("matrix.tsv");

    log_line(&log_path, "Starting all-binaries integration test");

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

    // Non-overlapping BED interval for this tiny fixture read.
    fs::write(&bed_file, "chr1\t100\t101\n").expect("failed to write bed file");
    fs::write(
        &discount_file,
        "base_change\tread_num\tread_position\tdiscount_count\nA>T\t1\t5\t1\n",
    )
    .expect("failed to write discount file");
    fs::write(&matrix_file, "1\t4\tA\tT\t0.5\n").expect("failed to write matrix file");

    // 1) tasmanian-mismatch with broad option coverage.
    let mismatch_output = temp_dir.join("mismatch_normalized.tsv");
    let mismatch_bin = env!("CARGO_BIN_EXE_tasmanian-mismatch");
    log_command(
        &log_path,
        mismatch_bin,
        &[
            "-t",
            "1",
            "-r",
            "1000",
            "-q",
            "0",
            "-m",
            "0",
            "-f",
            "0",
            "-F",
            "0",
            "-G",
            "0",
            "--overlap-mode",
            "stretch",
            "--position-mode",
            "insert",
            "--min-fragment-length",
            "0",
            "--max-fragment-length",
            "1000",
            "--methylation-mode",
            "--cpg-only",
            "-b",
            &bed_file.to_string_lossy(),
            "--bed-filter-mode",
            "mask",
            "--discount-table",
            &discount_file.to_string_lossy(),
            "--normalize",
            "-o",
            &mismatch_output.to_string_lossy(),
            &fixture_bam.to_string_lossy(),
            &reference_fa.to_string_lossy(),
        ],
    );
    let mismatch_run = Command::new(mismatch_bin)
        .arg("-t")
        .arg("1")
        .arg("-r")
        .arg("1000")
        .arg("-q")
        .arg("0")
        .arg("-m")
        .arg("0")
        .arg("-f")
        .arg("0")
        .arg("-F")
        .arg("0")
        .arg("-G")
        .arg("0")
        .arg("--overlap-mode")
        .arg("stretch")
        .arg("--position-mode")
        .arg("insert")
        .arg("--min-fragment-length")
        .arg("0")
        .arg("--max-fragment-length")
        .arg("1000")
        .arg("--methylation-mode")
        .arg("--cpg-only")
        .arg("-b")
        .arg(&bed_file)
        .arg("--bed-filter-mode")
        .arg("mask")
        .arg("--discount-table")
        .arg(&discount_file)
        .arg("--normalize")
        .arg("-o")
        .arg(&mismatch_output)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .output()
        .expect("failed to execute tasmanian-mismatch");
    log_line(
        &log_path,
        &format!("mismatch status: {}", mismatch_run.status),
    );
    log_line(
        &log_path,
        &format!(
            "mismatch stdout:\n{}",
            String::from_utf8_lossy(&mismatch_run.stdout)
        ),
    );
    log_line(
        &log_path,
        &format!(
            "mismatch stderr:\n{}",
            String::from_utf8_lossy(&mismatch_run.stderr)
        ),
    );
    assert!(
        mismatch_run.status.success(),
        "tasmanian-mismatch command failed"
    );
    let mismatch_text =
        fs::read_to_string(&mismatch_output).expect("failed to read mismatch output");
    log_line(
        &log_path,
        &format!("mismatch_normalized.tsv:\n{}", mismatch_text),
    );
    assert!(
        mismatch_text
            .lines()
            .next()
            .unwrap_or("")
            .contains("normalized_frequency"),
        "expected normalized mismatch header, got:\n{}",
        mismatch_text
    );

    // 2) tasmanian-diagnostics with broad option coverage.
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");
    let discounts_tsv = temp_dir.join("discounts.tsv");
    let diagnostics_bin = env!("CARGO_BIN_EXE_tasmanian-diagnostics");
    log_command(
        &log_path,
        diagnostics_bin,
        &[
            "-t",
            "1",
            "-r",
            "1000",
            "--softclip-threshold",
            "0.5",
            "-q",
            "0",
            "--min-map-quality",
            "0",
            "-m",
            "--cpg-only",
            "--use-read-len-max",
            "8",
            "--use-insert-mode",
            "--genomic-threshold",
            "1",
            "--genomic-depth-threshold",
            "1",
            "-f",
            "0",
            "-F",
            "0",
            "-G",
            "0",
            "-b",
            &bed_file.to_string_lossy(),
            "--bed-filter-mode",
            "filter",
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
    let diagnostics_run = Command::new(diagnostics_bin)
        .arg("-t")
        .arg("1")
        .arg("-r")
        .arg("1000")
        .arg("--softclip-threshold")
        .arg("0.5")
        .arg("-q")
        .arg("0")
        .arg("--min-map-quality")
        .arg("0")
        .arg("-m")
        .arg("--cpg-only")
        .arg("--use-read-len-max")
        .arg("8")
        .arg("--use-insert-mode")
        .arg("--genomic-threshold")
        .arg("1")
        .arg("--genomic-depth-threshold")
        .arg("1")
        .arg("-f")
        .arg("0")
        .arg("-F")
        .arg("0")
        .arg("-G")
        .arg("0")
        .arg("-b")
        .arg(&bed_file)
        .arg("--bed-filter-mode")
        .arg("filter")
        .arg("--variants-output")
        .arg(&variants_tsv)
        .arg("--inconsistencies-output")
        .arg(&inconsistencies_tsv)
        .arg("--discount-output")
        .arg(&discounts_tsv)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .output()
        .expect("failed to execute tasmanian-diagnostics");
    log_line(
        &log_path,
        &format!("diagnostics status: {}", diagnostics_run.status),
    );
    log_line(
        &log_path,
        &format!(
            "diagnostics stdout:\n{}",
            String::from_utf8_lossy(&diagnostics_run.stdout)
        ),
    );
    log_line(
        &log_path,
        &format!(
            "diagnostics stderr:\n{}",
            String::from_utf8_lossy(&diagnostics_run.stderr)
        ),
    );
    assert!(
        diagnostics_run.status.success(),
        "tasmanian-diagnostics command failed"
    );
    assert!(variants_tsv.exists(), "missing variants output");
    assert!(
        inconsistencies_tsv.exists(),
        "missing inconsistencies output"
    );
    assert!(discounts_tsv.exists(), "missing discount output");

    // 3) tasmanian-rescale-quality with explicit options.
    let rescaled_bam = temp_dir.join("rescaled_from_file_matrix.bam");
    let rescale_bin = env!("CARGO_BIN_EXE_tasmanian-rescale-quality");
    log_command(
        &log_path,
        rescale_bin,
        &[
            &fixture_bam.to_string_lossy(),
            &reference_fa.to_string_lossy(),
            &matrix_file.to_string_lossy(),
            "-t",
            "1",
            "-r",
            "1000",
            "-o",
            &rescaled_bam.to_string_lossy(),
        ],
    );
    let rescale_run = Command::new(rescale_bin)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .arg(&matrix_file)
        .arg("-t")
        .arg("1")
        .arg("-r")
        .arg("1000")
        .arg("-o")
        .arg(&rescaled_bam)
        .output()
        .expect("failed to execute tasmanian-rescale-quality");
    log_line(
        &log_path,
        &format!("rescale(file) status: {}", rescale_run.status),
    );
    log_line(
        &log_path,
        &format!(
            "rescale(file) stdout:\n{}",
            String::from_utf8_lossy(&rescale_run.stdout)
        ),
    );
    log_line(
        &log_path,
        &format!(
            "rescale(file) stderr:\n{}",
            String::from_utf8_lossy(&rescale_run.stderr)
        ),
    );
    assert!(
        rescale_run.status.success(),
        "tasmanian-rescale-quality command failed"
    );
    assert!(rescaled_bam.exists(), "missing rescaled BAM output");

    // 4) Full 3-binary pipe:
    // diagnostics (--discount-output -) -> mismatch (--discount-table -, --emit-rescaling-matrix)
    // -> rescale-quality (matrix_file = -)
    let piped_variants = temp_dir.join("piped_variants.tsv");
    let piped_incons = temp_dir.join("piped_incons.tsv");
    let piped_rescaled_bam = temp_dir.join("rescaled_from_pipe.bam");

    let mut diagnostics_child = Command::new(diagnostics_bin)
        .arg("-q")
        .arg("0")
        .arg("--min-map-quality")
        .arg("0")
        .arg("--genomic-threshold")
        .arg("1")
        .arg("--genomic-depth-threshold")
        .arg("1")
        .arg("--variants-output")
        .arg(&piped_variants)
        .arg("--inconsistencies-output")
        .arg(&piped_incons)
        .arg("--discount-output")
        .arg("-")
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn diagnostics child");

    let diagnostics_stdout = diagnostics_child
        .stdout
        .take()
        .expect("failed to capture diagnostics stdout");

    let mut mismatch_child = Command::new(mismatch_bin)
        .arg("-q")
        .arg("0")
        .arg("-m")
        .arg("0")
        .arg("--position-mode")
        .arg("read")
        .arg("--discount-table")
        .arg("-")
        .arg("--emit-rescaling-matrix")
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .stdin(Stdio::from(diagnostics_stdout))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn mismatch child");

    let mismatch_stdout = mismatch_child
        .stdout
        .take()
        .expect("failed to capture mismatch stdout");

    let mut rescale_child = Command::new(rescale_bin)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .arg("-")
        .arg("-t")
        .arg("1")
        .arg("-r")
        .arg("1000")
        .arg("-o")
        .arg(&piped_rescaled_bam)
        .stdin(Stdio::from(mismatch_stdout))
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn rescale child");

    let rescale_pipe_status = rescale_child
        .wait()
        .expect("failed waiting for rescale child");
    let mismatch_pipe_status = mismatch_child
        .wait()
        .expect("failed waiting for mismatch child");
    let diagnostics_pipe_status = diagnostics_child
        .wait()
        .expect("failed waiting for diagnostics child");

    let mut diagnostics_stderr = String::new();
    if let Some(mut stderr) = diagnostics_child.stderr.take() {
        stderr
            .read_to_string(&mut diagnostics_stderr)
            .expect("failed to read diagnostics stderr");
    }
    let mut mismatch_stderr = String::new();
    if let Some(mut stderr) = mismatch_child.stderr.take() {
        stderr
            .read_to_string(&mut mismatch_stderr)
            .expect("failed to read mismatch stderr");
    }
    let mut rescale_stderr = String::new();
    if let Some(mut stderr) = rescale_child.stderr.take() {
        stderr
            .read_to_string(&mut rescale_stderr)
            .expect("failed to read rescale stderr");
    }

    log_line(
        &log_path,
        &format!("pipe diagnostics status: {}", diagnostics_pipe_status),
    );
    log_line(
        &log_path,
        &format!("pipe diagnostics stderr:\n{}", diagnostics_stderr),
    );
    log_line(
        &log_path,
        &format!("pipe mismatch status: {}", mismatch_pipe_status),
    );
    log_line(
        &log_path,
        &format!("pipe mismatch stderr:\n{}", mismatch_stderr),
    );
    log_line(
        &log_path,
        &format!("pipe rescale status: {}", rescale_pipe_status),
    );
    log_line(
        &log_path,
        &format!("pipe rescale stderr:\n{}", rescale_stderr),
    );

    assert!(
        diagnostics_pipe_status.success(),
        "diagnostics in pipe failed"
    );
    assert!(mismatch_pipe_status.success(), "mismatch in pipe failed");
    assert!(rescale_pipe_status.success(), "rescale in pipe failed");
    assert!(
        piped_rescaled_bam.exists(),
        "missing piped rescaled BAM output"
    );
}
