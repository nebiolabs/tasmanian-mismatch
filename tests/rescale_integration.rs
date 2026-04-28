mod test_utils;

use rust_htslib::bam;
use rust_htslib::bam::index;
use rust_htslib::bam::{Format, Header, HeaderView, Read, Record, Writer};
use rust_htslib::bam::header::HeaderRecord;
use std::fs;
use std::io::Read as IoRead;
use std::path::Path;
use std::process::{Command, Stdio};
use test_utils::{log_line, repo_log_path, unique_temp_dir};


fn write_test_bam(path: &Path) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", &"chr1");
    sq.push_tag(b"LN", &8);
    header.push_record(&sq);

    let header_view = HeaderView::from_header(&header);
    let mut writer = Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");

    // Reference is ACGTACGT. This read has one mismatch at read position 4 (A->T).
    let sam_line = b"read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tNM:i:1";
    let record = Record::from_sam(&header_view, sam_line).expect("failed to parse SAM line");
    writer.write(&record).expect("failed to write BAM record");
}

#[test]
fn integration_rescale_matrix_produces_rescaled_sam() {
    let temp_dir = unique_temp_dir("rescale_integration");
    let log_path = repo_log_path("rescale_integration");
    let input_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let matrix_tsv = temp_dir.join("weights.tsv");
    let output_bam = temp_dir.join("rescaled.bam");
    let output_sam = temp_dir.join("rescaled.sam");

    log_line(&log_path, "Starting rescale integration test (file matrix)");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
    fs::write(&matrix_tsv, "1\t4\tA\tT\t0.5\n").expect("failed to write matrix");
    write_test_bam(&input_bam);

    index::build(&input_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let binary = env!("CARGO_BIN_EXE_tasmanian-rescale-quality");
    let output = Command::new(binary)
        .arg(&input_bam)
        .arg(&reference_fa)
        .arg(&matrix_tsv)
        .arg("-o")
        .arg(&output_bam)
        .output()
        .expect("failed to execute rescale binary");
    log_line(&log_path, &format!("Command status: {}", output.status));
    log_line(
        &log_path,
        &format!("Command stdout:\n{}", String::from_utf8_lossy(&output.stdout)),
    );
    log_line(
        &log_path,
        &format!("Command stderr:\n{}", String::from_utf8_lossy(&output.stderr)),
    );

    assert!(output.status.success(), "rescale command failed");
    assert!(output_bam.exists(), "rescaled BAM file was not created");

    let mut bam_reader = bam::Reader::from_path(&output_bam).expect("failed to open rescaled BAM");
    let mut records = bam_reader.records(); // records is Option<Result<Record>> panics are 1.Option 2.Result
    let record = records
        .next()
        .expect("expected one record in rescaled BAM")
        .expect("failed to read rescaled record");

    // Original 'I' Phred is 40. With 0.5 scaling, position 4 should become 20.
    assert_eq!(record.qual()[4], 20, "quality score at mismatch position was not rescaled");

    // Also write a SAM file from the rescaled BAM to validate end-to-end output artifact.
    let out_header = Header::from_template(bam_reader.header());
    let mut sam_writer = Writer::from_path(&output_sam, &out_header, Format::Sam)
        .expect("failed to open SAM writer");

    let mut second_reader = bam::Reader::from_path(&output_bam).expect("failed to re-open rescaled BAM");
    for rec in second_reader.records() {
        sam_writer.write(&rec.expect("failed to read record for SAM conversion")).expect("failed to write SAM record");
    }
    drop(sam_writer);

    assert!(output_sam.exists(), "rescaled SAM file was not created");
    let sam_text = fs::read_to_string(&output_sam).expect("failed to read rescaled SAM");
    let data_line = sam_text
        .lines()
        .find(|line| !line.starts_with('@'))
        .expect("expected a data line in SAM output");
    let fields: Vec<&str> = data_line.split('\t').collect();
    assert!(fields.len() >= 11, "SAM line did not contain required fields");
    assert_eq!(
        fields[10],
        "IIII5III",
        "expected rescaled quality string in SAM QUAL column"
    );
    log_line(&log_path, &format!("SAM output line:\n{}", data_line));
}

#[test]
fn integration_rescale_can_consume_matrix_from_mismatch_stdout() {
    let temp_dir = unique_temp_dir("rescale_pipeline_integration");
    let log_path = repo_log_path("rescale_pipeline_integration");
    let input_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");
    let output_bam = temp_dir.join("rescaled_from_pipe.bam");

    log_line(&log_path, "Starting rescale integration test (stdin matrix)");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
    write_test_bam(&input_bam);
    index::build(&input_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let mismatch_bin = env!("CARGO_BIN_EXE_tasmanian-mismatch");
    let mismatch_output = Command::new(mismatch_bin)
        .arg("-q")
        .arg("0")
        .arg("-m")
        .arg("0")
        .arg("--position-mode")
        .arg("read")
        .arg("--emit-rescaling-matrix")
        .arg(&input_bam)
        .arg(&reference_fa)
        .output()
        .expect("failed to execute mismatch binary");
    log_line(
        &log_path,
        &format!("mismatch status: {}", mismatch_output.status),
    );
    log_line(
        &log_path,
        &format!("mismatch stdout:\n{}", String::from_utf8_lossy(&mismatch_output.stdout)),
    );
    log_line(
        &log_path,
        &format!("mismatch stderr:\n{}", String::from_utf8_lossy(&mismatch_output.stderr)),
    );

    assert!(mismatch_output.status.success(), "mismatch command failed");
    assert!(
        !mismatch_output.stdout.is_empty(),
        "mismatch matrix stdout was empty"
    );

    let rescale_bin = env!("CARGO_BIN_EXE_tasmanian-rescale-quality");
    let mut child = Command::new(rescale_bin)
        .arg(&input_bam)
        .arg(&reference_fa)
        .arg("-")
        .arg("-o")
        .arg(&output_bam)
        .stdin(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to execute rescale binary");

    {
        let stdin = child.stdin.as_mut().expect("failed to open child stdin");
        use std::io::Write;
        stdin
            .write_all(&mismatch_output.stdout)
            .expect("failed to feed matrix stdin");
    }

    let status = child.wait().expect("failed waiting for rescale process");
    let mut rescale_stderr = String::new();
    if let Some(mut stderr) = child.stderr.take() {
        stderr
            .read_to_string(&mut rescale_stderr)
            .expect("failed to read rescale stderr");
    }
    log_line(&log_path, &format!("rescale status: {}", status));
    log_line(&log_path, &format!("rescale stderr:\n{}", rescale_stderr));
    assert!(status.success(), "rescale command failed");
    assert!(output_bam.exists(), "rescaled BAM file was not created");

    let mut bam_reader = bam::Reader::from_path(&output_bam).expect("failed to open rescaled BAM");
    let mut records = bam_reader.records();
    let record = records
        .next()
        .expect("expected one record in rescaled BAM")
        .expect("failed to read rescaled record");

    // Placeholder matrix emission currently uses 1.0 scaling for every key, so qualities stay unchanged.
    assert_eq!(record.qual()[4], 40, "quality should remain unchanged with 1.0 scaling");
}