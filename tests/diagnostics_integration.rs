use std::fs;
use std::path::PathBuf;
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

#[test]
fn integration_diagnostics_fixture_bam_produces_expected_outputs() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixture_bam = manifest_dir.join("tests/test.bam");
    let fixture_bai = manifest_dir.join("tests/test.bam.bai");

    assert!(fixture_bam.exists(), "missing BAM fixture at {}", fixture_bam.display());
    assert!(fixture_bai.exists(), "missing BAM index at {}", fixture_bai.display());

    let temp_dir = unique_temp_dir("diagnostics_integration");
    let reference_fa = temp_dir.join("reference.fa");
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");
    let discounts_tsv = temp_dir.join("discounts.tsv");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");

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
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixture_bam = manifest_dir.join("tests/test.bam");
    let fixture_bai = manifest_dir.join("tests/test.bam.bai");

    assert!(fixture_bam.exists(), "missing BAM fixture at {}", fixture_bam.display());
    assert!(fixture_bai.exists(), "missing BAM index at {}", fixture_bai.display());

    let temp_dir = unique_temp_dir("diagnostics_stdout_integration");
    let reference_fa = temp_dir.join("reference.fa");
    let variants_tsv = temp_dir.join("variants.tsv");
    let inconsistencies_tsv = temp_dir.join("inconsistencies.tsv");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");

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