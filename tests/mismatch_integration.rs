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
fn integration_mismatch_fixture_bam_produces_expected_counts() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixture_bam = manifest_dir.join("tests/test.bam");
    let fixture_bai = manifest_dir.join("tests/test.bam.bai");

    assert!(fixture_bam.exists(), "missing BAM fixture at {}", fixture_bam.display());
    assert!(fixture_bai.exists(), "missing BAM index at {}", fixture_bai.display());

    let temp_dir = unique_temp_dir("mismatch_integration");
    let reference_fa = temp_dir.join("reference.fa");
    let output_tsv = temp_dir.join("mismatch.tsv");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");

    let binary = env!("CARGO_BIN_EXE_tasmanian-mismatch");
    let status = Command::new(binary)
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
        .status()
        .expect("failed to execute mismatch binary");

    assert!(status.success(), "mismatch command failed");

    let output = fs::read_to_string(&output_tsv).expect("failed to read mismatch output");
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