use std::fs;
use std::path::PathBuf;
use std::process::{Command, Stdio};
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
fn integration_all_binaries_with_options_and_three_way_pipe() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixture_bam = manifest_dir.join("tests/test.bam");
    let fixture_bai = manifest_dir.join("tests/test.bam.bai");
    assert!(fixture_bam.exists(), "missing BAM fixture at {}", fixture_bam.display());
    assert!(fixture_bai.exists(), "missing BAM index at {}", fixture_bai.display());

    let temp_dir = PathBuf::from("./test_output");
    let _ = fs::remove_dir_all(&temp_dir);  // Clean previous run
    fs::create_dir_all(&temp_dir).unwrap();
    let reference_fa = temp_dir.join("reference.fa");
    let bed_file = temp_dir.join("regions.bed");
    let discount_file = temp_dir.join("discounts_in.tsv");
    let matrix_file = temp_dir.join("matrix.tsv");

    fs::write(&reference_fa, ">chr1\nACGTACGT\n").expect("failed to write reference");
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
    let mismatch_status = Command::new(mismatch_bin)
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
        .status()
        .expect("failed to execute tasmanian-mismatch");
    assert!(mismatch_status.success(), "tasmanian-mismatch command failed");
    let mismatch_text = fs::read_to_string(&mismatch_output).expect("failed to read mismatch output");
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
    let diagnostics_status = Command::new(diagnostics_bin)
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
        .status()
        .expect("failed to execute tasmanian-diagnostics");
    assert!(diagnostics_status.success(), "tasmanian-diagnostics command failed");
    assert!(variants_tsv.exists(), "missing variants output");
    assert!(inconsistencies_tsv.exists(), "missing inconsistencies output");
    assert!(discounts_tsv.exists(), "missing discount output");

    // 3) tasmanian-rescale-quality with explicit options.
    let rescaled_bam = temp_dir.join("rescaled_from_file_matrix.bam");
    let rescale_bin = env!("CARGO_BIN_EXE_tasmanian-rescale-quality");
    let rescale_status = Command::new(rescale_bin)
        .arg(&fixture_bam)
        .arg(&reference_fa)
        .arg(&matrix_file)
        .arg("-t")
        .arg("1")
        .arg("-r")
        .arg("1000")
        .arg("-o")
        .arg(&rescaled_bam)
        .status()
        .expect("failed to execute tasmanian-rescale-quality");
    assert!(rescale_status.success(), "tasmanian-rescale-quality command failed");
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
        .spawn()
        .expect("failed to spawn rescale child");

    let rescale_pipe_status = rescale_child.wait().expect("failed waiting for rescale child");
    let mismatch_pipe_status = mismatch_child.wait().expect("failed waiting for mismatch child");
    let diagnostics_pipe_status = diagnostics_child.wait().expect("failed waiting for diagnostics child");

    assert!(diagnostics_pipe_status.success(), "diagnostics in pipe failed");
    assert!(mismatch_pipe_status.success(), "mismatch in pipe failed");
    assert!(rescale_pipe_status.success(), "rescale in pipe failed");
    assert!(piped_rescaled_bam.exists(), "missing piped rescaled BAM output");
}
