mod test_utils;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::index;
use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
use std::fs;
use std::path::Path;
use std::process::Command;
use test_utils::{log_command, log_line, repo_log_path, unique_temp_dir, write_two_chrom_bam};

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
        .arg("--opportunities-output")
        .arg(inconsistencies_tsv.with_extension("opportunities.tsv"))
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
        "r1_position\tr2_position\tdiscordance_type\tcount\n"
    );

    let discounts = fs::read_to_string(&discounts_tsv).expect("failed to read discount output");
    log_line(&log_path, &format!("discounts.tsv:\n{}", discounts));
    assert!(discounts.contains("mismatch_type\tread_num\tread_position\tdiscount_count"));
    assert!(
        // Same position tasmanian-mismatch reports for this mismatch (mismatch_integration.rs).
        discounts.lines().any(|line| line == "A>T\t1\t5\t1"),
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
        .arg("--opportunities-output")
        .arg(inconsistencies_tsv.with_extension("opportunities.tsv"))
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
        // Same position tasmanian-mismatch reports for this mismatch (mismatch_integration.rs).
        stdout.lines().any(|line| line == "A>T\t1\t5\t1"),
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
        .arg("--opportunities-output")
        .arg(inconsistencies_tsv.with_extension("opportunities.tsv"))
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
    assert!(inconsistencies.starts_with("r1_position\tr2_position\tdiscordance_type\tcount\n"));

    // Overlap is ref 4..8: read1 (forward) stored 4..7 pairs with read2 (reverse) stored 0..3.
    // On the 8-bp read axis, counted from each read's 5' end, that is R1 5..8 vs R2 8..5.
    // All Q40 and concordant, one comparison each.
    let opportunities = fs::read_to_string(inconsistencies_tsv.with_extension("opportunities.tsv"))
        .expect("failed to read opportunities output");
    assert_eq!(
        opportunities,
        "r1_position\tr2_position\tr1_base\tmin_quality_bin\tcompared\tdiscordant\n\
         5\t8\tA\t40\t1\t0\n\
         6\t7\tC\t40\t1\t0\n\
         7\t6\tG\t40\t1\t0\n\
         8\t5\tT\t40\t1\t0\n"
    );

    let discounts = fs::read_to_string(&discounts_tsv).expect("failed to read discounts output");
    assert!(discounts.contains("mismatch_type\tread_num\tread_position\tdiscount_count"));
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
        let opportunities_tsv = inconsistencies_tsv.with_extension("opportunities.tsv");
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
            "--opportunities-output",
            &opportunities_tsv.to_string_lossy(),
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

/// Deterministic pseudo-random reference of `len` bases.
fn synthetic_reference(len: usize) -> Vec<u8> {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn substitute(base: u8) -> u8 {
    match base {
        b'A' => b'G',
        b'C' => b'T',
        b'G' => b'A',
        _ => b'C',
    }
}

/// Overlapping 100-bp pairs tiled across a 6 kb contig, so every region boundary cuts
/// through pairs and sites. Shared "variant" sites are mutated in every read covering
/// them; each read also carries one private error, plus one R1/R2 disagreement per pair.
fn write_boundary_test_bam(path: &Path, reference: &[u8]) {
    const READ_LEN: usize = 100;
    const MATE_OFFSET: usize = 60;
    let variant_sites = [950usize, 1000, 1001, 1050, 2990, 3000, 4500];

    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", reference.len());
    header.push_record(&sq);
    let header_view = HeaderView::from_header(&header);

    let read_seq = |start: usize, private_error: usize| -> String {
        let mut seq = reference[start..start + READ_LEN].to_vec();
        for &site in &variant_sites {
            if (start..start + READ_LEN).contains(&site) {
                seq[site - start] = substitute(reference[site]);
            }
        }
        seq[private_error] = substitute(seq[private_error]);
        String::from_utf8(seq).unwrap()
    };

    let mut lines: Vec<(usize, String)> = Vec::new();
    let mut start = 0usize;
    let mut i = 0usize;
    while start + MATE_OFFSET + READ_LEN <= reference.len() {
        let mate = start + MATE_OFFSET;
        let template_len = MATE_OFFSET + READ_LEN;
        let qual = "I".repeat(READ_LEN);
        lines.push((
            start,
            format!(
                "p{i}\t99\tchr1\t{}\t60\t{READ_LEN}M\t=\t{}\t{template_len}\t{}\t{qual}\tMC:Z:{READ_LEN}M",
                start + 1,
                mate + 1,
                read_seq(start, (i * 7) % READ_LEN)
            ),
        ));
        lines.push((
            mate,
            format!(
                "p{i}\t147\tchr1\t{}\t60\t{READ_LEN}M\t=\t{}\t-{template_len}\t{}\t{qual}\tMC:Z:{READ_LEN}M",
                mate + 1,
                start + 1,
                read_seq(mate, (i * 13 + 3) % READ_LEN)
            ),
        ));
        start += 19;
        i += 1;
    }
    lines.sort_by_key(|(pos, _)| *pos);

    let mut writer =
        Writer::from_path(path, &header, Format::Bam).expect("failed to open BAM writer");
    for (_, line) in lines {
        let record = Record::from_sam(&header_view, line.as_bytes()).expect("bad SAM line");
        writer.write(&record).expect("failed to write record");
    }
}

fn sorted_lines(path: &Path) -> Vec<String> {
    let mut lines: Vec<String> = fs::read_to_string(path)
        .expect("failed to read output")
        .lines()
        .map(str::to_owned)
        .collect();
    lines.sort();
    lines
}

/// Regions are a unit of parallel work only: pairs and genomic sites that straddle region
/// boundaries must give the same tables whatever the region size.
#[test]
fn integration_diagnostics_outputs_do_not_depend_on_region_size() {
    let temp_dir = unique_temp_dir("diagnostics_region_size");
    let fixture_bam = temp_dir.join("input.bam");
    let reference_fa = temp_dir.join("reference.fa");

    let reference = synthetic_reference(6_000);
    fs::write(
        &reference_fa,
        format!(">chr1\n{}\n", String::from_utf8_lossy(&reference)),
    )
    .expect("failed to write reference");
    write_boundary_test_bam(&fixture_bam, &reference);
    index::build(&fixture_bam, None, index::Type::Bai, 1).expect("failed to build BAM index");

    let binary = env!("CARGO_BIN_EXE_tasmanian-diagnostics");
    let tables = ["variants", "inconsistencies", "opportunities", "discounts"];
    let run = |region_size: &str| -> Vec<Vec<String>> {
        let out = temp_dir.join(format!("r{region_size}"));
        fs::create_dir_all(&out).unwrap();
        let output = Command::new(binary)
            .args([
                "-t",
                "4",
                "-r",
                region_size,
                "-q",
                "0",
                "--min-map-quality",
                "0",
            ])
            .args(["--genomic-threshold", "2", "--genomic-depth-threshold", "1"])
            .arg("--variants-output")
            .arg(out.join("variants.tsv"))
            .arg("--inconsistencies-output")
            .arg(out.join("inconsistencies.tsv"))
            .arg("--opportunities-output")
            .arg(out.join("opportunities.tsv"))
            .arg("--discount-output")
            .arg(out.join("discounts.tsv"))
            .arg(&fixture_bam)
            .arg(&reference_fa)
            .output()
            .expect("failed to execute diagnostics binary");
        assert!(
            output.status.success(),
            "diagnostics -r {region_size} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        tables
            .iter()
            .map(|t| sorted_lines(&out.join(format!("{t}.tsv"))))
            .collect()
    };

    let whole = run("100000");
    // Shared sites pass the threshold; ones near 1000 and 3000 straddle the cuts below.
    let variant_positions: Vec<&str> = whole[0]
        .iter()
        .filter_map(|l| l.split('\t').nth(1))
        .collect();
    for site in ["950", "1000", "1001", "1050", "2990", "3000", "4500"] {
        assert!(variant_positions.contains(&site), "missing variant {site}");
    }
    assert!(whole[1].len() > 1, "expected R1/R2 inconsistencies");

    for region_size in ["1000", "1003", "250", "37"] {
        let split = run(region_size);
        for (t, (a, b)) in tables.iter().zip(whole.iter().zip(split.iter())) {
            assert_eq!(a, b, "{t} differs between -r 100000 and -r {region_size}");
        }
    }
}
