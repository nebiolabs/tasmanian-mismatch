use clap::Parser;
use rayon::prelude::*;
use std::path::Path;
use std::sync::Arc;
use tasmanian_mismatch::bam::{BamReader, BamWriter, IndexedBam, Record, RecordExt};
use tasmanian_mismatch::*;

#[derive(Parser, Debug)]
#[command(name = "tasmanian-rescale-quality")]
#[command(version, about, long_about = None)]
struct Args {
    /// Input BAM file (must be indexed)
    bam_file: String,

    /// Reference FASTA file
    reference_fasta: String,

    /// Rescaling matrix file (tab-separated: read_num, position, ref_base, read_base, scaling_factor).
    /// Use '-' to read matrix rows from stdin.
    matrix_file: String,

    #[arg(short = 'r', long, default_value_t = 10_000_000)]
    /// Region size for parallel processing (bp)
    region_size: u64,

    #[arg(short = 't', long, default_value_t = 0)]
    /// Number of threads (0 = auto-detect)
    threads: usize,

    #[arg(short = 'o', long)]
    /// Output BAM file (default: stdout)
    output_file: Option<String>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();
    let bam_path = &args.bam_file;
    let reference_path = &args.reference_fasta;
    let matrix_path = &args.matrix_file;
    let num_threads = args.threads;
    let region_size = args.region_size;
    let output_path = &args.output_file;

    // Set default thread count
    let actual_threads = if num_threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        num_threads
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(actual_threads)
        .build_global()?;

    log::info!("Rescaling quality scores in BAM file: {}", bam_path);
    log::info!("Using reference: {}", reference_path);
    log::info!("Using rescaling matrix: {}", matrix_path);
    log::info!("Using {} threads", actual_threads);
    log::info!("Region size: {} bp", region_size);

    // Load reference genome
    let reference = load_reference_genome(reference_path);
    log::info!("Loaded {} contigs from reference", reference.len());

    // Load rescaling matrix
    let rescaling_matrix = load_rescaling_matrix(matrix_path)?;
    log::info!("Loaded {} rescaling entries", rescaling_matrix.len());

    // Open BAM to read header and create regions
    let bam = IndexedBam::open(bam_path)?;
    let (tid_to_name, regions) = build_tid_map_and_regions(bam.header(), region_size as usize);
    log::info!("Created {} regions to process", regions.len());

    // Process regions in parallel
    let matrix_arc = Arc::new(rescaling_matrix);
    let reference_arc = Arc::new(reference);

    let processed_records: Vec<Vec<Record>> = regions
        .par_iter()
        .with_max_len(1)
        .map(|region| {
            let mut region_records = Vec::new();
            let result = bam.for_each_in_region(region, |mut record| {
                // Avoid boundary duplicates from region-based fetch.
                let rec_start = record.pos();
                if rec_start < region.start || rec_start >= region.end {
                    return;
                }

                rescale_phred_scores(&mut record, &reference_arc, &tid_to_name, &matrix_arc);
                region_records.push(record);
            });
            if let Err(e) = result {
                log::warn!(
                    "Failed to fetch region {}:{}-{}: {}",
                    tid_to_name[&region.tid],
                    region.start,
                    region.end,
                    e
                );
            }

            region_records
        })
        .collect();

    // Write all records in order
    log::info!("Writing rescaled records to output...");
    let mut writer =
        BamWriter::create(output_path.as_deref().map(Path::new), bam.header().clone())?;

    let mut total_records = 0usize;
    for region_records in processed_records {
        for record in region_records {
            writer.write(&record)?;
            total_records += 1;
        }
    }

    // Include unmapped records to preserve whole-BAM output behavior.
    let mut full_reader = BamReader::open(bam_path)?;
    for record_result in full_reader.records() {
        let mut record = record_result?;
        if record.tid() < 0 {
            rescale_phred_scores(&mut record, &reference_arc, &tid_to_name, &matrix_arc);
            writer.write(&record)?;
            total_records += 1;
        }
    }
    writer.finish()?;

    log::info!("Finished rescaling {} quality scores.", total_records);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(prefix: &str, ext: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!("{}_{}.{}", prefix, nanos, ext))
    }

    #[test]
    fn parse_args_defaults_stdout_and_threads() {
        let args = Args::parse_from([
            "tasmanian-rescale-quality",
            "input.bam",
            "ref.fa",
            "matrix.tsv",
        ]);

        assert_eq!(args.bam_file, "input.bam");
        assert_eq!(args.reference_fasta, "ref.fa");
        assert_eq!(args.matrix_file, "matrix.tsv");
        assert_eq!(args.threads, 0);
        assert_eq!(args.region_size, 10_000_000);
        assert!(args.output_file.is_none());
    }

    #[test]
    fn load_rescaling_matrix_reads_valid_rows_and_skips_short_lines() {
        let path = temp_path("rescaling_matrix", "tsv");
        let content = ["1\t10\tC\tT\t0.5", "2\t42\tG\tA\t1.25", "bad\tline", ""].join("\n");
        fs::write(&path, content).expect("failed to write temp matrix file");

        let matrix = load_rescaling_matrix(path.to_str().expect("invalid temp path"))
            .expect("failed to load matrix");

        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix.get(&(1, 10, 'C', 'T')).copied(), Some(0.5));
        assert_eq!(matrix.get(&(2, 42, 'G', 'A')).copied(), Some(1.25));

        let _ = fs::remove_file(path);
    }
}
