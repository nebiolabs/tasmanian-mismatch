use rust_htslib::bam::{FetchDefinition, Format, Header, IndexedReader, Read, Reader, Record, Writer};
use std::collections::HashMap;
use std::sync::Arc;
use rustmanian_mismatch::*;
use clap::Parser;
use rayon::prelude::*;

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
    let args = Args::parse();
    let bam_path = &args.bam_file;
    let reference_path = &args.reference_fasta;
    let matrix_path = &args.matrix_file;
    let num_threads = args.threads;
    let region_size = args.region_size;
    let output_path = &args.output_file;

    // Set default thread count
    let actual_threads = if num_threads == 0 {
        num_cpus::get()
    } else {
        num_threads
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(actual_threads)
        .build_global()?;

    eprintln!("Rescaling quality scores in BAM file: {}", bam_path);
    eprintln!("Using reference: {}", reference_path);
    eprintln!("Using rescaling matrix: {}", matrix_path);
    eprintln!("Using {} threads", actual_threads);
    eprintln!("Region size: {} bp", region_size);

    // Load reference genome
    let reference = load_reference_genome(reference_path);
    eprintln!("Loaded {} contigs from reference", reference.len());

    // Load rescaling matrix
    let rescaling_matrix = load_rescaling_matrix(matrix_path)?;
    eprintln!("Loaded {} rescaling entries", rescaling_matrix.len());

    // Open BAM to read header and create regions
    let bam = Reader::from_path(bam_path)?;
    let header_view = bam.header();
    let out_header = Header::from_template(header_view);

    // Create thread-safe chromosome name mapping
    let tid_to_name: Arc<HashMap<i32, String>> = Arc::new(
        (0..header_view.target_count())
            .map(|i| {
                (
                    i as i32,
                    String::from_utf8_lossy(header_view.tid2name(i as u32)).to_string(),
                )
            })
            .collect(),
    );

    // Create regions to process in parallel
    let mut regions = Vec::new();
    for (&tid, chr_name) in tid_to_name.iter() {
        let chr_len = header_view.target_len(tid as u32).unwrap() as i64;

        let mut start = 0i64;
        while start < chr_len {
            let end = std::cmp::min(start + region_size as i64, chr_len);
            regions.push((tid, start, end, chr_name.clone()));
            start = end;
        }
    }

    // Keep write order deterministic and coordinate-like.
    regions.sort_by_key(|(tid, start, _, _)| (*tid, *start));

    eprintln!("Created {} regions to process", regions.len());
    drop(bam); // Close the initial BAM reader

    // Process regions in parallel
    let bam_path_arc = Arc::new(bam_path.to_string());
    let matrix_arc = Arc::new(rescaling_matrix);
    let tid_to_name_arc = Arc::clone(&tid_to_name);
    let reference_arc = Arc::new(reference);

    let processed_records: Vec<Vec<Record>> = regions
        .par_iter()
        .map(|(tid, start, end, chr_name)| {
            let mut region_records = Vec::new();

            match IndexedReader::from_path(bam_path_arc.as_str()) {
                Ok(mut bam) => {
                    if let Ok(_) = bam.fetch(FetchDefinition::Region(*tid, *start, *end)) {
                        for record_result in bam.records() {
                            if let Ok(mut record) = record_result {
                                // Avoid boundary duplicates from region-based fetch.
                                let rec_start = record.pos();
                                if rec_start < *start || rec_start >= *end {
                                    continue;
                                }

                                rescale_phred_scores(
                                    &mut record,
                                    &reference_arc,
                                    &tid_to_name_arc,
                                    &matrix_arc,
                                );
                                region_records.push(record);
                            }
                        }
                    }
                }
                Err(e) => {
                    eprintln!(
                        "Warning: Failed to fetch region {}:{}-{}: {}",
                        chr_name, start, end, e
                    );
                }
            }

            region_records
        })
        .collect();

    // Write all records in order
    eprintln!("Writing rescaled records to output...");
    let mut writer = if let Some(output_path) = output_path {
        Writer::from_path(output_path, &out_header, Format::Bam)?
    } else {
        Writer::from_stdout(&out_header, Format::Bam)?
    };

    let mut total_records = 0usize;
    for region_records in processed_records {
        for record in region_records {
            writer.write(&record)?;
            total_records += 1;
        }
    }

    // Include unmapped records to preserve whole-BAM output behavior.
    let mut full_reader = Reader::from_path(bam_path)?;
    for record_result in full_reader.records() {
        let mut record = record_result?;
        if record.tid() < 0 {
            rescale_phred_scores(
                &mut record,
                &reference_arc,
                &tid_to_name_arc,
                &matrix_arc,
            );
            writer.write(&record)?;
            total_records += 1;
        }
    }

    eprintln!("Finished rescaling {} quality scores.", total_records);
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
        let content = [
            "1\t10\tC\tT\t0.5",
            "2\t42\tG\tA\t1.25",
            "bad\tline",
            "",
        ]
        .join("\n");
        fs::write(&path, content).expect("failed to write temp matrix file");

        let matrix = load_rescaling_matrix(path.to_str().expect("invalid temp path"))
            .expect("failed to load matrix");

        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix.get(&(1, 10, 'C', 'T')).copied(), Some(0.5));
        assert_eq!(matrix.get(&(2, 42, 'G', 'A')).copied(), Some(1.25));

        let _ = fs::remove_file(path);
    }
}