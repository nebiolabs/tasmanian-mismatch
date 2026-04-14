//! Input helpers for reference FASTA and BAM-derived metadata.

use crate::types::GenomicMismatchKey;
use crate::types::ReferenceGenome;
use bio::io::fasta;
use rust_htslib::bam::{Read, Reader};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process;

/// Load a reference genome from a FASTA file.
///
/// # Arguments
/// * `fasta_path` - Path to the reference FASTA file.
///
/// # Returns
/// * A [`ReferenceGenome`] mapping chromosome names to base sequences.
///
/// # Panics
/// * If the FASTA file cannot be opened.
/// * If any FASTA record cannot be read.
pub fn load_reference_genome(fasta_path: &str) -> ReferenceGenome {
    eprintln!("Loading reference genome from: {}", fasta_path);
    let reader = fasta::Reader::from_file(fasta_path).expect("Failed to open reference FASTA file");

    let mut genome: ReferenceGenome = HashMap::new();
    for result in reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let chr_name = record.id().to_string();
        let sequence = record.seq().to_vec(); // Vec<u8> = byte, not UTF-8 char (overhead)
        genome.insert(chr_name, sequence);
    }

    eprintln!("Loaded {} chromosome(s)", genome.len());
    genome
}

/// Compute the mode (most common) read length from a sample of a BAM file.
///
/// # Arguments
/// * `bam_path` - Path to the input BAM file.
/// * `sample_size` - Number of records to sample from the start of the BAM stream.
///
/// # Returns
/// * The maximum read length (ideally the most common) in the sampled records.
/// * Returns a `default_len` if no valid records are observed in the sample.
///
/// # Panics
/// * If the BAM file cannot be opened.
pub fn compute_read_len_max_from_sample_bam(bam_path: &str, sample_size: usize) -> usize {
    let mut bam =
        Reader::from_path(bam_path).expect("Failed to open BAM file for read length sampling");

    let mut max_length: usize = 0;

    for result in bam.records().take(sample_size) {
        match result {
            Err(_) => {
                eprintln!(
                    "ERROR: Failed to sample BAM file for read length estimation. Exiting."
                );
                process::exit(1);
            }
            Ok(record) => {
                let read_len = record.seq().len();
                if read_len > max_length {
                    max_length = read_len;
                }
            }
        }
    }
    
    eprintln!("Computed max read length from sample: {}", max_length);
    max_length
}

/// Write genomic mismatch counts to a TSV file for potential variant review.
///
/// # Arguments
/// * `genomic_counts` - Genomic mismatch counts keyed by chromosome/position/mismatch.
/// * `output_path` - Output TSV path.
///
/// # Returns
/// * `Ok(())` on success.
/// * Any I/O error encountered while creating or writing the file.
pub fn write_potential_variants_tsv(
    genomic_counts: &HashMap<GenomicMismatchKey, (usize, usize)>,
    output_path: &str,
) -> std::io::Result<()> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth"
    )?;

    // Keep unsorted iteration for speed on large datasets.
    for (key, count) in genomic_counts.iter() {
        let parts: Vec<&str> = key.mismatch_type.split('>').collect();
        if parts.len() == 2 {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                key.chromosome, key.genomic_position, parts[0], parts[1], count.0, count.1
            )?;
        }
    }

    eprintln!(
        "Wrote {} genomic positions to {}",
        genomic_counts.len(),
        output_path
    );

    Ok(())
}

/// Print read pair inconsistency counts in CSV table format.
///
/// # Arguments
/// * `incons_types` - Sorted discordance type columns.
/// * `incons_positions` - Sorted `(read1_pos, read2_pos)` rows.
/// * `incons_position_map` - Mapping from position pairs to discordance counts.
pub fn print_read_pair_inconsistency_table(
    incons_types: &[String],
    incons_positions: &[(usize, usize)],
    incons_position_map: &HashMap<(usize, usize), HashMap<String, usize>>,
) {
    println!("\n# Read Pair Inconsistencies");
    print!("Read1_Pos,Read2_Pos");
    for incons_type in incons_types {
        print!(",{}", incons_type);
    }
    println!();

    for &(r1_pos, r2_pos) in incons_positions {
        print!("{},{}", r1_pos, r2_pos);
        if let Some(incons_map) = incons_position_map.get(&(r1_pos, r2_pos)) {
            for incons_type in incons_types {
                let count = incons_map.get(incons_type).unwrap_or(&0);
                print!(",{}", count);
            }
        }
        println!();
    }
}
