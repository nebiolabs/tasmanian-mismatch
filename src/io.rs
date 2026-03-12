use crate::types::ReferenceGenome;
use bio::io::fasta;
use rust_htslib::bam::{Read, Reader};
use std::collections::HashMap;

/// Load reference genome from FASTA file into memory
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

/// Compute mode (most common) read length from sample of BAM file
pub fn compute_read_len_mode_from_sample_bam(bam_path: &str, sample_size: usize) -> usize {
    let mut bam =
        Reader::from_path(bam_path).expect("Failed to open BAM file for read length sampling");

    let mut length_counts: HashMap<usize, usize> = HashMap::new();

    for (_i, result) in bam.records().enumerate().take(sample_size) {
        if let Ok(record) = result {
            let read_len = record.seq().len();
            *length_counts.entry(read_len).or_insert(0) += 1;
        }
    }

    // Find mode of read lengths
    let mode_length = length_counts
        .iter()
        .max_by_key(|&(_, count)| count)
        .map(|(&key, _)| key)
        .unwrap_or(0);

    eprintln!("Computed mode read length from sample: {}", mode_length);
    mode_length
}
