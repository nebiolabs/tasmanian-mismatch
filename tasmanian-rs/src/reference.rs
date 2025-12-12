//! FASTA reference genome loading

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

/// Reference genome sequences indexed by chromosome name
#[derive(Debug)]
pub struct Reference {
    sequences: HashMap<String, Vec<u8>>,
}

impl Reference {
    /// Load a FASTA reference file into memory
    ///
    /// The entire reference is loaded into memory for fast access.
    /// This matches the behavior of the Python implementation.
    pub fn load(path: &Path) -> Result<Self> {
        let start = Instant::now();

        let file = File::open(path)
            .with_context(|| format!("Failed to open reference file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
        let mut current_name: Option<String> = None;
        let mut current_seq: Vec<u8> = Vec::new();

        for line_result in reader.lines() {
            let line = line_result.context("Failed to read line from reference")?;

            if line.is_empty() {
                continue;
            }

            if line.starts_with('>') {
                // Save previous sequence if any
                if let Some(name) = current_name.take() {
                    sequences.insert(name, std::mem::take(&mut current_seq));
                }

                // Parse new sequence name (first word after '>')
                let name = line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                current_name = Some(name);
                current_seq = Vec::with_capacity(100_000_000); // Pre-allocate for typical chromosome
            } else if current_name.is_some() {
                // Append sequence data
                current_seq.extend(line.trim().as_bytes());
            }
        }

        // Save last sequence
        if let Some(name) = current_name {
            sequences.insert(name, current_seq);
        }

        let elapsed = start.elapsed();
        eprintln!(
            "Loaded reference in {:.2} seconds ({} sequences)",
            elapsed.as_secs_f64(),
            sequences.len()
        );

        Ok(Reference { sequences })
    }

    /// Get a slice of the reference sequence
    ///
    /// Returns None if chromosome doesn't exist or coordinates are out of bounds
    pub fn get(&self, chrom: &str, start: u64, end: u64) -> Option<&[u8]> {
        let seq = self.sequences.get(chrom)?;
        let start = start as usize;
        let end = end as usize;

        if start >= seq.len() || end > seq.len() {
            return None;
        }

        Some(&seq[start..end])
    }

    /// Get a single base from the reference
    pub fn get_base(&self, chrom: &str, pos: u64) -> Option<u8> {
        let seq = self.sequences.get(chrom)?;
        seq.get(pos as usize).copied()
    }

    /// Check if a chromosome exists in the reference
    pub fn has_chromosome(&self, chrom: &str) -> bool {
        self.sequences.contains_key(chrom)
    }

    /// Get the length of a chromosome
    pub fn chromosome_length(&self, chrom: &str) -> Option<usize> {
        self.sequences.get(chrom).map(|s| s.len())
    }

    /// Get the number of sequences in the reference
    pub fn num_sequences(&self) -> usize {
        self.sequences.len()
    }

    /// Get iterator over chromosome names
    pub fn chromosomes(&self) -> impl Iterator<Item = &String> {
        self.sequences.keys()
    }
}

/// Get the complement of a DNA base
#[inline]
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => base,
    }
}

/// Check if a base is valid (A, C, G, T in any case)
#[inline]
pub fn is_valid_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
    )
}

/// Convert a base to uppercase
#[inline]
pub fn to_uppercase(base: u8) -> u8 {
    match base {
        b'a' => b'A',
        b'c' => b'C',
        b'g' => b'G',
        b't' => b'T',
        _ => base,
    }
}

/// Convert a base to an index (0-3 for A, C, G, T)
#[inline]
pub fn base_to_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Convert an index (0-3) back to a base
#[inline]
pub fn index_to_base(index: usize) -> Option<u8> {
    match index {
        0 => Some(b'A'),
        1 => Some(b'C'),
        2 => Some(b'G'),
        3 => Some(b'T'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_temp_fasta(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_load_fasta() {
        let content = ">chr1\nACGT\nTGCA\n>chr2\nAAAA\n";
        let file = create_temp_fasta(content);
        let reference = Reference::load(file.path()).unwrap();

        assert_eq!(reference.num_sequences(), 2);
        assert!(reference.has_chromosome("chr1"));
        assert!(reference.has_chromosome("chr2"));
        assert!(!reference.has_chromosome("chr3"));
    }

    #[test]
    fn test_get_sequence() {
        let content = ">chr1\nACGTTGCA\n";
        let file = create_temp_fasta(content);
        let reference = Reference::load(file.path()).unwrap();

        assert_eq!(reference.get("chr1", 0, 4), Some(b"ACGT".as_slice()));
        assert_eq!(reference.get("chr1", 4, 8), Some(b"TGCA".as_slice()));
        assert_eq!(reference.get_base("chr1", 0), Some(b'A'));
        assert_eq!(reference.get_base("chr1", 1), Some(b'C'));
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'a'), b't');
    }

    #[test]
    fn test_base_to_index() {
        assert_eq!(base_to_index(b'A'), Some(0));
        assert_eq!(base_to_index(b'C'), Some(1));
        assert_eq!(base_to_index(b'G'), Some(2));
        assert_eq!(base_to_index(b'T'), Some(3));
        assert_eq!(base_to_index(b'a'), Some(0));
        assert_eq!(base_to_index(b'N'), None);
    }
}
