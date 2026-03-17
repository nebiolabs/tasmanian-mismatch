//! Shared public types used across mismatch processing.

use std::collections::{HashMap, HashSet};

/// Key for counting mismatches at specific read positions.
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct MismatchKey {
    /// Mismatch class such as `A>G` or `C>T`.
    pub mismatch_type: String, // e.g., "A>G", "C>T"
    /// Position within the read in 0-based coordinates.
    pub read_position: usize,  // Position in the read (0-based)
    /// Read number within a pair.
    pub read_num: u8,          // 1 or 2 for paired-end reads
}

/// Key for counting read-pair inconsistencies inside overlap regions.
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct InconsistencyKey {
    /// Encoded discordant base combination, for example `R1:A_R2:G`.
    pub discordance_type: String, // e.g., "R1:A_R2:G"
    /// Read 1 position in 0-based coordinates.
    pub read1_position: usize,
    /// Read 2 position in 0-based coordinates.
    pub read2_position: usize,
}

/// Minimal read placement information for overlap detection.
#[derive(Debug, Clone)]
pub struct ReadInfo {
    /// Target identifier from the BAM header.
    pub tid: i32,
    /// Leftmost alignment position in 0-based genomic coordinates.
    pub pos: i64,
}

/// Reference genome sequences keyed by chromosome or contig name.
pub type ReferenceGenome = HashMap<String, Vec<u8>>;

/// Key for genomic mismatch aggregation.
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct GenomicMismatchKey {
    /// Chromosome or contig name.
    pub chromosome: String,    // Chromosome/contig name
    /// Mismatch class such as `A>G` or `C>T`.
    pub mismatch_type: String, // e.g., "A>G", "C>T"
    /// Genomic position associated with the event.
    pub genomic_position: i64, // (1-based)
}

/// Aggregated details for a genomic mismatch site.
#[derive(Debug, Clone)]
pub struct GenomicMismatchValue {
    /// Unique read-position mismatch keys observed at this genomic site.
    pub mismatch_keys: HashSet<MismatchKey>, // Set of all MismatchKey combinations that contributed
    /// Total number of unique mismatch keys contributing to the site.
    pub count: usize,                        // Total count (should match mismatch_keys.len())
}
