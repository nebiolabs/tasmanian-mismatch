use std::collections::HashMap;

/// Key for tracking mismatches at specific read positions
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct MismatchKey {
    pub mismatch_type: String,  // e.g., "A>G", "C>T"
    pub read_position: usize,   // Position in the read (0-based)
    pub read_num: u8,           // 1 or 2 for paired-end reads
}

/// Key for tracking inconsistencies between read1 and read2 in overlap regions
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct InconsistencyKey {
    pub discordance_type: String,  // e.g., "R1:A_R2:G"
    pub read1_position: usize,
    pub read2_position: usize,
}

/// Store read information for overlap detection
#[derive(Debug, Clone)]
pub struct ReadInfo {
    pub tid: i32,
    pub pos: i64,
}

/// Store reference sequences in memory
pub type ReferenceGenome = HashMap<String, Vec<u8>>;

/// Key for tracking genomic mismatches
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct GenomicMismatchKey {
    pub mismatch_type: String,  // e.g., "A>G", "C>T"
    pub genomic_position: i64,  // (1-based)
    pub chromosome: String,
}