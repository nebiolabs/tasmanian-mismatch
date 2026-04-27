//! Shared public types used across mismatch processing.

use clap::Parser;
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

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct InsertKey {
    pub base_change: String,
    pub read_num: u8,
    pub base_position: usize,
    /// 1 = first read in reference coordinates, 2 = second
    pub reference_order: u8,
}

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct DiscountKey {
    pub base_change: String,
    pub read_num: u8,
    pub base_position: usize,
}

#[derive(Debug, Clone)]
pub struct SoftclipComparison {
    pub read_pos: usize,
    pub ref_pos: usize,
    pub read_base: char,
    pub ref_base: char,
}

#[derive(Debug, Clone, clap::ValueEnum, PartialEq)]
pub enum OverlapMode {
    /// Split the overlap at its midpoint; read 1 keeps [ov_start, mid), read 2 keeps [mid, ov_end)
    Cut,
    /// Scale each read's positions to fill its half of [1, 2*max_read_len] via a cubic spline; no bases are dropped
    Stretch,
}

#[derive(Debug, Clone, Copy, clap::ValueEnum, PartialEq, Eq)]
pub enum PositionMode {
    Read,
    Insert,
}

#[derive(Parser, Debug)]
#[command(name = "tasmanian-mismatch")]
#[command(version, about = "Parallel indexed BAM mismatch caller by read and position")]
pub struct Args {
    /// Coordinate-sorted BAM file path (requires .bai)
    pub bam_path: String,

    /// Reference FASTA path
    pub reference_path: String,

    /// Number of threads to use (0 keeps rayon default)
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,

    /// Region size in bp for indexed chunking
    #[arg(short = 'r', long, default_value_t = 1_000_000)]
    pub region_size: usize,

    /// Minimum base quality
    #[arg(short = 'q', long, default_value_t = 20)]
    pub min_base_quality: u8,

    /// Minimum mapping quality
    #[arg(short = 'm', long, default_value_t = 10)]
    pub min_map_quality: u8,

    /// SAM flags that must be present
    #[arg(short = 'f', default_value_t = 0)]
    pub required_flags: u16,

    /// SAM flags that, if present, skip read
    #[arg(short = 'F', default_value_t = 0)]
    pub filter_flags: u16,

    /// SAM flags that, if all present, skip read
    #[arg(short = 'G', default_value_t = 0)]
    pub excl_flags: u16,

    /// How to handle read-pair overlapping positions: cut (split at midpoint) or stretch (cubic spline scale)
    #[arg(long, value_enum, default_value = "cut")]
    pub overlap_mode: OverlapMode,

    /// Write output to file instead of stdout
    #[arg(short = 'o', long)]
    pub output_file: Option<String>,

    /// Methylation mode for C>T and G>A changes: if true, treat them as C>C and G>G to avoid confounding by bisulfite conversion
    #[arg(long, default_value_t = false)]
    pub methylation_mode: bool,

    /// Restrict methylation collapsing to CpG context only.
    #[arg(long, default_value_t = false)]
    pub cpg_only: bool,

    /// minimum fragment length to be accepted for insert mode counting;
    #[arg(long, default_value_t = 0)]
    pub min_fragment_length: usize,

    /// maximum fragment length for insert mode counting
    #[arg(long, default_value_t = 1500)]
    pub max_fragment_length: usize,

    /// read position or fragment (insert) position mode
    #[arg(long, value_enum, default_value = "insert")]
    pub position_mode: PositionMode,

    /// Optional discount table emitted by tasmanian-diagnostics (variant_discounts.tsv).
    #[arg(long)]
    pub discount_table: Option<String>,

    /// BED file with regions to filter/mask (optional)
    #[arg(short = 'b', long)]
    pub bed_file: Option<String>,

    /// Filter mode: 'mask' (skip individual bases in BED regions) or 'filter' (skip whole reads overlapping BED regions)
    #[arg(long, default_value = "mask")]
    pub bed_filter_mode: String,
}
