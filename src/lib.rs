// Tasmanian Mismatch - BAM file mismatch analysis tool
//
// This library provides functionality for analyzing mismatches between
// sequencing reads and a reference genome, with support for methylation-aware
// analysis and paired-end read overlap detection.

// Module declarations
pub mod types;
pub mod utils;
pub mod io;
pub mod methylation;
pub mod processing;

// Re-export commonly used items
pub use types::{MismatchKey, InconsistencyKey, ReadInfo, ReferenceGenome, GenomicMismatchKey, GenomicMismatchValue, GenomicPositionCounts};
pub use utils::{complement, base_to_char, correct_read_len_with_mode, calculate_end_pos, parse_md_tag};
pub use io::{load_reference_genome, compute_read_len_mode_from_sample_bam};
pub use methylation::adjust_methylation_base;
pub use processing::{create_mismatch_key, compare_and_count, get_overlap_region, process_overlap_region, process_record};
