// Tasmanian Mismatch - BAM file mismatch analysis tool
//
// This library provides functionality for analyzing mismatches between
// sequencing reads and a reference genome, with support for methylation-aware
// analysis and paired-end read overlap detection.

// Module declarations
pub mod bed;
pub mod io;
pub mod methylation;
pub mod processing;
pub mod types;
pub mod utils;

// Re-export commonly used items
pub use bed::{
    filter_bed_for_region, mask_reference_with_bed, parse_bed_file, position_overlaps_intervals,
    BedInterval, BedRegions,
};
pub use io::{compute_read_len_mode_from_sample_bam, load_reference_genome};
pub use methylation::adjust_methylation_base;
pub use processing::{
    compare_and_count, create_mismatch_key, get_overlap_region, process_overlap_region,
    process_paired_reads_with_overlap, process_record, process_single_record,
    ProcessingConfig, ProcessingContext,
};
pub use types::{
    GenomicMismatchKey, GenomicMismatchValue, InconsistencyKey, MismatchKey, ReadInfo,
    ReferenceGenome,
};
pub use utils::{
    base_to_char, calculate_end_pos, complement, correct_read_len_with_mode, parse_md_tag,
    print_main_output, print_position_table,
};
