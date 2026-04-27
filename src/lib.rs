//! Tasmanian Mismatch library crate.
//!
//! This crate provides reusable building blocks for mismatch analysis on BAM
//! files, including reference loading, BED masking/filtering, methylation-aware
//! base normalization, and paired-read overlap handling.

/// BED parsing and region masking/filtering utilities.
pub mod bed;
/// FASTA and BAM input helpers.
pub mod io;
/// Insert-mode-specific processing and CLI support.
pub mod insert_mode;
/// Methylation-aware base adjustment logic.
pub mod methylation;
/// Core record-processing and mismatch-counting routines.
pub mod processing;
/// Shared types used throughout the library API.
pub mod types;
/// General helper functions for sequence and output handling.
pub mod utils;

// Re-export commonly used items
pub use bed::{
    filter_bed_for_region, mask_reference_with_bed, parse_bed_file, position_overlaps_intervals,
    BedInterval, BedRegions,
};
pub use io::{
    compute_read_len_max_from_sample_bam, load_reference_genome,
    load_rescaling_matrix, print_read_pair_inconsistency_table,
    write_inconsistencies_tsv, write_mismatch_discounts_tsv, write_potential_variants_tsv,
};
pub use methylation::adjust_methylation_base;
pub use processing::{
    compare_and_count, create_mismatch_key, get_overlap_region, process_overlap_region,
    process_paired_reads_with_overlap, process_record, process_single_record,
    rescale_phred_scores, merge_reads_into_insert_position_mode,
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
