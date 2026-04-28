//! Tasmanian Mismatch library crate.
//!
//! This crate provides reusable building blocks for mismatch analysis on BAM
//! files, including reference loading, BED masking/filtering, methylation-aware
//! base normalization, and paired-read overlap handling.

/// BED parsing and region masking/filtering utilities.
pub mod bed;
/// FASTA and BAM input helpers.
pub mod io;
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
    apply_external_discounts, compute_read_len_max_from_sample_bam, frequencies_to_rescaling_matrix,
    load_discount_table, load_discount_table_from_reader, load_reference_genome,
    load_rescaling_matrix, load_rescaling_matrix_from_reader, normalize_mismatch_counts,
    print_read_pair_inconsistency_table, write_inconsistencies_tsv, write_mismatch_discounts_to_writer,
    write_mismatch_discounts_tsv,
    write_normalized_output, write_output, write_potential_variants_tsv,
    write_rescaling_matrix_output,
};
pub use processing::{
    base_position_for_mode,
    build_base_change, compare_record_to_reference, estimated_fragment_length, mc_mate_end,
    cigar_reference_span, qualifying_softclip_comparisons, softclip_identity,
    softclip_side_comparisons, build_tid_map_and_regions, configure_thread_pool,
    insert_mode_read_position, read_mode_read_position,
    overlap_interval, process_region, read_is_first_in_reference, record_read_num,
    should_skip_record, should_skip_whole_read_for_bed,
    compare_and_count, create_mismatch_key, get_overlap_region, process_overlap_region,
    process_paired_reads_with_overlap, process_record, process_single_record,
    rescale_phred_scores, merge_reads_into_insert_position_mode,
    ProcessingConfig, ProcessingContext,
};
pub use methylation::adjust_methylation_base;
pub use types::{
    Args, DiscountKey, GenomicMismatchKey, GenomicMismatchValue, InconsistencyKey, InsertKey,
    MismatchKey, OverlapMode, PositionMode, ReadInfo, ReferenceGenome, SoftclipComparison,
};
pub use utils::{
    base_to_char, calculate_end_pos, complement, correct_read_len_with_mode, parse_md_tag,
    print_main_output, print_position_table,
};
