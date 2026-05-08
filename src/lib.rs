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
    filter_bed_for_region, mask_reference_with_bed, maybe_parse_bed_file, parse_bed_file,
    position_overlaps_intervals, BedFilter, BedInterval, BedRegions,
};
pub use io::{
    apply_external_discounts, compute_read_len_max_from_sample_bam,
    frequencies_to_rescaling_matrix, load_discount_table, load_discount_table_from_reader,
    load_reference_genome, load_rescaling_matrix, load_rescaling_matrix_from_reader,
    normalize_mismatch_counts, print_read_pair_inconsistency_table, write_inconsistencies_tsv,
    write_mismatch_discounts_to_writer, write_mismatch_discounts_tsv, write_normalized_output,
    write_output, write_potential_variants_tsv, write_rescaling_matrix_output,
};
pub use methylation::adjust_methylation_base;
pub use processing::{
    base_position_for_mode, build_base_change, build_tid_map_and_regions, cigar_reference_span,
    compare_and_count, compare_record_to_reference, configure_thread_pool, create_mismatch_key,
    estimated_fragment_length, get_overlap_region, insert_mode_read_position, mc_mate_end,
    merge_reads_into_insert_position_mode, overlap_interval, process_overlap_region,
    process_paired_reads_with_overlap, process_record, process_region, process_single_record,
    qualifying_softclip_comparisons, read_is_first_in_reference, read_mode_read_position,
    record_read_num, rescale_phred_scores, should_skip_record, should_skip_whole_read_for_bed,
    softclip_identity, softclip_side_comparisons, ProcessingContext, ReadContext,
};
pub use types::{
    Args, DiscountKey, GenomicMismatchKey, GenomicMismatchValue, GenomicRegion, InconsistencyKey,
    InsertKey, MismatchKey, OverlapMode, PositionMode, ProcessingConfig, ReadInfo, ReferenceGenome,
    SoftclipComparison,
};
pub use utils::{
    base_to_char, calculate_end_pos, complement, correct_read_len_with_mode, parse_md_tag,
    print_main_output, print_position_table,
};
