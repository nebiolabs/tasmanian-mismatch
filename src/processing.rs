//! Core mismatch-processing routines.
//!
//! This module contains the main logic for turning aligned BAM records into
//! mismatch, overlap, inconsistency, and genomic summary counts.

use crate::bed::position_overlaps_intervals;
use crate::methylation::adjust_methylation_base;
use crate::types::{
    GenomicMismatchKey, GenomicMismatchValue, GenomicRegion, InconsistencyKey, InsertKey,
    MismatchKey, OverlapMode, PositionMode, ProcessingConfig, ReadInfo, ReferenceGenome,
    ReferenceOrder, RescalingMatrix, SoftclipComparison,
};
use crate::utils::{base_to_char, calculate_end_pos, complement, correct_read_len_with_mode};
use rayon;
use rust_htslib::bam::{record::Aux, FetchDefinition, IndexedReader, Read, Reader, Record};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// Maps BAM target ID (tid) to chromosome name.
type TidNameMap = Arc<HashMap<i32, String>>;
/// List of contiguous genomic regions to process.
type RegionList = Vec<GenomicRegion>;

/// Shared references needed while processing records within a region.
pub struct ProcessingContext<'a> {
    /// Reference genome sequences keyed by contig name.
    pub reference: &'a ReferenceGenome,
    /// Mapping from BAM target ID to chromosome or contig name.
    pub tid_to_name: &'a HashMap<i32, String>,
    /// BED intervals relevant to the current processing region.
    pub bed_intervals: &'a [crate::bed::BedInterval],
}

/// Per-read state needed when comparing bases against the reference.
pub struct ReadContext<'a> {
    /// Encoded read sequence.
    pub seq: &'a rust_htslib::bam::record::Seq<'a>,
    /// Base quality scores aligned to `seq`.
    pub qual: &'a [u8],
    /// Reference sequence for the current contig.
    pub ref_seq: &'a [u8],
    /// Whether the read aligns to the reverse strand.
    pub is_reverse: bool,
    /// Read number within the pair (1 or 2).
    pub read_num: u8,
}

/// Mutable count accumulators used when processing overlapping read pairs.
pub struct OverlapCounts<'a> {
    /// Non-overlap mismatch counts.
    pub local_counts: &'a mut HashMap<MismatchKey, usize>,
    /// Overlap-only mismatch counts.
    pub overlap_counts: &'a mut HashMap<MismatchKey, usize>,
    /// Read-pair inconsistency counts within the overlap.
    pub inconsistency_counts: &'a mut HashMap<InconsistencyKey, usize>,
    /// Optional genomic mismatch summary.
    pub genomic_counts: Option<&'a mut HashMap<GenomicMismatchKey, GenomicMismatchValue>>,
    /// Optional per-position depth counter.
    pub genomic_depth: Option<&'a mut HashMap<i64, usize>>,
}

/// Canonicalize a mismatch by ensuring the reference base is earlier in the A<C<G<T order.
///
/// This ensures complement pairs like C>T and G>A both report as C>T.
/// Returns the canonical mismatch in ref>alt form.
fn canonicalize_mismatch(ref_base: char, read_base: char) -> String {
    let complement_ref = complement(ref_base);
    let complement_read = complement(read_base);

    let base_order = |b: char| match b {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => 4,
    };

    if base_order(ref_base) <= base_order(complement_ref) {
        format!("{}>{}", ref_base, read_base)
    } else {
        format!("{}>{}", complement_ref, complement_read)
    }
}

/// Build a [`MismatchKey`] after strand and methylation normalization.
///
/// # Arguments
/// * `read_base` - Observed base in the read.
/// * `ref_base` - Reference base at the aligned position.
/// * `r_pos` - Position within the read.
/// * `genome_pos` - Genomic position within `ref_seq`.
/// * `read_ctx` - Per-read state (sequence, quality, reference, strand, read number).
/// * `config` - Shared processing configuration.
///
/// # Returns
/// * A normalized [`MismatchKey`] suitable for counting.
pub fn create_mismatch_key(
    read_base: char,
    ref_base: char,
    r_pos: usize,
    read_ctx: &ReadContext,
    config: &ProcessingConfig,
) -> MismatchKey {
    let seq_len = read_ctx.seq.len();
    // Apply reverse complement if read is mapped to reverse strand
    let strand_adjusted_read_base = if read_ctx.is_reverse {
        complement(read_base)
    } else {
        read_base
    };
    let strand_adjusted_ref_base = if read_ctx.is_reverse {
        complement(ref_base)
    } else {
        ref_base
    };

    // Apply methylation-aware base conversion if enabled
    let meth_and_strand_adjusted_read_base = adjust_methylation_base(
        strand_adjusted_read_base,
        strand_adjusted_ref_base,
        read_ctx.read_num,
        config.is_methylation,
    );
    let adjusted_r_pos = if read_ctx.is_reverse {
        seq_len - r_pos - 1
    } else {
        r_pos
    };
    let mode_adjusted_r_pos = correct_read_len_with_mode(
        adjusted_r_pos,
        seq_len,
        config.mode_len,
        config.use_insert_mode,
        read_ctx.read_num,
    );

    MismatchKey {
        mismatch_type: format!(
            "{}>{}",
            strand_adjusted_ref_base, meth_and_strand_adjusted_read_base
        ),
        read_position: mode_adjusted_r_pos,
        read_num: read_ctx.read_num,
    }
}

/// Compare a read base to the reference and update mismatch aggregates.
///
/// # Arguments
/// * `read_ctx` - Per-read state (sequence, quality, reference, strand, read number).
/// * `r_pos` - Position within the read.
/// * `genome_pos` - Position within the reference sequence.
/// * `config` - Shared processing configuration.
/// * `local_counts` - Per-region mismatch counts to update.
/// * `genomic_counts` - Optional genomic mismatch summary to update.
/// * `chromosome` - Chromosome or contig name, used for genomic summaries.
///
/// # Panics
/// * This function does not panic intentionally, but it emits errors to stderr
///   and returns early when coordinates fall outside the supplied sequence data.
pub fn compare_and_count(
    read_ctx: &ReadContext,
    r_pos: usize,
    genome_pos: usize,
    config: &ProcessingConfig,
    local_counts: &mut HashMap<MismatchKey, usize>,
    genomic_counts: Option<&mut HashMap<GenomicMismatchKey, GenomicMismatchValue>>,
    chromosome: &str,
) {
    let seq_len = read_ctx.seq.len();
    let ref_len = read_ctx.ref_seq.len();

    if r_pos >= seq_len {
        log::error!(
            "Read position {} exceeds sequence length {} at chromosome {}:{}. This indicates a bug or corrupted data.",
            r_pos, seq_len, chromosome, genome_pos
        );
        return;
    }

    if genome_pos >= ref_len {
        log::error!(
            "Genomic position {} exceeds reference length {} on chromosome {}. This may indicate reference/BAM mismatch or corrupted data.",
            genome_pos, ref_len, chromosome
        );
        return;
    }

    // Check base quality score (skip if below threshold)
    if r_pos >= read_ctx.qual.len() || read_ctx.qual[r_pos] < config.min_base_quality {
        return;
    }

    let Some(read_base) = base_to_char(read_ctx.seq[r_pos]) else {
        return;
    };
    let Some(ref_base) = base_to_char(read_ctx.ref_seq[genome_pos]) else {
        return;
    };

    // Count both matches and mismatches
    let key = create_mismatch_key(read_base, ref_base, r_pos, read_ctx, config);
    *local_counts.entry(key.clone()).or_insert(0) += 1;

    // Track genomic position for mismatches only (if genomic_counts is provided)
    if let Some(genomic_counts) = genomic_counts {
        if read_base != ref_base {
            let strand_adjusted_read_base = if read_ctx.is_reverse {
                complement(read_base)
            } else {
                read_base
            };
            let strand_adjusted_ref_base = if read_ctx.is_reverse {
                complement(ref_base)
            } else {
                ref_base
            };
            let meth_adjusted_read_base = adjust_methylation_base(
                strand_adjusted_read_base,
                strand_adjusted_ref_base,
                read_ctx.read_num,
                config.is_methylation,
            );

            let genomic_key = GenomicMismatchKey {
                // e.g. key:value = {"chr1", "C>T", 123456}:  {{"C>T", 5, 1}, 23} vals=(_, read_position, readnum), counts
                chromosome: chromosome.to_string(),
                mismatch_type: format!("{}>{}", strand_adjusted_ref_base, meth_adjusted_read_base),
                genomic_position: genome_pos as i64,
            };

            let genomic_values =
                genomic_counts
                    .entry(genomic_key)
                    .or_insert_with(|| GenomicMismatchValue {
                        mismatch_keys: HashSet::new(),
                        count: 0,
                    });
            genomic_values.mismatch_keys.insert(key);
            genomic_values.count = genomic_values.mismatch_keys.len();
        }
    }
}

/// Determine the genomic overlap between two reads.
///
/// # Arguments
/// * `read1` - Placement summary for the first read.
/// * `read2` - Placement summary for the second read.
/// * `record1` - Full BAM record for the first read.
/// * `record2` - Full BAM record for the second read.
///
/// # Returns
/// * `Some((start, end))` when the reads overlap on the same contig.
/// * `None` when they are on different contigs or do not overlap.
pub fn get_overlap_region(
    read1: &ReadInfo,
    read2: &ReadInfo,
    record1: &Record,
    record2: &Record,
) -> Option<(i64, i64)> {
    if read1.tid != read2.tid {
        return None;
    }

    let end1 = calculate_end_pos(read1.pos, &record1.cigar());
    let end2 = calculate_end_pos(read2.pos, &record2.cigar());

    let overlap_start = read1.pos.max(read2.pos);
    let overlap_end = end1.min(end2);

    if overlap_start < overlap_end {
        Some((overlap_start, overlap_end))
    } else {
        None
    }
}

/// Iterate over read/reference-aligned bases for Match/Equal/Diff CIGAR operations.
/// If `range` is provided, only emit positions in `[start, end)` genomic coordinates.
fn for_each_aligned_base_in_range<F>(record: &Record, range: Option<(usize, usize)>, mut on_base: F)
where
    F: FnMut(usize, usize),
{
    let mut read_pos = 0usize;
    let mut ref_pos = record.pos() as usize;

    use rust_htslib::bam::record::Cigar::*;
    for cigar_op in record.cigar().iter() {
        match cigar_op {
            Match(len) | Equal(len) | Diff(len) => {
                for i in 0..*len {
                    let r_pos = read_pos + i as usize;
                    let genome_pos = ref_pos + i as usize;

                    if let Some((start, end)) = range {
                        if genome_pos < start || genome_pos >= end {
                            continue;
                        }
                    }

                    on_base(r_pos, genome_pos);
                }
                read_pos += *len as usize;
                ref_pos += *len as usize;
            }
            Ins(len) => {
                read_pos += *len as usize;
            }
            Del(len) | RefSkip(len) => {
                ref_pos += *len as usize;
            }
            _ => {}
        }
    }
}

#[inline]
fn is_bed_masked(bed_intervals: &[crate::bed::BedInterval], genome_pos: usize) -> bool {
    !bed_intervals.is_empty() && position_overlaps_intervals(bed_intervals, genome_pos as i64)
}

/// Count eligible soft-clip mismatches for a single soft-clip CIGAR operation.
///
/// Mismatches are only added when the soft-clipped segment meets `softclip_threshold`
/// match rate against the reference.
fn count_softclip_mismatches(
    read_ctx: &ReadContext,
    read_pos: usize,
    ref_pos: usize,
    softclip_len: u32,
    config: &ProcessingConfig,
    bed_intervals: &[crate::bed::BedInterval],
    local_counts: &mut HashMap<MismatchKey, usize>,
) {
    let mut total_bases = 0usize;
    let mut matching_bases = 0usize;
    let mut temp_mismatch_keys = Vec::<MismatchKey>::new();
    let seq_len = read_ctx.seq.len();

    for i in 0..softclip_len {
        let r_pos = read_pos + i as usize;
        let genome_pos = if read_pos == 0 {
            // Soft-clip at beginning of read.
            if ref_pos >= (softclip_len - i) as usize {
                ref_pos - (softclip_len - i) as usize
            } else {
                continue;
            }
        } else {
            // Soft-clip at end of read.
            ref_pos + i as usize
        };

        if r_pos >= seq_len || genome_pos >= read_ctx.ref_seq.len() {
            continue;
        }

        if is_bed_masked(bed_intervals, genome_pos) {
            continue;
        }

        if r_pos >= read_ctx.qual.len() || read_ctx.qual[r_pos] < config.min_base_quality {
            continue;
        }

        let read_base = match read_ctx.seq[r_pos] {
            b'A' => 'A',
            b'C' => 'C',
            b'G' => 'G',
            b'T' => 'T',
            _ => continue,
        };

        let ref_base = match read_ctx.ref_seq[genome_pos] {
            b'A' | b'a' => 'A',
            b'C' | b'c' => 'C',
            b'G' | b'g' => 'G',
            b'T' | b't' => 'T',
            _ => continue,
        };

        total_bases += 1;
        if read_base == ref_base {
            matching_bases += 1;
        } else {
            let key = create_mismatch_key(read_base, ref_base, r_pos, read_ctx, config);
            temp_mismatch_keys.push(key);
        }
    }

    if total_bases > 0 && (matching_bases as f64 / total_bases as f64) >= config.softclip_threshold
    {
        for key in temp_mismatch_keys {
            *local_counts.entry(key).or_insert(0) += 1;
        }
    }
}

#[inline]
fn should_skip_record_core(
    record: &Record,
    min_map_quality: u8,
    required_flags: u16,
    filter_flags: u16,
    excl_flags: u16,
    skip_secondary_and_supplementary: bool,
) -> bool {
    let flags = record.flags();

    // samtools-like flag filtering
    if required_flags != 0 && (flags & required_flags) != required_flags {
        return true;
    }
    if filter_flags != 0 && (flags & filter_flags) != 0 {
        return true;
    }
    if excl_flags != 0 && (flags & excl_flags) == excl_flags {
        return true;
    }

    // Standard read quality/status filtering
    if record.is_unmapped() || record.mapq() < min_map_quality {
        return true;
    }
    if skip_secondary_and_supplementary && (record.is_secondary() || record.is_supplementary()) {
        return true;
    }

    false
}

/// Process the overlapping segment of a read pair.
///
/// This counts overlap mismatches and read-pair inconsistencies while honoring
/// mapping-quality, flag, methylation, and BED-masking settings.
///
/// # Arguments
/// * `record1` - First BAM record in the pair.
/// * `record2` - Second BAM record in the pair.
/// * `overlap_start` - Inclusive overlap start position.
/// * `overlap_end` - Exclusive overlap end position.
/// * `local_counts` - Mismatch counts to update.
/// * `inconsistency_counts` - Overlap inconsistency counts to update.
/// * `context` - Shared reference, target-name, and BED context.
/// * `config` - Shared processing configuration.
pub fn process_overlap_region(
    record1: &Record,
    record2: &Record,
    overlap: (i64, i64),
    local_counts: &mut HashMap<MismatchKey, usize>,
    inconsistency_counts: &mut HashMap<InconsistencyKey, usize>,
    context: &ProcessingContext,
    config: &ProcessingConfig,
) {
    let (overlap_start, overlap_end) = overlap;
    // Build genome_pos -> (read_pos, base, qual) mapping for both reads
    let mut read1_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();
    let mut read2_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();

    for (record, map) in [(record1, &mut read1_map), (record2, &mut read2_map)] {
        if should_skip_record_core(
            record,
            config.min_map_quality,
            config.required_flags,
            config.filter_flags,
            config.excl_flags,
            false,
        ) {
            continue;
        }

        let seq = record.seq();
        let qual = record.qual();

        for_each_aligned_base_in_range(
            record,
            Some((overlap_start as usize, overlap_end as usize)),
            |r_pos, genome_pos| {
                if r_pos < seq.len() && r_pos < qual.len() {
                    map.insert(genome_pos, (r_pos, seq[r_pos], qual[r_pos]));
                }
            },
        );
    }

    // Compare bases at overlapping genomic positions
    for (genome_pos, (r1_pos, r1_base, r1_qual)) in read1_map.iter() {
        if let Some((r2_pos, r2_base, r2_qual)) = read2_map.get(genome_pos) {
            // Both reads cover this position - check for inconsistency
            if r1_qual >= &config.min_base_quality
                && r2_qual >= &config.min_base_quality
                && r1_base != r2_base
            {
                // Bases disagree - record inconsistency
                let r1_char = base_to_char(*r1_base).unwrap_or('N');
                let r2_char = base_to_char(*r2_base).unwrap_or('N');

                let key = InconsistencyKey {
                    discordance_type: format!("R1:{}_R2:{}", r1_char, r2_char),
                    read1_position: *r1_pos,
                    read2_position: *r2_pos,
                };
                *inconsistency_counts.entry(key).or_insert(0) += 1;
            }
        }
    }

    // Process mismatch counts for both reads
    for record in [record1, record2] {
        if should_skip_record_core(
            record,
            config.min_map_quality,
            config.required_flags,
            config.filter_flags,
            config.excl_flags,
            false,
        ) {
            continue;
        }

        let tid = record.tid();
        let Some(chr_name) = context.tid_to_name.get(&tid) else {
            continue;
        };
        let Some(ref_seq) = context.reference.get(chr_name.as_str()) else {
            continue;
        };
        let seq = record.seq();
        let qual = record.qual();
        let read_num = if record.is_first_in_template() {
            1
        } else if record.is_last_in_template() {
            2
        } else {
            1
        };

        let read_ctx = ReadContext {
            seq: &seq,
            qual,
            ref_seq,
            is_reverse: record.is_reverse(),
            read_num,
        };

        for_each_aligned_base_in_range(
            record,
            Some((overlap_start as usize, overlap_end as usize)),
            |r_pos, genome_pos| {
                // Skip if position overlaps BED region (mask mode)
                if is_bed_masked(context.bed_intervals, genome_pos) {
                    return;
                }

                compare_and_count(
                    &read_ctx,
                    r_pos,
                    genome_pos,
                    config,
                    local_counts,
                    None,
                    chr_name,
                );
            },
        );
    }
}

/// Process a single BAM record and update mismatch-related counts.
///
/// This handles aligned positions and eligible soft clips, optionally updating
/// both per-position mismatch counts and per-genomic-site summaries.
///
/// # Arguments
/// * `record` - BAM record to process.
/// * `local_counts` - Mismatch counts to update.
/// * `genomic_counts` - Optional genomic mismatch summary to update.
/// * `context` - Shared reference, target-name, and BED context.
/// * `config` - Shared processing configuration.
/// * `genomic_depth` - Optional per-position depth counter.
pub fn process_record(
    record: &Record,
    local_counts: &mut HashMap<MismatchKey, usize>,
    mut genomic_counts: Option<&mut HashMap<GenomicMismatchKey, GenomicMismatchValue>>,
    context: &ProcessingContext,
    config: &ProcessingConfig,
    mut genomic_depth: Option<&mut HashMap<i64, usize>>,
) {
    if should_skip_record_core(
        record,
        config.min_map_quality,
        config.required_flags,
        config.filter_flags,
        config.excl_flags,
        true,
    ) {
        return;
    }

    // Get chromosome name (target_id - defined from bam header) and other read info
    let tid = record.tid();
    if tid < 0 {
        return;
    } // unmapped
    let Some(chr_name) = context.tid_to_name.get(&tid) else {
        return;
    };
    let Some(ref_seq) = context.reference.get(chr_name.as_str()) else {
        return;
    };
    let ref_start = record.pos() as usize; // (0-based)
    let seq = record.seq();
    let qual = record.qual();
    let cigar = record.cigar();
    let read_num = if record.is_first_in_template() {
        1
    } else if record.is_last_in_template() {
        2
    } else {
        1
    };

    let read_ctx = ReadContext {
        seq: &seq,
        qual,
        ref_seq,
        is_reverse: record.is_reverse(),
        read_num,
    };

    // Single CIGAR pass: handle aligned bases and soft clips together.
    let mut read_pos = 0usize;
    let mut ref_pos = ref_start;
    for cigar_op in cigar.iter() {
        use rust_htslib::bam::record::Cigar::*;
        match cigar_op {
            Match(len) | Equal(len) | Diff(len) => {
                for i in 0..*len {
                    let r_pos = read_pos + i as usize;
                    let genome_pos = ref_pos + i as usize;

                    // Skip if position overlaps BED region (mask mode)
                    if is_bed_masked(context.bed_intervals, genome_pos) {
                        continue;
                    }

                    compare_and_count(
                        &read_ctx,
                        r_pos,
                        genome_pos,
                        config,
                        local_counts,
                        genomic_counts.as_deref_mut(),
                        chr_name,
                    );

                    // Update genomic depth counts.
                    if let Some(genomic_depth) = genomic_depth.as_mut() {
                        *genomic_depth.entry(genome_pos as i64).or_insert(0) += 1;
                    }
                }

                read_pos += *len as usize;
                ref_pos += *len as usize;
            }
            SoftClip(len) => {
                count_softclip_mismatches(
                    &read_ctx,
                    read_pos,
                    ref_pos,
                    *len,
                    config,
                    context.bed_intervals,
                    local_counts,
                );
                read_pos += *len as usize;
            }
            Ins(len) => {
                read_pos += *len as usize;
            }
            Del(len) | RefSkip(len) => {
                ref_pos += *len as usize;
            }
            HardClip(_) | Pad(_) => {}
        }
    }
}

/// Process a single record using a shared [`ProcessingContext`] and [`ProcessingConfig`].
#[inline]
pub fn process_single_record(
    record: &Record,
    region_counts: &mut HashMap<MismatchKey, usize>,
    genomic_region_counts: &mut HashMap<GenomicMismatchKey, GenomicMismatchValue>,
    genomic_position_depth: &mut HashMap<i64, usize>,
    context: &ProcessingContext,
    config: ProcessingConfig,
) {
    process_record(
        record,
        region_counts,
        Some(genomic_region_counts),
        context,
        &config,
        Some(genomic_position_depth),
    );
}

/// Process an overlapping read pair in a single pass.
///
/// This avoids redundant CIGAR traversal in the shared region while updating
/// standard mismatch counts, overlap-specific counts, inconsistencies, genomic
/// summaries, and optional depth information.
///
/// # Arguments
/// * `record1` - First BAM record in the pair.
/// * `record2` - Second BAM record in the pair.
/// * `overlap` - Inclusive start and exclusive end of the overlap region.
/// * `counts` - Mutable count accumulators for mismatches, overlaps, and inconsistencies.
/// * `context` - Shared reference, target-name, and BED context.
/// * `config` - Shared processing configuration.
pub fn process_paired_reads_with_overlap(
    record1: &Record,
    record2: &Record,
    overlap: (i64, i64),
    counts: &mut OverlapCounts,
    context: &ProcessingContext,
    config: ProcessingConfig,
) {
    let (overlap_start, overlap_end) = overlap;
    // Build genome_pos -> (read_pos, base, qual) mapping for overlap detection
    let mut read1_overlap_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();
    let mut read2_overlap_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();

    // Process both reads
    for (record, is_first_read, overlap_map) in [
        (record1, true, &mut read1_overlap_map),
        (record2, false, &mut read2_overlap_map),
    ] {
        if should_skip_record_core(
            record,
            config.min_map_quality,
            config.required_flags,
            config.filter_flags,
            config.excl_flags,
            true,
        ) {
            continue;
        }

        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let Some(chr_name) = context.tid_to_name.get(&tid) else {
            continue;
        };
        let Some(ref_seq) = context.reference.get(chr_name.as_str()) else {
            continue;
        };

        let ref_start = record.pos() as usize;
        let seq = record.seq();
        let qual = record.qual();
        let cigar = record.cigar();
        let read_num = if record.is_first_in_template() {
            1
        } else if record.is_last_in_template() {
            2
        } else {
            1
        };

        let read_ctx = ReadContext {
            seq: &seq,
            qual,
            ref_seq,
            is_reverse: record.is_reverse(),
            read_num,
        };

        let mut read_pos = 0;
        let mut ref_pos = ref_start;

        // Walk through CIGAR operations
        for cigar_op in cigar.iter() {
            use rust_htslib::bam::record::Cigar::*;
            match cigar_op {
                Match(len) | Equal(len) | Diff(len) => {
                    for i in 0..*len {
                        let r_pos = read_pos + i as usize;
                        let genome_pos = ref_pos + i as usize;

                        // Skip if position overlaps BED region (mask mode)
                        if is_bed_masked(context.bed_intervals, genome_pos) {
                            continue;
                        }

                        let is_in_overlap = genome_pos >= overlap_start as usize
                            && genome_pos < overlap_end as usize;

                        if is_in_overlap {
                            // Store for inconsistency detection
                            if r_pos < seq.len() && r_pos < qual.len() {
                                overlap_map.insert(genome_pos, (r_pos, seq[r_pos], qual[r_pos]));
                            }

                            // Count mismatches in overlap (only process once - use first read)
                            if is_first_read {
                                compare_and_count(
                                    &read_ctx,
                                    r_pos,
                                    genome_pos,
                                    &config,
                                    counts.overlap_counts,
                                    None, // Don't track genomic counts in overlap
                                    chr_name,
                                );

                                // Update genomic depth for overlap positions
                                if let Some(depth) = counts.genomic_depth.as_mut() {
                                    // Count both reads covering this position
                                    depth
                                        .entry(genome_pos as i64)
                                        .and_modify(|d| *d += 2)
                                        .or_insert(2);
                                }
                            }
                        } else {
                            // Not in overlap - count normally
                            compare_and_count(
                                &read_ctx,
                                r_pos,
                                genome_pos,
                                &config,
                                counts.local_counts,
                                if is_first_read {
                                    counts.genomic_counts.as_deref_mut()
                                } else {
                                    None
                                }, // Only track genomic counts once
                                chr_name,
                            );

                            // Update genomic depth
                            if let Some(depth) = counts.genomic_depth.as_mut() {
                                depth
                                    .entry(genome_pos as i64)
                                    .and_modify(|d| *d += 1)
                                    .or_insert(1);
                            }
                        }
                    }
                    read_pos += *len as usize;
                    ref_pos += *len as usize;
                }
                SoftClip(len) => {
                    count_softclip_mismatches(
                        &read_ctx,
                        read_pos,
                        ref_pos,
                        *len,
                        &config,
                        context.bed_intervals,
                        counts.local_counts,
                    );
                    read_pos += *len as usize;
                }
                Ins(len) => {
                    read_pos += *len as usize;
                }
                Del(len) | RefSkip(len) => {
                    ref_pos += *len as usize;
                }
                HardClip(_) | Pad(_) => {}
            }
        }
    }

    // Detect inconsistencies in the overlap region
    for (genome_pos, (r1_pos, r1_base, r1_qual)) in read1_overlap_map.iter() {
        if let Some((r2_pos, r2_base, r2_qual)) = read2_overlap_map.get(genome_pos) {
            if r1_qual >= &config.min_base_quality
                && r2_qual >= &config.min_base_quality
                && r1_base != r2_base
            {
                let r1_char = base_to_char(*r1_base).unwrap_or('N');
                let r2_char = base_to_char(*r2_base).unwrap_or('N');

                let key = InconsistencyKey {
                    discordance_type: format!("R1:{}_R2:{}", r1_char, r2_char),
                    read1_position: *r1_pos,
                    read2_position: *r2_pos,
                };
                *counts.inconsistency_counts.entry(key).or_insert(0) += 1;
            }
        }
    }
}

/// Collect quality-score rescaling adjustments for mismatch positions.
///
/// # Arguments
/// * `record` - BAM record whose quality scores should be examined.
/// * `reference` - Reference genome sequences.
/// * `tid_to_name` - BAM target ID to contig-name mapping.
/// * `rescaling_matrix` - Per-read-position mismatch scaling factors.
///
/// # Notes
/// The current implementation computes quality-score modifications but does not
/// write them back into the record.
pub fn rescale_phred_scores(
    record: &mut Record,
    reference: &ReferenceGenome,
    tid_to_name: &HashMap<i32, String>,
    rescaling_matrix: &RescalingMatrix,
) {
    let tid = record.tid();
    let Some(chr_name) = tid_to_name.get(&tid) else {
        return;
    };
    let Some(ref_seq) = reference.get(chr_name.as_str()) else {
        return;
    };
    let ref_start = record.pos() as usize;
    let seq = record.seq();
    let read_num = if record.is_first_in_template() {
        1
    } else if record.is_last_in_template() {
        2
    } else {
        1
    };
    let cigar = record.cigar(); // There could be INDELS...
    let qual = record.qual();
    let is_reverse = record.is_reverse();
    let mut read_pos = 0usize;
    let mut ref_pos = ref_start;

    // Collect quality score modifications
    let mut qual_modifications: Vec<(usize, u8)> = Vec::new();

    for cigar_op in cigar.iter() {
        use rust_htslib::bam::record::Cigar::*;
        match cigar_op {
            Match(len) | Equal(len) | Diff(len) => {
                for i in 0..*len {
                    let r_pos = read_pos + i as usize;
                    let genome_pos = ref_pos + i as usize;

                    let Some(read_base) = base_to_char(seq[r_pos]) else {
                        continue;
                    };
                    let Some(ref_base) = base_to_char(ref_seq[genome_pos]) else {
                        continue;
                    };

                    let phred_score = qual[r_pos];

                    if read_base != ref_base {
                        let strand_adjusted_read_base = if is_reverse {
                            complement(read_base)
                        } else {
                            read_base
                        };
                        let strand_adjusted_ref_base = if is_reverse {
                            complement(ref_base)
                        } else {
                            ref_base
                        };

                        let key = (
                            read_num,
                            r_pos as u16,
                            strand_adjusted_ref_base,
                            strand_adjusted_read_base,
                        );

                        if let Some(&scaling_factor) = rescaling_matrix.get(&key) {
                            let new_phred = (phred_score as f32 * scaling_factor).round() as u8;
                            qual_modifications.push((r_pos, new_phred.min(40)));
                            // Cap at max Phred score of 40
                        }
                    }
                }
                read_pos += *len as usize;
                ref_pos += *len as usize;
            }
            Ins(len) => {
                read_pos += *len as usize;
            }
            Del(len) | RefSkip(len) => {
                ref_pos += *len as usize;
            }
            _ => {}
        }
    }

    // Rebuild variable-length record data with the rescaled quality values.
    if !qual_modifications.is_empty() {
        let mut new_qual = qual.to_vec();
        for (pos, new_score) in qual_modifications {
            new_qual[pos] = new_score;
        }

        let qname = record.qname().to_vec();
        let cigar = record.cigar().take();
        let seq = record.seq().as_bytes();

        record.set(&qname, Some(&cigar), &seq, &new_qual);
    }
}

pub fn merge_reads_into_insert_position_mode(
    region_counts: &HashMap<MismatchKey, usize>,
    max_len: usize, // read_max_len
) -> HashMap<(String, usize), usize> {
    let mut insert_position_counts: HashMap<(String, usize), usize> = HashMap::new();

    // example key:val = MismatchKey { mismatch_type: "C>T", read_position: 12, read_num: 1 }: 21
    for (key, count) in region_counts.iter() {
        let insert_pos = if key.read_num == 1 {
            key.read_position
        } else {
            max_len * 2 + 10 - key.read_position // 10 is arbitrary separation between reads
        };

        log::debug!(
            "Processing key: {:?}, count: {}, insert_pos: {}",
            key,
            count,
            insert_pos
        );

        // Canonicalize the mismatch type (e.g., both "C>T" and "G>A" become "C>T")
        let canonical_mismatch =
            if key.mismatch_type.len() == 3 && key.mismatch_type.chars().nth(1) == Some('>') {
                let chars: Vec<char> = key.mismatch_type.chars().collect();
                canonicalize_mismatch(chars[0], chars[2])
            } else {
                key.mismatch_type.clone()
            };

        let insert_key = (canonical_mismatch, insert_pos);
        *insert_position_counts.entry(insert_key).or_insert(0) += count;
    }
    insert_position_counts
}

pub fn insert_mode_read_position(
    order: ReferenceOrder,
    read_pos: usize,
    seq_len: usize,
    max_read_len: usize,
    stretch: bool,
    fragment_len: Option<usize>,
) -> usize {
    let short_fragment = fragment_len.is_some_and(|fl| fl <= max_read_len - 10);
    let trailing_bases = seq_len - (read_pos + 1);

    if stretch && seq_len > 1 {
        if order == ReferenceOrder::First {
            // Cubic smooth-step: s = t^2(3 - 2t), maps [0,1] to [0,1].let t = read_pos as f64 / (seq_len - 1) as f64;
            let t = read_pos as f64 / (seq_len - 1) as f64;
            let s = t * t * (3.0 - 2.0 * t);
            1 + (s * (max_read_len - 1) as f64).round() as usize
        } else {
            let trailing = seq_len - 1 - read_pos;
            let t = trailing as f64 / (seq_len - 1) as f64;
            let s = t * t * (3.0 - 2.0 * t);
            2 * max_read_len - (s * (max_read_len - 1) as f64).round() as usize
        }
    } else if order == ReferenceOrder::First {
        if read_pos >= seq_len / 2 && short_fragment {
            2 * max_read_len - trailing_bases
        } else {
            read_pos + 1
        }
    } else if read_pos >= seq_len / 2 && short_fragment {
        read_pos + 1
    } else {
        2 * max_read_len - trailing_bases
    }
}

pub fn read_mode_read_position(
    pos: usize,
    seq_len: usize,
    is_reverse: bool,
    max_read_len: usize,
) -> usize {
    let half = seq_len.div_ceil(2);

    match (is_reverse, pos <= half) {
        (true, true) => max_read_len - pos,
        (true, false) => seq_len - pos,
        (false, false) => pos + (max_read_len - seq_len) + 1,
        (false, true) => pos + 1,
    }
}

pub fn base_position_for_mode(
    config: &ProcessingConfig,
    read_pos: usize,
    seq_len: usize,
    is_reverse: bool,
    reference_order: ReferenceOrder,
    stretch: bool,
    fragment_len: Option<usize>,
) -> usize {
    match config.position_mode {
        PositionMode::Read => {
            read_mode_read_position(read_pos, seq_len, is_reverse, config.mode_len)
        }
        PositionMode::Insert => insert_mode_read_position(
            reference_order,
            read_pos,
            seq_len,
            config.mode_len,
            stretch,
            fragment_len,
        ),
    }
}

pub fn configure_thread_pool(threads: usize) {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Failed to configure rayon thread pool");
    }
}

pub fn build_tid_map_and_regions(bam_path: &str, region_size: usize) -> (TidNameMap, RegionList) {
    let bam = Reader::from_path(bam_path).expect("Failed to open BAM file");

    let tid_to_name: TidNameMap = Arc::new(
        (0..bam.header().target_count())
            .map(|i| {
                (
                    i as i32,
                    String::from_utf8_lossy(bam.header().tid2name(i)).to_string(),
                )
            })
            .collect(),
    );

    let mut regions = Vec::new();
    for &tid in tid_to_name.keys() {
        let chr_len = bam.header().target_len(tid as u32).unwrap_or(0) as i64;
        let mut start = 0i64;
        while start < chr_len {
            let end = std::cmp::min(start + region_size as i64, chr_len);
            regions.push(GenomicRegion { tid, start, end });
            start = end;
        }
    }

    (tid_to_name, regions)
}

pub fn cigar_reference_span(cigar: &str) -> Option<usize> {
    let mut ref_bases = 0usize;
    let mut run_len = 0usize;

    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            run_len = run_len * 10 + (b - b'0') as usize;
            continue;
        }

        match b {
            b'M' | b'=' | b'X' | b'D' | b'N' => {
                ref_bases += run_len;
            }
            b'I' | b'S' | b'H' | b'P' => {}
            _ => return None,
        }
        run_len = 0;
    }

    if run_len != 0 {
        return None;
    }

    Some(ref_bases)
}

pub fn softclip_side_comparisons(
    record: &Record,
    ref_seq: &[u8],
    read_start: usize,
    clip_len: usize,
    ref_start: usize,
) -> Vec<SoftclipComparison> {
    if clip_len == 0
        || read_start + clip_len > record.seq().len()
        || ref_start + clip_len > ref_seq.len()
    {
        return Vec::new();
    }

    let read_seq = record.seq();
    let mut comparisons = Vec::with_capacity(clip_len);

    for offset in 0..clip_len {
        let read_index = read_start + offset;
        let ref_index = ref_start + offset;

        let Some(read_base) = base_to_char(read_seq[read_index]) else {
            continue;
        };
        let Some(ref_base) = base_to_char(ref_seq[ref_index].to_ascii_uppercase()) else {
            continue;
        };

        comparisons.push(SoftclipComparison {
            read_pos: read_index,
            ref_pos: ref_index,
            read_base,
            ref_base,
        });
    }

    comparisons
}

pub fn softclip_identity(comparisons: &[SoftclipComparison]) -> Option<f64> {
    if comparisons.is_empty() {
        None
    } else {
        let matches = comparisons
            .iter()
            .filter(|comparison| comparison.read_base == comparison.ref_base)
            .count();
        Some(matches as f64 / comparisons.len() as f64)
    }
}

pub fn qualifying_softclip_comparisons(
    record: &Record,
    ref_seq: &[u8],
    min_identity: f64,
) -> Vec<SoftclipComparison> {
    let cigar = record.cigar();
    let aligned_start = record.pos();
    let aligned_end = cigar.end_pos();
    let mut comparisons = Vec::new();

    if let Some(rust_htslib::bam::record::Cigar::SoftClip(len)) = cigar.iter().next() {
        if aligned_start >= *len as i64 {
            let side = softclip_side_comparisons(
                record,
                ref_seq,
                0,
                *len as usize,
                (aligned_start - *len as i64) as usize,
            );
            if softclip_identity(&side).is_some_and(|identity| identity >= min_identity) {
                comparisons.extend(side);
            }
        }
    }

    if let Some(rust_htslib::bam::record::Cigar::SoftClip(len)) = cigar.iter().last() {
        let side = softclip_side_comparisons(
            record,
            ref_seq,
            record.seq().len().saturating_sub(*len as usize),
            *len as usize,
            aligned_end as usize,
        );
        if softclip_identity(&side).is_some_and(|identity| identity >= min_identity) {
            comparisons.extend(side);
        }
    }

    comparisons
}

/// Parse the MC tag and return the mate's exclusive end position, or `None`
/// when the MC tag is absent, the record is unpaired / mate-unmapped, or the
/// reads are on different contigs.
pub fn mc_mate_end(record: &Record) -> Option<i64> {
    if !record.is_paired() || record.is_mate_unmapped() {
        return None;
    }
    if record.tid() < 0 || record.mtid() < 0 || record.tid() != record.mtid() {
        return None;
    }
    let mate_start = record.mpos();
    let mate_span = match record.aux(b"MC".as_ref()) {
        Ok(Aux::String(mc)) => i64::try_from(cigar_reference_span(mc)?).ok()?,
        _ => return None,
    };
    Some(mate_start + mate_span)
}

pub fn estimated_fragment_length(record: &Record, mate_end: Option<i64>) -> Option<usize> {
    let mate_end = mate_end?;

    let read_start = record.pos();
    let read_end = record.cigar().end_pos();
    let mate_start = record.mpos();

    let fragment_start = std::cmp::min(read_start, mate_start);
    let fragment_end = std::cmp::max(read_end, mate_end);
    if fragment_end <= fragment_start {
        return None;
    }

    usize::try_from(fragment_end - fragment_start).ok()
}

pub fn should_skip_record(record: &Record, config: ProcessingConfig) -> bool {
    if record.is_unmapped() || record.mapq() < config.min_map_quality {
        return true;
    }
    let flags = record.flags();
    let missing_required =
        config.required_flags != 0 && (flags & config.required_flags) != config.required_flags;
    let has_filtered = config.filter_flags != 0 && (flags & config.filter_flags) != 0;
    let has_excluded = config.excl_flags != 0 && (flags & config.excl_flags) == config.excl_flags;
    missing_required || has_filtered || has_excluded
}

pub fn should_skip_whole_read_for_bed(
    record: &Record,
    bed_filter_whole_reads: bool,
    chunk_bed_intervals: &[crate::bed::BedInterval],
    bed_cursor: &mut usize,
) -> bool {
    if !bed_filter_whole_reads || chunk_bed_intervals.is_empty() {
        return false;
    }

    let read_start = record.pos();
    let read_end = calculate_end_pos(read_start, &record.cigar());

    // Reads are fetched in coordinate order, so we can advance a cursor
    // and never revisit intervals that end before this read starts.
    while *bed_cursor < chunk_bed_intervals.len()
        && chunk_bed_intervals[*bed_cursor].end < read_start
    {
        *bed_cursor += 1;
    }

    let mut idx = *bed_cursor;
    while idx < chunk_bed_intervals.len() {
        let interval = &chunk_bed_intervals[idx];

        if interval.start > read_end {
            break;
        }

        if read_start <= interval.end && read_end >= interval.start {
            return true;
        }

        idx += 1;
    }

    false
}

pub fn record_read_num(record: &Record) -> u8 {
    if record.is_first_in_template() {
        1
    } else if record.is_last_in_template() {
        2
    } else {
        1
    }
}

pub fn read_is_first_in_reference(record: &Record) -> ReferenceOrder {
    if !record.is_paired() || record.is_mate_unmapped() {
        return ReferenceOrder::First;
    }

    let read_coord = (record.tid(), record.pos());
    let mate_coord = (record.mtid(), record.mpos());

    if read_coord < mate_coord {
        ReferenceOrder::First
    } else if read_coord > mate_coord {
        ReferenceOrder::Second
    } else if record_read_num(record) == 1 {
        ReferenceOrder::First
    } else {
        ReferenceOrder::Second
    }
}

pub fn overlap_interval(record: &Record, mate_end: Option<i64>) -> Option<(usize, usize)> {
    let mate_end = mate_end?;

    let read_start = record.pos();
    let read_end = record.cigar().end_pos();
    let mate_start = record.mpos();

    let ov_start = std::cmp::max(read_start, mate_start);
    let ov_end = std::cmp::min(read_end, mate_end);
    if ov_start < ov_end {
        Some((ov_start as usize, ov_end as usize))
    } else {
        None
    }
}

pub fn build_base_change(
    read_num: u8,
    ref_base: char,
    read_base: char,
    is_reverse: bool,
    methylation_mode: bool,
) -> String {
    let (strand_ref, strand_read) = if is_reverse {
        (complement(ref_base), complement(read_base))
    } else {
        (ref_base, read_base)
    };

    if !methylation_mode {
        return format!("{}>{}", strand_ref, strand_read);
    }

    match (read_num, ref_base, read_base, is_reverse) {
        (1, 'C', 'T', false) | (2, 'G', 'A', false) | (1, 'G', 'A', true) | (2, 'C', 'T', true) => {
            format!("{}>{}", strand_ref, strand_ref)
        }
        _ => format!("{}>{}", strand_ref, strand_read),
    }
}

pub fn compare_record_to_reference(
    record: &Record,
    context: &ProcessingContext,
    config: ProcessingConfig,
    mate_end: Option<i64>,
    local_counts: &mut HashMap<InsertKey, usize>,
) {
    let tid = record.tid();
    if tid < 0 {
        return;
    }

    let Some(chr_name) = context.tid_to_name.get(&tid) else {
        return;
    };
    let Some(ref_seq) = context.reference.get(chr_name.as_str()) else {
        return;
    };

    let seq = record.seq();
    let qual = record.qual();
    let seq_len = seq.len();
    let read_num = record_read_num(record);
    let reference_order = read_is_first_in_reference(record);
    let overlap = overlap_interval(record, mate_end);
    let softclip_comparisons = qualifying_softclip_comparisons(record, ref_seq, 0.66);
    let stretch = config.overlap_mode == OverlapMode::Stretch;

    let mut read_pos = 0usize;
    let mut ref_pos = record.pos() as usize;

    use rust_htslib::bam::record::Cigar::*;
    for op in record.cigar().iter() {
        match op {
            Match(len) | Equal(len) | Diff(len) => {
                for i in 0..*len as usize {
                    let rp = read_pos + i;
                    let gp = ref_pos + i;

                    if config.overlap_mode == OverlapMode::Cut {
                        if let Some((ov_start, ov_end)) = overlap {
                            if gp >= ov_start && gp < ov_end && read_num == 2 {
                                continue;
                            }
                        }
                    }

                    if rp >= seq_len || rp >= qual.len() || gp >= ref_seq.len() {
                        continue;
                    }
                    if qual[rp] < config.min_base_quality {
                        continue;
                    }

                    let Some(read_base) = base_to_char(seq[rp]) else {
                        continue;
                    };
                    let Some(ref_base) = base_to_char(ref_seq[gp].to_ascii_uppercase()) else {
                        continue;
                    };

                    let ref_base_change = build_base_change(
                        read_num,
                        ref_base,
                        read_base,
                        record.is_reverse(),
                        config.is_methylation,
                    );

                    let base_position = base_position_for_mode(
                        &config,
                        rp,
                        seq_len,
                        record.is_reverse(),
                        reference_order,
                        stretch,
                        estimated_fragment_length(record, mate_end),
                    );

                    *local_counts
                        .entry(InsertKey {
                            base_change: ref_base_change,
                            read_num,
                            base_position,
                            reference_order,
                        })
                        .or_insert(0) += 1;
                }

                read_pos += *len as usize;
                ref_pos += *len as usize;
            }
            Ins(len) | SoftClip(len) => {
                read_pos += *len as usize;
            }
            Del(len) | RefSkip(len) => {
                ref_pos += *len as usize;
            }
            HardClip(_) | Pad(_) => {}
        }
    }

    let frag_len = estimated_fragment_length(record, mate_end);
    for comparison in softclip_comparisons {
        if comparison.read_pos >= qual.len() || qual[comparison.read_pos] < config.min_base_quality
        {
            continue;
        }

        let ref_base_change = build_base_change(
            read_num,
            comparison.ref_base,
            comparison.read_base,
            record.is_reverse(),
            config.is_methylation,
        );
        let base_position = base_position_for_mode(
            &config,
            comparison.read_pos,
            seq_len,
            record.is_reverse(),
            reference_order,
            stretch,
            frag_len,
        );
        *local_counts
            .entry(InsertKey {
                base_change: ref_base_change,
                read_num,
                base_position,
                reference_order,
            })
            .or_insert(0) += 1;
    }
}

pub fn process_region(
    bam_path: &str,
    region: &GenomicRegion,
    context: &ProcessingContext,
    config: ProcessingConfig,
    bed_filter: &crate::bed::BedFilter<'_>,
) -> (HashMap<InsertKey, usize>, usize) {
    let chr_name = context
        .tid_to_name
        .get(&region.tid)
        .expect("tid not found in tid_to_name");
    let mut bam = IndexedReader::from_path(bam_path).expect("Failed to open indexed BAM");
    if let Err(error) = bam.fetch(FetchDefinition::Region(
        region.tid,
        region.start,
        region.end,
    )) {
        log::warn!(
            "Failed to fetch region {}:{}-{}: {}",
            chr_name,
            region.start,
            region.end,
            error
        );
        return (HashMap::new(), 0);
    }

    let mut local_counts: HashMap<InsertKey, usize> = HashMap::new();
    let mut local_read_count = 0usize;
    let chunk_bed_intervals = bed_filter
        .regions
        .map(|bed| crate::filter_bed_for_region(bed, chr_name, region.start, region.end))
        .unwrap_or_default();
    let mut bed_cursor = 0usize;

    for result in bam.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        // Ensure each read is counted once across chunks.
        let read_start = record.pos();
        if read_start < region.start || read_start >= region.end {
            continue;
        }

        if should_skip_whole_read_for_bed(
            &record,
            bed_filter.filter_whole_reads,
            &chunk_bed_intervals,
            &mut bed_cursor,
        ) {
            continue;
        }

        let mate_end = mc_mate_end(&record);

        if should_skip_record(&record, config) {
            continue;
        }

        compare_record_to_reference(&record, context, config, mate_end, &mut local_counts);
        local_read_count += 1;
    }

    (local_counts, local_read_count)
}
