use crate::bed::position_overlaps_intervals;
use crate::methylation::adjust_methylation_base;
use crate::types::{
    GenomicMismatchKey, GenomicMismatchValue, InconsistencyKey, MismatchKey, ReadInfo,
    ReferenceGenome,
};
use crate::utils::{base_to_char, calculate_end_pos, complement, correct_read_len_with_mode};
use rust_htslib::bam::Record;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy)]
pub struct ProcessingConfig {
    pub softclip_threshold: f64,
    pub min_base_quality: u8,
    pub is_methylation: bool,
    pub cpg_only: bool,
    pub mode_len: usize,
    pub min_map_quality: u8,
    pub required_flags: u16,
    pub filter_flags: u16,
    pub excl_flags: u16,
}

pub struct ProcessingContext<'a> {
    pub reference: &'a ReferenceGenome,
    pub tid_to_name: &'a HashMap<i32, String>,
    pub bed_intervals: &'a [crate::bed::BedInterval],
}

/// Create a MismatchKey with strand and methylation adjustments
pub fn create_mismatch_key(
    read_base: char,
    ref_base: char,
    r_pos: usize,
    seq_len: usize,
    is_reverse: bool,
    read_num: u8,
    is_methylation: bool,
    cpg_only: bool,
    ref_seq: &[u8],
    genome_pos: usize,
    mode_len: usize,
) -> MismatchKey {
    // Apply reverse complement if read is mapped to reverse strand
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

    // Apply methylation-aware base conversion if enabled
    let meth_and_strand_adjusted_read_base = adjust_methylation_base(
        strand_adjusted_read_base,
        strand_adjusted_ref_base,
        read_num,
        is_methylation,
        cpg_only,
        ref_seq,
        genome_pos,
    );
    let adjusted_r_pos = if is_reverse {
        seq_len - r_pos - 1
    } else {
        r_pos
    };
    let mode_adjusted_r_pos = correct_read_len_with_mode(adjusted_r_pos, seq_len, mode_len);

    MismatchKey {
        mismatch_type: format!(
            "{}>{}",
            strand_adjusted_ref_base, meth_and_strand_adjusted_read_base
        ),
        read_position: mode_adjusted_r_pos,
        read_num,
    }
}

/// Compare read base to reference and update mismatch counts
pub fn compare_and_count(
    seq: &rust_htslib::bam::record::Seq,
    qual: &[u8],
    ref_seq: &[u8],
    r_pos: usize,
    genome_pos: usize,
    is_reverse: bool,
    read_num: u8,
    min_base_quality: u8,
    is_methylation: bool,
    cpg_only: bool,
    local_counts: &mut HashMap<MismatchKey, usize>,
    genomic_counts: Option<&mut HashMap<GenomicMismatchKey, GenomicMismatchValue>>,
    chromosome: &str,
    mode_len: usize,
) {
    let seq_len = seq.len();
    let ref_len = ref_seq.len();

    if r_pos >= seq_len {
        eprintln!(
            "ERROR: Read position {} exceeds sequence length {} at chromosome {}:{}. This indicates a bug or corrupted data.",
            r_pos, seq_len, chromosome, genome_pos
        );
        return;
    }

    if genome_pos >= ref_len {
        eprintln!(
            "ERROR: Genomic position {} exceeds reference length {} on chromosome {}. This may indicate reference/BAM mismatch or corrupted data.",
            genome_pos, ref_len, chromosome
        );
        return;
    }

    // Check base quality score (skip if below threshold)
    if r_pos >= qual.len() || qual[r_pos] < min_base_quality {
        return;
    }

    let Some(read_base) = base_to_char(seq[r_pos]) else {
        return;
    };
    let Some(ref_base) = base_to_char(ref_seq[genome_pos]) else {
        return;
    };

    // Count both matches and mismatches
    let key = create_mismatch_key(
        read_base,
        ref_base,
        r_pos,
        seq_len,
        is_reverse,
        read_num,
        is_methylation,
        cpg_only,
        ref_seq,
        genome_pos,
        mode_len,
    );
    *local_counts.entry(key.clone()).or_insert(0) += 1;

    // Track genomic position for mismatches only (if genomic_counts is provided)
    if let Some(genomic_counts) = genomic_counts {
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
            let meth_adjusted_read_base = adjust_methylation_base(
                strand_adjusted_read_base,
                strand_adjusted_ref_base,
                read_num,
                is_methylation,
                cpg_only,
                ref_seq,
                genome_pos,
            );

            let genomic_key = GenomicMismatchKey {
                // e.g. key:value = {"chr1", "C>T", 123456}:  {{"C>T", 5, 1}, 23} vals=(_, read_position, counts)
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

/// Check if two reads overlap genomically
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
fn for_each_aligned_base_in_range<F>(
    record: &Record,
    range: Option<(usize, usize)>,
    mut on_base: F,
) where
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

#[inline]
fn should_skip_record(
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
    if record.is_unmapped() || record.mapq() <= min_map_quality {
        return true;
    }
    if skip_secondary_and_supplementary
        && (record.is_secondary() || record.is_supplementary())
    {
        return true;
    }

    false
}

/// Process overlap region between two reads
/// Detects both mismatches and inconsistencies between read pairs
/// If bed_intervals is provided, bases overlapping those intervals will be masked (skipped)
pub fn process_overlap_region(
    record1: &Record,
    record2: &Record,
    overlap_start: i64,
    overlap_end: i64,
    local_counts: &mut HashMap<MismatchKey, usize>,
    inconsistency_counts: &mut HashMap<InconsistencyKey, usize>,
    reference: &ReferenceGenome,
    tid_to_name: &HashMap<i32, String>,
    min_base_quality: u8,
    is_methylation: bool,
    cpg_only: bool,
    mode_len: usize,
    min_map_quality: u8,
    required_flags: u16,
    filter_flags: u16,
    excl_flags: u16,
    bed_intervals: &[crate::bed::BedInterval],
) {
    // Build genome_pos -> (read_pos, base, qual) mapping for both reads
    let mut read1_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();
    let mut read2_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();

    for (record, map) in [(record1, &mut read1_map), (record2, &mut read2_map)] {
        if should_skip_record(
            record,
            min_map_quality,
            required_flags,
            filter_flags,
            excl_flags,
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
            if r1_qual >= &min_base_quality && r2_qual >= &min_base_quality {
                if r1_base != r2_base {
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
    }

    // Process mismatch counts for both reads
    for record in [record1, record2] {
        if should_skip_record(
            record,
            min_map_quality,
            required_flags,
            filter_flags,
            excl_flags,
            false,
        ) {
            continue;
        }

        let tid = record.tid();
        let Some(chr_name) = tid_to_name.get(&tid) else {
            continue;
        };
        let Some(ref_seq) = reference.get(chr_name.as_str()) else {
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

        for_each_aligned_base_in_range(
            record,
            Some((overlap_start as usize, overlap_end as usize)),
            |r_pos, genome_pos| {
                // Skip if position overlaps BED region (mask mode)
                if is_bed_masked(bed_intervals, genome_pos) {
                    return;
                }

                compare_and_count(
                    &seq,
                    &qual,
                    ref_seq,
                    r_pos,
                    genome_pos,
                    record.is_reverse(),
                    read_num,
                    min_base_quality,
                    is_methylation,
                    cpg_only,
                    local_counts,
                    None,
                    chr_name,
                    mode_len,
                );
            },
        );
    }
}

/// Process a single BAM record and update mismatch counts
/// If bed_intervals is provided, bases overlapping those intervals will be masked (skipped)
pub fn process_record(
    record: &Record,
    local_counts: &mut HashMap<MismatchKey, usize>,
    mut genomic_counts: Option<&mut HashMap<GenomicMismatchKey, GenomicMismatchValue>>,
    reference: &ReferenceGenome,
    tid_to_name: &HashMap<i32, String>,
    softclip_threshold: f64,
    min_base_quality: u8,
    is_methylation: bool,
    cpg_only: bool,
    mode_len: usize,
    min_map_quality: u8,
    mut genomic_depth: Option<&mut HashMap<i64, usize>>,
    required_flags: u16,
    filter_flags: u16,
    excl_flags: u16,
    bed_intervals: &[crate::bed::BedInterval],
) {
    if should_skip_record(
        record,
        min_map_quality,
        required_flags,
        filter_flags,
        excl_flags,
        true,
    ) {
        return;
    }

    // Get chromosome name (target_id - defined from bam header) and other read info
    let tid = record.tid();
    if tid < 0 {
        return;
    } // unmapped
    let Some(chr_name) = tid_to_name.get(&tid) else {
        return;
    };
    let Some(ref_seq) = reference.get(chr_name.as_str()) else {
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
                    if is_bed_masked(bed_intervals, genome_pos) {
                        continue;
                    }

                    compare_and_count(
                        &seq,
                        &qual,
                        ref_seq,
                        r_pos,
                        genome_pos,
                        record.is_reverse(),
                        read_num,
                        min_base_quality,
                        is_methylation,
                        cpg_only,
                        local_counts,
                        genomic_counts.as_deref_mut(),
                        chr_name,
                        mode_len,
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
                // Soft-clipped bases - only count mismatches if at least 66% of bases match reference
                let mut total_bases = 0;
                let mut matching_bases = 0;
                let mut temp_mismatch_keys = Vec::<MismatchKey>::new();
                let seq_len = seq.len();

                for i in 0..*len {
                    let r_pos = read_pos + i as usize;
                    let genome_pos = if read_pos == 0 {
                        // beginning of read
                        if ref_pos >= (*len - i) as usize {
                            ref_pos - (*len - i) as usize
                        } else {
                            continue;
                        }
                    } else {
                        // end of read
                        ref_pos + i as usize
                    };

                    if r_pos >= seq_len || genome_pos >= ref_seq.len() {
                        continue;
                    }

                    // Skip if position overlaps BED region (mask mode)
                    if is_bed_masked(bed_intervals, genome_pos) {
                        continue;
                    }

                    // Skip bases with low quality
                    if r_pos >= qual.len() || qual[r_pos] < min_base_quality {
                        continue;
                    }

                    let read_base = match seq[r_pos] {
                        b'A' => 'A',
                        b'C' => 'C',
                        b'G' => 'G',
                        b'T' => 'T',
                        _ => continue,
                    };

                    let ref_base = match ref_seq[genome_pos] {
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
                        // Create MismatchKey using helper function
                        let key = create_mismatch_key(
                            read_base,
                            ref_base,
                            r_pos,
                            seq_len,
                            record.is_reverse(),
                            read_num,
                            is_methylation,
                            cpg_only,
                            ref_seq,
                            genome_pos,
                            mode_len,
                        );
                        temp_mismatch_keys.push(key);
                    }
                }

                // Only add mismatches to counts if at least 66% (default) match
                if total_bases > 0
                    && (matching_bases as f64 / total_bases as f64) >= softclip_threshold
                {
                    for key in temp_mismatch_keys {
                        *local_counts.entry(key).or_insert(0) += 1;
                    }
                }
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

/// Wrapper for processing a single record with standard parameters
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
        context.reference,
        context.tid_to_name,
        config.softclip_threshold,
        config.min_base_quality,
        config.is_methylation,
        config.cpg_only,
        config.mode_len,
        config.min_map_quality,
        Some(genomic_position_depth),
        config.required_flags,
        config.filter_flags,
        config.excl_flags,
        context.bed_intervals,
    );
}

/// Efficiently process a pair of overlapping reads in a single pass
/// Avoids redundant CIGAR walking for overlap positions
pub fn process_paired_reads_with_overlap(
    record1: &Record,
    record2: &Record,
    overlap_start: i64,
    overlap_end: i64,
    local_counts: &mut HashMap<MismatchKey, usize>,
    overlap_counts: &mut HashMap<MismatchKey, usize>,
    inconsistency_counts: &mut HashMap<InconsistencyKey, usize>,
    mut genomic_counts: Option<&mut HashMap<GenomicMismatchKey, GenomicMismatchValue>>,
    mut genomic_depth: Option<&mut HashMap<i64, usize>>,
    context: &ProcessingContext,
    config: ProcessingConfig,
) {
    // Build genome_pos -> (read_pos, base, qual) mapping for overlap detection
    let mut read1_overlap_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();
    let mut read2_overlap_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();

    // Process both reads
    for (record, is_first_read, overlap_map) in [
        (record1, true, &mut read1_overlap_map),
        (record2, false, &mut read2_overlap_map),
    ] {
        if should_skip_record(
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
                                    &seq,
                                    &qual,
                                    ref_seq,
                                    r_pos,
                                    genome_pos,
                                    record.is_reverse(),
                                    read_num,
                                    config.min_base_quality,
                                    config.is_methylation,
                                    config.cpg_only,
                                    overlap_counts,
                                    None, // Don't track genomic counts in overlap
                                    chr_name,
                                    config.mode_len,
                                );

                                // Update genomic depth for overlap positions
                                if let Some(depth) = genomic_depth.as_mut() {
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
                                &seq,
                                &qual,
                                ref_seq,
                                r_pos,
                                genome_pos,
                                record.is_reverse(),
                                read_num,
                                config.min_base_quality,
                                config.is_methylation,
                                config.cpg_only,
                                local_counts,
                                if is_first_read {
                                    genomic_counts.as_deref_mut()
                                } else {
                                    None
                                }, // Only track genomic counts once
                                chr_name,
                                config.mode_len,
                            );

                            // Update genomic depth
                            if let Some(depth) = genomic_depth.as_mut() {
                                depth.entry(genome_pos as i64).and_modify(|d| *d += 1).or_insert(1);
                            }
                        }
                    }
                    read_pos += *len as usize;
                    ref_pos += *len as usize;
                }
                SoftClip(len) => {
                    // Handle softclips (simplified - not processing to avoid complexity)
                    // In overlap regions, softclips are unlikely to be reliable anyway
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
            if r1_qual >= &config.min_base_quality && r2_qual >= &config.min_base_quality {
                if r1_base != r2_base {
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
    }
}

pub fn rescale_phred_scores(
    record: &mut Record,
    reference: &ReferenceGenome,
    tid_to_name: &HashMap<i32, String>,
    rescaling_matrix: &HashMap<(u8, u16, char, char), f32>,
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

    // Apply quality score modifications
    if !qual_modifications.is_empty() {
        let mut new_qual = qual.to_vec();
        for (pos, new_score) in qual_modifications {
            new_qual[pos] = new_score;
        }
        // Note: rust-htslib may not support modifying quality scores in-place
        // This would require reconstructing the record with modified quality scores
        // For now, we collect the modifications but cannot apply them without record.set()
    }
}
