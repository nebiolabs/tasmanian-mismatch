//! Core artifact analysis engine
//!
//! Processes SAM records and counts mismatches against the reference.

use crate::error_table::ErrorTables;
use crate::reference::{complement, to_uppercase, Reference};
use crate::sam::record::{parse_cigar, TasmanianRecord};
use crate::sam::flags::{ReadNumber, Strand};
use anyhow::Result;
use log::debug;
use noodles::sam::Header;
use std::collections::HashMap;
use std::time::Instant;

/// Configuration for artifact analysis
#[derive(Debug, Clone)]
pub struct AnalysisConfig {
    pub min_base_quality: u8,
    pub min_mapping_quality: u8,
    pub filter_indels: bool,
    pub min_read_length: u32,
    pub max_read_length: u32,
    pub min_fragment_length: u32,
    pub max_fragment_length: u32,
    pub softclip_bypass: u8,
    pub unmask_genome: bool,
    pub confidence_threshold: u32,
    pub ont_mode: bool,
    pub mask_methyl_c: bool,
    pub mask_methyl_cpg: bool,
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            min_base_quality: 20,
            min_mapping_quality: 20,
            filter_indels: false,
            min_read_length: 0,
            max_read_length: 350,
            min_fragment_length: 0,
            max_fragment_length: 10000,
            softclip_bypass: 0,
            unmask_genome: false,
            confidence_threshold: 20,
            ont_mode: false,
            mask_methyl_c: false,
            mask_methyl_cpg: false,
        }
    }
}

/// Statistics for skipped reads
#[derive(Debug, Default)]
pub struct SkipStats {
    pub indel: u64,
    pub map_quality: u64,
    pub frag_length: u64,
    pub read_length: u64,
    pub other: u64,
    pub softclips: u64,
}

/// Statistics for skipped bases
#[derive(Debug, Default)]
pub struct BaseSkipStats {
    pub quality: u64,
    pub softclips: u64,
}

/// Run artifact analysis on SAM records
pub fn analyze<R: std::io::BufRead>(
    reader: &mut noodles::sam::io::Reader<R>,
    header: &Header,
    reference: &Reference,
    config: &AnalysisConfig,
) -> Result<(ErrorTables, usize)> {
    let start_time = Instant::now();

    let mut tables = ErrorTables::new(config.ont_mode);
    let mut skip_reads = SkipStats::default();
    let mut skip_bases = BaseSkipStats::default();
    let mut read_lengths: Vec<usize> = Vec::with_capacity(100);
    let mut has_tasmanian_tag = false;
    let mut tag_check_count = 0;

    let phred_threshold = config.min_base_quality + 33;

    for result in reader.records() {
        let record = result?;

        // Parse into our record type
        let tas_record = match TasmanianRecord::from_noodles(&record, header) {
            Ok(r) => r,
            Err(e) => {
                debug!("Failed to parse record: {}", e);
                skip_reads.other += 1;
                continue;
            }
        };

        // Check for tasmanian tag in first 10 records
        if tag_check_count < 10 && !has_tasmanian_tag {
            has_tasmanian_tag = tas_record.has_tasmanian_tag();
            tag_check_count += 1;
        }

        // Skip invalid records
        if !tas_record.has_valid_cigar() || tas_record.mapq == 255 {
            skip_reads.other += 1;
            continue;
        }

        if !reference.has_chromosome(&tas_record.chrom) {
            skip_reads.other += 1;
            continue;
        }

        // Filter by indels
        if config.filter_indels && tas_record.has_indels() {
            skip_reads.indel += 1;
            continue;
        }

        // Filter by read length
        let seq_len = tas_record.seq_len();
        if seq_len < config.min_read_length as usize || seq_len > config.max_read_length as usize {
            skip_reads.read_length += 1;
            continue;
        }

        // Filter by mapping quality
        if tas_record.mapq < config.min_mapping_quality {
            skip_reads.map_quality += 1;
            continue;
        }

        // Filter by fragment length
        let abs_tlen = tas_record.abs_tlen();
        if abs_tlen < config.min_fragment_length as u64
            || abs_tlen > config.max_fragment_length as u64
        {
            skip_reads.frag_length += 1;
            continue;
        }

        // Get read number and strand
        let read_num = match tas_record.read_number(config.ont_mode) {
            Some(rn) => rn,
            None => continue,
        };
        let strand = match tas_record.strand(config.ont_mode) {
            Some(s) => s,
            None => continue,
        };

        // Track read lengths for mode calculation
        if read_lengths.len() < 100 && !config.ont_mode {
            read_lengths.push(seq_len);
        } else if config.ont_mode && seq_len > read_lengths.first().copied().unwrap_or(0) {
            if read_lengths.is_empty() {
                read_lengths.push(seq_len);
            } else {
                read_lengths[0] = seq_len;
            }
        }

        // Get reference sequence
        let ref_seq = match reference.get(&tas_record.chrom, tas_record.start, tas_record.end()) {
            Some(s) => s,
            None => {
                debug!(
                    "Reference out of bounds: {}:{}-{}",
                    tas_record.chrom,
                    tas_record.start,
                    tas_record.end()
                );
                continue;
            }
        };

        // Process CIGAR for indels
        let (seq, ref_aligned) = apply_cigar_alignment(&tas_record, ref_seq, config)?;

        if seq.len() != ref_aligned.len() {
            debug!(
                "Sequence/reference length mismatch after CIGAR: {} vs {}",
                seq.len(),
                ref_aligned.len()
            );
            continue;
        }

        // Confidence value from tasmanian tag
        let confidence_value = tas_record.confidence_bases.unwrap_or(0);

        // Process each base
        for pos in 0..seq.len() {
            let base = seq[pos];
            let ref_base = ref_aligned[pos];
            let qual = tas_record.qual.get(pos).copied().unwrap_or(0);

            // Skip N bases and low quality
            if base == b'N' || ref_base == b'N' || qual < phred_threshold {
                skip_bases.quality += 1;
                continue;
            }

            // Skip masked reference bases unless unmasking
            if !config.unmask_genome && ref_base.is_ascii_lowercase() {
                continue;
            }

            // Calculate read position (account for strand)
            let read_pos = match strand {
                Strand::Forward => pos,
                Strand::Reverse => seq_len - pos - 1,
            };

            // Transform bases for reverse strand
            let (final_ref, mut final_base) = match strand {
                Strand::Forward => (to_uppercase(ref_base), to_uppercase(base)),
                Strand::Reverse => (complement(ref_base), complement(base)),
            };

            // Apply methylation masking
            if config.mask_methyl_c || config.mask_methyl_cpg {
                final_base = apply_methylation_mask(
                    final_base,
                    final_ref,
                    read_num,
                    pos,
                    &ref_aligned,
                    config.mask_methyl_cpg,
                );
            }

            // Determine which table to use
            let read_num_val = read_num.as_number();

            if !has_tasmanian_tag || tas_record.tm_tag.is_none() {
                // No intersection info - use unrelated table
                tables
                    .unrelated
                    .increment(read_num_val, read_pos, final_ref, final_base);
            } else if tas_record.is_in_intersection(pos) {
                // Base is in intersection region
                if confidence_value >= config.confidence_threshold {
                    tables
                        .intersection_confident
                        .increment(read_num_val, read_pos, final_ref, final_base);
                } else {
                    tables
                        .intersection
                        .increment(read_num_val, read_pos, final_ref, final_base);
                }
            } else {
                // Base is in complement region
                if confidence_value >= config.confidence_threshold {
                    tables
                        .complement_confident
                        .increment(read_num_val, read_pos, final_ref, final_base);
                } else {
                    tables
                        .complement
                        .increment(read_num_val, read_pos, final_ref, final_base);
                }
            }
        }
    }

    // Calculate mode read length
    let mode_length = if config.ont_mode {
        read_lengths.first().copied().unwrap_or(1000)
    } else {
        calculate_mode(&read_lengths).unwrap_or(150)
    };

    // Trim tables to mode length
    tables.trim_all(mode_length);

    let elapsed = start_time.elapsed();

    // Report statistics to stderr
    eprintln!("Analysis completed in {:.2} seconds", elapsed.as_secs_f64());
    eprintln!("Reads skipped:");
    eprintln!("  indel:       {}", skip_reads.indel);
    eprintln!("  mapquality:  {}", skip_reads.map_quality);
    eprintln!("  fragLength:  {}", skip_reads.frag_length);
    eprintln!("  readlength:  {}", skip_reads.read_length);
    eprintln!("  other:       {}", skip_reads.other);
    eprintln!("Bases skipped:");
    eprintln!("  quality:     {}", skip_bases.quality);
    eprintln!("  softclips:   {}", skip_bases.softclips);

    Ok((tables, mode_length))
}

/// Apply CIGAR alignment to get aligned sequence and reference
fn apply_cigar_alignment(
    record: &TasmanianRecord,
    ref_seq: &[u8],
    config: &AnalysisConfig,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let cigar_ops = parse_cigar(&record.cigar_str);

    if cigar_ops.is_empty() {
        return Ok((record.seq.clone(), ref_seq.to_vec()));
    }

    let mut seq_out = Vec::with_capacity(record.seq.len());
    let mut ref_out = Vec::with_capacity(record.seq.len());

    let mut seq_pos = 0usize;
    let mut ref_pos = 0usize;

    for (op, len) in cigar_ops {
        let len = len as usize;

        use noodles::sam::alignment::record::cigar::op::Kind;
        match op {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Both sequence and reference advance
                for _ in 0..len {
                    if seq_pos < record.seq.len() && ref_pos < ref_seq.len() {
                        seq_out.push(record.seq[seq_pos]);
                        ref_out.push(ref_seq[ref_pos]);
                    }
                    seq_pos += 1;
                    ref_pos += 1;
                }
            }
            Kind::Insertion => {
                // Sequence advances, reference doesn't
                // Skip inserted bases (not in reference)
                seq_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                // Reference advances, sequence doesn't
                // Skip deleted positions
                ref_pos += len;
            }
            Kind::SoftClip => {
                // Handle softclips based on bypass mode
                match config.softclip_bypass {
                    1 => {
                        // Skip softclipped bases
                        seq_pos += len;
                    }
                    2 => {
                        // Force use softclipped bases
                        for _ in 0..len {
                            if seq_pos < record.seq.len() && ref_pos < ref_seq.len() {
                                seq_out.push(record.seq[seq_pos]);
                                ref_out.push(ref_seq[ref_pos]);
                            }
                            seq_pos += 1;
                            ref_pos += 1;
                        }
                    }
                    _ => {
                        // Mode 0: check if garbage
                        // For simplicity, skip them for now (TODO: implement garbage check)
                        seq_pos += len;
                    }
                }
            }
            Kind::HardClip | Kind::Pad => {
                // Don't consume sequence or reference
            }
        }
    }

    Ok((seq_out, ref_out))
}

/// Apply methylation masking to a base
fn apply_methylation_mask(
    base: u8,
    ref_base: u8,
    read_num: ReadNumber,
    pos: usize,
    ref_seq: &[u8],
    cpg_only: bool,
) -> u8 {
    if cpg_only {
        // Only mask at CpG sites
        let is_cpg = match read_num {
            ReadNumber::First => {
                ref_base == b'C' && ref_seq.get(pos + 1).map_or(false, |&b| b == b'G' || b == b'g')
            }
            ReadNumber::Second => {
                ref_base == b'G' && ref_seq.get(pos + 1).map_or(false, |&b| b == b'C' || b == b'c')
            }
        };

        if is_cpg {
            match (read_num, ref_base, base) {
                (ReadNumber::First, b'C', b'T') => b'C',
                (ReadNumber::Second, b'G', b'A') => b'G',
                _ => base,
            }
        } else {
            base
        }
    } else {
        // Mask all C->T (R1) and G->A (R2)
        match (read_num, ref_base, base) {
            (ReadNumber::First, b'C', b'T') => b'C',
            (ReadNumber::Second, b'G', b'A') => b'G',
            _ => base,
        }
    }
}

/// Calculate the mode of a vector of values
fn calculate_mode(values: &[usize]) -> Option<usize> {
    if values.is_empty() {
        return None;
    }

    let mut counts: HashMap<usize, usize> = HashMap::new();
    for &v in values {
        *counts.entry(v).or_insert(0) += 1;
    }

    counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(value, _)| value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_mode() {
        assert_eq!(calculate_mode(&[100, 100, 100, 150, 150]), Some(100));
        assert_eq!(calculate_mode(&[]), None);
        assert_eq!(calculate_mode(&[42]), Some(42));
    }

    #[test]
    fn test_methylation_mask() {
        // C->T on read 1 should become C
        assert_eq!(
            apply_methylation_mask(b'T', b'C', ReadNumber::First, 0, &[], false),
            b'C'
        );

        // G->A on read 2 should become G
        assert_eq!(
            apply_methylation_mask(b'A', b'G', ReadNumber::Second, 0, &[], false),
            b'G'
        );

        // Other mismatches unchanged
        assert_eq!(
            apply_methylation_mask(b'G', b'A', ReadNumber::First, 0, &[], false),
            b'G'
        );
    }
}
