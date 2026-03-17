//! BED file handling for filtering and masking genomic regions.
//!
//! This module provides helpers to parse sorted BED files, query interval
//! overlaps, and mask reference bases that fall inside the supplied regions.

use crate::types::ReferenceGenome;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Represents a single interval in a BED file.
#[derive(Debug, Clone)]
pub struct BedInterval {
    /// Inclusive 0-based interval start.
    pub start: i64,
    /// Exclusive 0-based interval end.
    pub end: i64,
}

/// BED intervals grouped by chromosome or contig name.
pub type BedRegions = HashMap<String, Vec<BedInterval>>;

/// Parse a BED file into chromosome-indexed intervals.
///
/// # Arguments
/// * `bed_path` - Path to the BED file.
///
/// # Returns
/// * A [`BedRegions`] map where each key is a chromosome name and each value is
///   the corresponding list of intervals.
///
/// # Errors
/// * If the BED file cannot be opened.
/// * If a start or end coordinate cannot be parsed.
///
/// # Panics
/// * If intervals for a chromosome are not sorted by start position.
pub fn parse_bed_file(bed_path: &str) -> Result<BedRegions, Box<dyn std::error::Error>> {
    let file = File::open(bed_path)?;
    let reader = BufReader::new(file);

    let mut regions: BedRegions = HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;

        // Skip empty lines and comments
        if line.trim().is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 3 {
            eprintln!(
                "Warning: Skipping malformed line {} in BED file (need at least 3 fields)",
                line_num + 1
            );
            continue;
        }

        let chrom = fields[0].to_string();
        let start: i64 = fields[1]
            .parse()
            .map_err(|_| format!("Invalid start coordinate at line {}", line_num + 1))?;
        let end: i64 = fields[2]
            .parse()
            .map_err(|_| format!("Invalid end coordinate at line {}", line_num + 1))?;

        if start >= end {
            eprintln!(
                "Warning: Skipping invalid interval at line {} (start >= end)",
                line_num + 1
            );
            continue;
        }

        regions
            .entry(chrom)
            .or_insert_with(Vec::new)
            .push(BedInterval { start, end });
    }

    // Validate that intervals are sorted by start position within each chromosome
    for (chrom, intervals) in &regions {
        for i in 1..intervals.len() {
            if intervals[i].start < intervals[i - 1].start {
                panic!(
                    "BED file is not sorted! Chromosome '{}' has unsorted intervals: \
                    interval at position {} (start={}) comes before interval at position {} (start={}). \
                    Please sort your BED file using 'sort -k1,1 -k2,2n'",
                    chrom, i - 1, intervals[i - 1].start, i, intervals[i].start
                );
            }
        }
    }

    let total_intervals: usize = regions.values().map(|v| v.len()).sum();
    eprintln!(
        "Loaded {} BED intervals across {} chromosomes",
        total_intervals,
        regions.len()
    );

    Ok(regions)
}

/// Check whether a genomic position overlaps any interval in a sorted slice.
///
/// # Arguments
/// * `intervals` - Pre-filtered intervals for the current genomic region.
/// * `pos` - Genomic position in 0-based coordinates.
///
/// # Returns
/// * `true` if the position overlaps an interval.
/// * `false` otherwise.
#[inline]
pub fn position_overlaps_intervals(intervals: &[BedInterval], pos: i64) -> bool {
    for interval in intervals {
        if pos < interval.start {
            break; // Early exit since sorted
        }
        if pos >= interval.start && pos < interval.end {
            return true;
        }
    }
    false
}

/// Return BED intervals that overlap a genomic region.
///
/// # Arguments
/// * `bed_regions` - BED intervals grouped by chromosome.
/// * `chrom` - Chromosome or contig name.
/// * `region_start` - Region start in 0-based inclusive coordinates.
/// * `region_end` - Region end in 0-based exclusive coordinates.
///
/// # Returns
/// * A vector of [`BedInterval`] values that overlap the requested region.
pub fn filter_bed_for_region(
    bed_regions: &BedRegions,
    chrom: &str,
    region_start: i64,
    region_end: i64,
) -> Vec<BedInterval> {
    if let Some(intervals) = bed_regions.get(chrom) {
        intervals
            .iter()
            .filter(|interval| {
                // Check if interval overlaps with region
                interval.start < region_end && interval.end > region_start
            })
            .cloned()
            .collect()
    } else {
        Vec::new()
    }
}

/// Mask BED-covered bases in the reference genome with `N`.
///
/// This is more efficient than checking BED overlap for every aligned base.
///
/// # Arguments
/// * `reference` - Mutable reference genome to update in place.
/// * `bed_regions` - BED intervals that should be masked.
///
/// # Returns
/// * The number of reference bases replaced with `N`.
pub fn mask_reference_with_bed(reference: &mut ReferenceGenome, bed_regions: &BedRegions) -> usize {
    let mut total_masked = 0;

    for (chrom, intervals) in bed_regions.iter() {
        if let Some(ref_seq) = reference.get_mut(chrom) {
            for interval in intervals {
                let start = interval.start.max(0) as usize;
                let end = (interval.end as usize).min(ref_seq.len());

                for pos in start..end {
                    ref_seq[pos] = b'N'; // Mask with 'N'
                    total_masked += 1;
                }
            }
        }
    }

    total_masked
}
