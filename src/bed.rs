// BED file handling for filtering/masking genomic regions
//
// This module provides functionality to parse BED files and check for
// overlaps between reads and BED regions.

use crate::types::ReferenceGenome;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Represents a single interval in a BED file
#[derive(Debug, Clone)]
pub struct BedInterval {
    pub start: i64,
    pub end: i64,
}

/// Container for BED intervals organized by chromosome
pub type BedRegions = HashMap<String, Vec<BedInterval>>;

/// Parse a BED file and return intervals organized by chromosome
///
/// # Arguments
/// * `bed_path` - Path to the BED file
///
/// # Returns
/// * A HashMap where keys are chromosome names and values are vectors of intervals
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

/// Fast position check using a pre-filtered list of intervals
///
/// # Arguments
/// * `intervals` - Pre-filtered intervals for the current genomic region
/// * `pos` - Genomic position (0-based)
///
/// # Returns
/// * `true` if the position overlaps with any interval, `false` otherwise
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

/// Filter BED intervals to only those overlapping a specific genomic region
///
/// # Arguments
/// * `bed_regions` - The BED regions to filter
/// * `chrom` - Chromosome name
/// * `region_start` - Start position of the region (0-based, inclusive)
/// * `region_end` - End position of the region (0-based, exclusive)
///
/// # Returns
/// * A vector of BED intervals that overlap with the specified region
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

/// Mask regions in the reference genome by converting bases to 'N'
/// This is much more efficient than checking BED overlap for every read position
///
/// # Arguments
/// * `reference` - Mutable reference to the reference genome
/// * `bed_regions` - BED regions to mask
///
/// # Returns
/// * Number of bases masked
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
