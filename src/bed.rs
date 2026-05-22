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

/// Bundles BED filtering state for a processing run.
#[derive(Clone, Copy)]
pub struct BedFilter<'a> {
    /// Per-chromosome BED intervals, or `None` when no BED file was supplied.
    pub regions: Option<&'a BedRegions>,
    /// When `true`, skip any read whose alignment overlaps a BED interval.
    /// When `false`, individual bases at BED positions are masked in the reference.
    pub filter_whole_reads: bool,
}

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
            log::warn!(
                "Skipping malformed line {} in BED file (need at least 3 fields)",
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
            log::warn!(
                "Skipping invalid interval at line {} (start >= end)",
                line_num + 1
            );
            continue;
        }

        regions
            .entry(chrom)
            .or_default()
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

    // Merge overlapping, nested, or adjacent intervals to enforce non-overlapping invariant.
    // This allows binary search in position_overlaps_intervals to be correct and fast.
    for intervals in regions.values_mut() {
        if intervals.is_empty() {
            continue;
        }
        let mut merged = Vec::with_capacity(intervals.len());
        merged.push(intervals[0].clone());
        for next in intervals.iter().skip(1) {
            let last = merged.last_mut().unwrap();
            if next.start <= last.end {
                // Overlapping or adjacent, merge them
                last.end = last.end.max(next.end);
            } else {
                merged.push(next.clone());
            }
        }
        *intervals = merged;
    }

    let total_intervals: usize = regions.values().map(|v| v.len()).sum();
    log::info!(
        "Loaded {} BED intervals across {} chromosomes",
        total_intervals,
        regions.len()
    );

    Ok(regions)
}

/// Optionally parse a BED file, returning `None` when no path is provided.
///
/// This is a thin convenience wrapper around [`parse_bed_file`] for call sites
/// that hold an `Option<&str>` path (e.g. a CLI flag that may not be set).
///
/// # Panics
/// * Panics with a descriptive message if a path is given but the file cannot
///   be opened or parsed.
pub fn maybe_parse_bed_file(path: Option<&str>) -> Option<BedRegions> {
    let path = path?;
    Some(
        parse_bed_file(path)
            .expect("Failed to load BED file. Please check the file path and format."),
    )
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
    // Binary search for the rightmost interval whose start <= pos, then
    // check whether pos falls before its end.  Assumes non-overlapping
    // intervals sorted by start (the invariant enforced by parse_bed_file).
    let idx = intervals.partition_point(|i| i.start <= pos);
    idx > 0 && pos < intervals[idx - 1].end
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

                if start >= end {
                    continue;
                }

                ref_seq[start..end].fill(b'N');
                total_masked += end - start;
            }
        }
    }

    total_masked
}
