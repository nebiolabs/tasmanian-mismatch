//! Intersection detection between reads and BED regions
//!
//! This module processes SAM records and detects overlaps with BED regions,
//! adding Tasmanian-specific tags to the output.
//!
//! NOTE: This is a simplified implementation for the MVP. Full intersection
//! processing with SAM output will be implemented in a future version.

use crate::bed::{BedIndex, BedScanner};
use anyhow::Result;
use noodles::sam::{
    alignment::Record,
    Header,
};
use std::collections::HashMap;
use std::io::{BufRead, Write};

/// Intersection category codes
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    /// No intersection with BED regions
    Unrelated = 1,
    /// Only read 1 intersects
    Read1Only = 2,
    /// Only read 2 intersects
    Read2Only = 3,
    /// Both reads intersect same region
    BothSame = 4,
    /// Both reads intersect different regions
    BothDifferent = 5,
}

/// Subcategory for intersection position
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubCategory {
    /// Completely contained
    Complete,
    /// Partial left (read extends before region)
    Left,
    /// Partial right (read extends after region)
    Right,
    /// Both sides (region contained in read)
    Both,
}

impl SubCategory {
    pub fn as_char(&self) -> char {
        match self {
            SubCategory::Complete => 'a',
            SubCategory::Left => 'b',
            SubCategory::Right => 'c',
            SubCategory::Both => 'd',
        }
    }

    /// Determine subcategory from intersection bounds
    pub fn from_bounds(read_start: u64, read_end: u64, int_start: u64, int_end: u64) -> Self {
        let extends_left = read_start < int_start;
        let extends_right = read_end > int_end;

        match (extends_left, extends_right) {
            (false, false) => SubCategory::Complete,
            (true, false) => SubCategory::Left,
            (false, true) => SubCategory::Right,
            (true, true) => SubCategory::Both,
        }
    }
}

/// Result of intersection detection for a single read
#[derive(Debug, Clone)]
pub struct IntersectionResult {
    /// Intersection bounds within read (0-based)
    pub bounds: (usize, usize),
    /// BED region index
    pub bed_index: usize,
    /// Subcategory
    pub subcategory: SubCategory,
}

/// Simplified intersection processor that works with raw SAM lines
///
/// This is a temporary implementation that processes SAM text directly
/// rather than using noodles RecordBuf (which has complex API requirements).
pub struct IntersectionProcessor {
    bed_index: BedIndex,
    scanner: BedScanner,
    /// Buffer for unpaired reads, keyed by read name
    read_buffer: HashMap<String, (String, Option<IntersectionResult>)>,
}

impl IntersectionProcessor {
    /// Create a new intersection processor
    pub fn new(bed_index: BedIndex) -> Self {
        Self {
            bed_index,
            scanner: BedScanner::new(),
            read_buffer: HashMap::new(),
        }
    }

    /// Process SAM input from a reader and write annotated output
    ///
    /// This is a line-based implementation that avoids noodles RecordBuf complexity.
    pub fn process_stream<R: BufRead, W: Write>(
        &mut self,
        reader: R,
        mut writer: W,
    ) -> Result<()> {
        for line_result in reader.lines() {
            let line = line_result?;

            // Pass through header lines
            if line.starts_with('@') {
                writeln!(writer, "{}", line)?;
                continue;
            }

            // Parse SAM line
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 11 {
                continue;
            }

            let read_name = fields[0];
            let flag: u16 = fields[1].parse().unwrap_or(0);
            let chrom = fields[2];
            let start: u64 = fields[3].parse::<u64>().unwrap_or(1) - 1; // Convert to 0-based
            let seq = fields[9];
            let seq_len = seq.len() as u64;
            let end = start + seq_len;

            // Only process proper pairs
            if !is_proper_pair(flag) {
                // Pass through unchanged with default tags
                writeln!(writer, "{}\ttm:Z:-1;\ttc:i:0", line)?;
                continue;
            }

            // Check for intersection
            let intersection = if self.bed_index.has_chromosome(chrom) {
                self.scanner
                    .find_intersection(&self.bed_index, chrom, start, end)
                    .map(|(bed_idx, int_start, int_end)| {
                        let read_int_start = (int_start.saturating_sub(start)) as usize;
                        let read_int_end = (int_end - start).min(seq_len) as usize;

                        IntersectionResult {
                            bounds: (read_int_start, read_int_end),
                            bed_index: bed_idx,
                            subcategory: SubCategory::from_bounds(start, end, int_start, int_end),
                        }
                    })
            } else {
                None
            };

            // Check if we have the mate
            if let Some((mate_line, mate_int)) = self.read_buffer.remove(read_name) {
                // Calculate complement bases
                let mate_fields: Vec<&str> = mate_line.split('\t').collect();
                let mate_seq_len = mate_fields.get(9).map_or(0, |s| s.len());

                let complement1 = if let Some(ref int) = mate_int {
                    mate_seq_len - (int.bounds.1 - int.bounds.0)
                } else {
                    mate_seq_len
                };

                let complement2 = if let Some(ref int) = intersection {
                    seq_len as usize - (int.bounds.1 - int.bounds.0)
                } else {
                    seq_len as usize
                };

                let total_complement = complement1 + complement2;

                // Output mate with tags
                let mate_junction = mate_int
                    .as_ref()
                    .map(|i| format!("{}.{};", i.bounds.0, i.bounds.1))
                    .unwrap_or_else(|| "-1;".to_string());

                // Apply lowercase masking to intersection region
                let mate_masked_seq = if let Some(ref int) = mate_int {
                    mask_sequence(mate_fields[9], int.bounds.0, int.bounds.1)
                } else {
                    mate_fields[9].to_string()
                };

                // Reconstruct mate line with masked sequence
                let mut mate_out_fields: Vec<String> = mate_fields.iter().map(|s| s.to_string()).collect();
                mate_out_fields[9] = mate_masked_seq;
                writeln!(
                    writer,
                    "{}\ttm:Z:{}\ttc:i:{}",
                    mate_out_fields.join("\t"),
                    mate_junction,
                    total_complement
                )?;

                // Output current with tags
                let current_junction = intersection
                    .as_ref()
                    .map(|i| format!("{}.{};", i.bounds.0, i.bounds.1))
                    .unwrap_or_else(|| "-1;".to_string());

                let current_masked_seq = if let Some(ref int) = intersection {
                    mask_sequence(seq, int.bounds.0, int.bounds.1)
                } else {
                    seq.to_string()
                };

                // Reconstruct current line with masked sequence
                let mut current_fields_owned: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
                current_fields_owned[9] = current_masked_seq;
                writeln!(
                    writer,
                    "{}\ttm:Z:{}\ttc:i:{}",
                    current_fields_owned.join("\t"),
                    current_junction,
                    total_complement
                )?;
            } else {
                // Buffer this read
                self.read_buffer
                    .insert(read_name.to_string(), (line.clone(), intersection));
            }
        }

        // Flush remaining unpaired reads
        for (_, (line, intersection)) in self.read_buffer.drain() {
            let junction = intersection
                .as_ref()
                .map(|i| format!("{}.{};", i.bounds.0, i.bounds.1))
                .unwrap_or_else(|| "-1;".to_string());

            let fields: Vec<&str> = line.split('\t').collect();
            let masked_seq = if let Some(ref int) = intersection {
                mask_sequence(fields.get(9).unwrap_or(&""), int.bounds.0, int.bounds.1)
            } else {
                fields.get(9).unwrap_or(&"").to_string()
            };

            let mut out_fields: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
            if out_fields.len() > 9 {
                out_fields[9] = masked_seq;
            }
            writeln!(writer, "{}\ttm:Z:{}\ttc:i:0", out_fields.join("\t"), junction)?;
        }

        Ok(())
    }

    /// Placeholder for noodles-based processing (not fully implemented)
    pub fn process<R: Record>(&mut self, _record: &R, _header: &Header) -> Result<Option<Vec<u8>>> {
        // Not implemented - use process_stream instead
        unimplemented!("Use process_stream for intersection processing")
    }

    /// Placeholder for flush (not needed with process_stream)
    pub fn flush(&mut self, _header: &Header) -> Result<Vec<Vec<u8>>> {
        Ok(Vec::new())
    }
}

/// Mask sequence at intersection region (convert to lowercase)
fn mask_sequence(seq: &str, start: usize, end: usize) -> String {
    let mut result = String::with_capacity(seq.len());
    for (i, c) in seq.chars().enumerate() {
        if i >= start && i < end {
            result.push(c.to_ascii_lowercase());
        } else {
            result.push(c);
        }
    }
    result
}

/// Check if a flag represents a proper pair
fn is_proper_pair(flag: u16) -> bool {
    matches!(flag, 99 | 163 | 83 | 147)
}

/// Check if this is read 1
fn is_read1(flag: u16) -> bool {
    matches!(flag, 99 | 83)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subcategory() {
        // Read completely contained in region
        assert_eq!(
            SubCategory::from_bounds(100, 200, 50, 250),
            SubCategory::Complete
        );

        // Read extends left
        assert_eq!(
            SubCategory::from_bounds(50, 150, 100, 200),
            SubCategory::Left
        );

        // Read extends right
        assert_eq!(
            SubCategory::from_bounds(150, 250, 100, 200),
            SubCategory::Right
        );

        // Read contains region
        assert_eq!(
            SubCategory::from_bounds(50, 250, 100, 200),
            SubCategory::Both
        );
    }

    #[test]
    fn test_mask_sequence() {
        assert_eq!(mask_sequence("ACGTACGT", 2, 6), "ACgtacGT");
        assert_eq!(mask_sequence("ACGT", 0, 4), "acgt");
        assert_eq!(mask_sequence("ACGT", 0, 0), "ACGT");
    }
}
