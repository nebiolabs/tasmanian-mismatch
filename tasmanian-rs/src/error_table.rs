//! Mismatch counting error tables
//!
//! Tracks mismatches between sequenced bases and reference bases,
//! organized by read number and position.

use crate::reference::{base_to_index, index_to_base};

/// Maximum read length to pre-allocate
const DEFAULT_MAX_LENGTH: usize = 1000;
const ONT_MAX_LENGTH: usize = 100000;

/// A single error table tracking mismatches
///
/// Structure: [read_number][position][ref_base][alt_base] = count
/// - read_number: 0 (read1) or 1 (read2)
/// - position: 0-based position in read
/// - ref_base: 0-3 for A,C,G,T
/// - alt_base: 0-3 for A,C,G,T
#[derive(Debug, Clone)]
pub struct ErrorTable {
    /// counts[read][position][ref_base][alt_base]
    counts: Vec<Vec<[[u64; 4]; 4]>>,
    max_length: usize,
}

impl ErrorTable {
    /// Create a new error table with specified max length
    pub fn new(max_length: usize) -> Self {
        let counts = vec![vec![[[0u64; 4]; 4]; max_length]; 2];
        Self { counts, max_length }
    }

    /// Create a standard error table (max length 1000)
    pub fn standard() -> Self {
        Self::new(DEFAULT_MAX_LENGTH)
    }

    /// Create an ONT error table (max length 100000)
    pub fn ont() -> Self {
        Self::new(ONT_MAX_LENGTH)
    }

    /// Increment the count for a specific mismatch
    ///
    /// - read_num: 1 or 2
    /// - position: 0-based position in read
    /// - ref_base: reference base (A/C/G/T)
    /// - alt_base: sequenced base (A/C/G/T)
    pub fn increment(&mut self, read_num: u8, position: usize, ref_base: u8, alt_base: u8) {
        let read_idx = (read_num - 1) as usize;
        if read_idx >= 2 || position >= self.max_length {
            return;
        }

        if let (Some(ref_idx), Some(alt_idx)) = (base_to_index(ref_base), base_to_index(alt_base)) {
            self.counts[read_idx][position][ref_idx][alt_idx] += 1;
        }
    }

    /// Get the count for a specific mismatch
    pub fn get(&self, read_num: u8, position: usize, ref_base: u8, alt_base: u8) -> u64 {
        let read_idx = (read_num - 1) as usize;
        if read_idx >= 2 || position >= self.max_length {
            return 0;
        }

        if let (Some(ref_idx), Some(alt_idx)) = (base_to_index(ref_base), base_to_index(alt_base)) {
            self.counts[read_idx][position][ref_idx][alt_idx]
        } else {
            0
        }
    }

    /// Trim the table to a specific length
    pub fn trim(&mut self, new_length: usize) {
        if new_length < self.max_length {
            for read_counts in &mut self.counts {
                read_counts.truncate(new_length);
            }
            self.max_length = new_length;
        }
    }

    /// Get the current max length
    pub fn max_length(&self) -> usize {
        self.max_length
    }

    /// Get counts for a specific read and position
    pub fn get_position_counts(&self, read_num: u8, position: usize) -> Option<&[[u64; 4]; 4]> {
        let read_idx = (read_num - 1) as usize;
        if read_idx >= 2 || position >= self.max_length {
            return None;
        }
        Some(&self.counts[read_idx][position])
    }

    /// Iterate over all positions for a read
    pub fn iter_read(&self, read_num: u8) -> impl Iterator<Item = (usize, &[[u64; 4]; 4])> {
        let read_idx = (read_num - 1) as usize;
        self.counts
            .get(read_idx)
            .into_iter()
            .flat_map(|v| v.iter().enumerate())
    }
}

/// Collection of error tables for different categories
#[derive(Debug)]
pub struct ErrorTables {
    /// Bases overlapping BED regions (all)
    pub intersection: ErrorTable,
    /// Non-overlapping bases from intersecting reads (all)
    pub complement: ErrorTable,
    /// All bases from reads with no BED intersection
    pub unrelated: ErrorTable,
    /// High-confidence intersection bases
    pub intersection_confident: ErrorTable,
    /// High-confidence complement bases
    pub complement_confident: ErrorTable,
}

impl ErrorTables {
    /// Create a new set of error tables
    pub fn new(ont_mode: bool) -> Self {
        let create = if ont_mode {
            ErrorTable::ont
        } else {
            ErrorTable::standard
        };

        Self {
            intersection: create(),
            complement: create(),
            unrelated: create(),
            intersection_confident: create(),
            complement_confident: create(),
        }
    }

    /// Trim all tables to a specific length
    pub fn trim_all(&mut self, length: usize) {
        self.intersection.trim(length);
        self.complement.trim(length);
        self.unrelated.trim(length);
        self.intersection_confident.trim(length);
        self.complement_confident.trim(length);
    }

    /// Get the current max length (from any table)
    pub fn max_length(&self) -> usize {
        self.unrelated.max_length()
    }
}

/// Iterator over mismatch types (for output)
pub struct MismatchTypeIter {
    ref_idx: usize,
    alt_idx: usize,
}

impl MismatchTypeIter {
    pub fn new() -> Self {
        Self {
            ref_idx: 0,
            alt_idx: 0,
        }
    }
}

impl Default for MismatchTypeIter {
    fn default() -> Self {
        Self::new()
    }
}

impl Iterator for MismatchTypeIter {
    type Item = (u8, u8); // (ref_base, alt_base)

    fn next(&mut self) -> Option<Self::Item> {
        if self.ref_idx >= 4 {
            return None;
        }

        let ref_base = index_to_base(self.ref_idx)?;
        let alt_base = index_to_base(self.alt_idx)?;

        self.alt_idx += 1;
        if self.alt_idx >= 4 {
            self.alt_idx = 0;
            self.ref_idx += 1;
        }

        Some((ref_base, alt_base))
    }
}

/// Get the column header for a mismatch type
pub fn mismatch_column_name(prefix: &str, ref_base: u8, alt_base: u8) -> String {
    format!(
        "{}{}_{}",
        prefix,
        (ref_base as char).to_ascii_lowercase(),
        (alt_base as char).to_ascii_lowercase()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_table_increment() {
        let mut table = ErrorTable::standard();
        table.increment(1, 0, b'A', b'T');
        table.increment(1, 0, b'A', b'T');
        table.increment(1, 0, b'C', b'G');

        assert_eq!(table.get(1, 0, b'A', b'T'), 2);
        assert_eq!(table.get(1, 0, b'C', b'G'), 1);
        assert_eq!(table.get(1, 0, b'A', b'A'), 0);
    }

    #[test]
    fn test_mismatch_iter() {
        let iter = MismatchTypeIter::new();
        let types: Vec<_> = iter.collect();

        assert_eq!(types.len(), 16);
        assert_eq!(types[0], (b'A', b'A'));
        assert_eq!(types[1], (b'A', b'C'));
        assert_eq!(types[4], (b'C', b'A'));
        assert_eq!(types[15], (b'T', b'T'));
    }

    #[test]
    fn test_trim() {
        let mut table = ErrorTable::standard();
        assert_eq!(table.max_length(), 1000);
        table.trim(100);
        assert_eq!(table.max_length(), 100);
    }
}
