//! BED file parsing and region indexing

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A single BED region with optional metadata
#[derive(Debug, Clone)]
pub struct BedRegion {
    pub start: u64,
    pub end: u64,
    pub strand: Option<char>,
    pub name: Option<String>,
    pub class: Option<String>,
    pub family: Option<String>,
}

/// Index of BED regions organized by chromosome
///
/// Regions are sorted by start position within each chromosome
/// for efficient sequential scanning.
#[derive(Debug)]
pub struct BedIndex {
    regions: HashMap<String, Vec<BedRegion>>,
    total_regions: usize,
}

impl BedIndex {
    /// Load a BED file and create an index
    ///
    /// Supports 3-column (chrom, start, end) and 7+ column formats
    /// (chrom, start, end, strand, name, class, family)
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open BED file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut regions: HashMap<String, Vec<BedRegion>> = HashMap::new();
        let mut total_regions = 0;

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result.with_context(|| format!("Failed to read line {}", line_num + 1))?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();

            // Need at least 3 columns
            if fields.len() < 3 {
                continue;
            }

            let chrom = fields[0].to_string();

            // Parse start and end, skip header lines
            let start = match fields[1].parse::<u64>() {
                Ok(v) => v,
                Err(_) => continue, // Likely header row
            };
            let end = match fields[2].parse::<u64>() {
                Ok(v) => v + 1, // Include upper bound
                Err(_) => continue,
            };

            // Parse optional fields
            let (strand, name, class, family) = if fields.len() >= 7 {
                (
                    fields[3].chars().next(),
                    Some(fields[4].to_string()),
                    Some(fields[5].to_string()),
                    Some(fields[6].to_string()),
                )
            } else {
                (None, None, None, None)
            };

            let region = BedRegion {
                start,
                end,
                strand,
                name,
                class,
                family,
            };

            regions.entry(chrom).or_default().push(region);
            total_regions += 1;
        }

        // Sort regions by start position within each chromosome
        for regions in regions.values_mut() {
            regions.sort_by_key(|r| r.start);
        }

        Ok(BedIndex {
            regions,
            total_regions,
        })
    }

    /// Get regions for a chromosome
    pub fn get_regions(&self, chrom: &str) -> Option<&[BedRegion]> {
        self.regions.get(chrom).map(|v| v.as_slice())
    }

    /// Check if a chromosome has any regions
    pub fn has_chromosome(&self, chrom: &str) -> bool {
        self.regions.contains_key(chrom)
    }

    /// Get total number of regions
    pub fn total_regions(&self) -> usize {
        self.total_regions
    }

    /// Get number of chromosomes with regions
    pub fn num_chromosomes(&self) -> usize {
        self.regions.len()
    }

    /// Get the number of regions for a specific chromosome
    pub fn num_regions(&self, chrom: &str) -> usize {
        self.regions.get(chrom).map_or(0, |v| v.len())
    }
}

/// Tracks position in BED file for efficient sequential scanning
#[derive(Debug, Default)]
pub struct BedScanner {
    current_chrom: String,
    current_index: usize,
    finished: bool,
}

impl BedScanner {
    pub fn new() -> Self {
        Self::default()
    }

    /// Find intersection with BED regions for a given genomic interval
    ///
    /// Returns (region_index, intersection_start, intersection_end) if found
    pub fn find_intersection(
        &mut self,
        bed: &BedIndex,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Option<(usize, u64, u64)> {
        // Handle chromosome change
        if chrom != self.current_chrom {
            self.current_chrom = chrom.to_string();
            self.current_index = 0;
            self.finished = false;
        }

        if self.finished {
            return None;
        }

        let regions = bed.get_regions(chrom)?;

        // Advance past regions that end before our start
        while self.current_index < regions.len() && regions[self.current_index].end <= start {
            self.current_index += 1;
        }

        if self.current_index >= regions.len() {
            self.finished = true;
            return None;
        }

        let region = &regions[self.current_index];

        // Check if read ends before region starts (no intersection)
        if end <= region.start {
            return None;
        }

        // We have an intersection
        let intersection_start = start.max(region.start);
        let intersection_end = end.min(region.end);

        Some((self.current_index, intersection_start, intersection_end))
    }

    /// Reset scanner for new chromosome
    pub fn reset(&mut self) {
        self.current_chrom.clear();
        self.current_index = 0;
        self.finished = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_temp_bed(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_load_bed_3col() {
        let content = "chr1\t100\t200\nchr1\t300\t400\nchr2\t500\t600\n";
        let file = create_temp_bed(content);
        let index = BedIndex::load(file.path()).unwrap();

        assert_eq!(index.total_regions(), 3);
        assert_eq!(index.num_chromosomes(), 2);
        assert_eq!(index.num_regions("chr1"), 2);
        assert_eq!(index.num_regions("chr2"), 1);
    }

    #[test]
    fn test_load_bed_7col() {
        let content = "chr1\t100\t200\t+\tL1P5\tLINE\tL1\n";
        let file = create_temp_bed(content);
        let index = BedIndex::load(file.path()).unwrap();

        let regions = index.get_regions("chr1").unwrap();
        assert_eq!(regions[0].name.as_deref(), Some("L1P5"));
        assert_eq!(regions[0].class.as_deref(), Some("LINE"));
    }

    #[test]
    fn test_scanner_intersection() {
        let content = "chr1\t100\t200\nchr1\t300\t400\n";
        let file = create_temp_bed(content);
        let index = BedIndex::load(file.path()).unwrap();

        let mut scanner = BedScanner::new();

        // No intersection (before first region)
        assert!(scanner.find_intersection(&index, "chr1", 0, 50).is_none());

        // Intersection with first region
        let result = scanner.find_intersection(&index, "chr1", 150, 250);
        assert!(result.is_some());
        let (idx, start, end) = result.unwrap();
        assert_eq!(idx, 0);
        assert_eq!(start, 150);
        assert_eq!(end, 201); // end is inclusive in BedRegion

        // Between regions
        assert!(scanner
            .find_intersection(&index, "chr1", 220, 280)
            .is_none());
    }
}
