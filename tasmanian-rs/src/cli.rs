//! CLI argument definitions using clap

use clap::{Parser, Subcommand};
use std::path::PathBuf;

/// Tasmanian - Analysis of reference mismatches in high throughput sequencing data
///
/// Unlike other tools, Tasmanian can evaluate portions of reads that overlap with
/// specified genomic regions (e.g., repeats) separately from non-overlapping portions.
#[derive(Parser, Debug)]
#[command(name = "tasmanian")]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Process intersections with BED regions
    ///
    /// Reads SAM from stdin, identifies intersections with BED regions,
    /// and outputs annotated SAM to stdout.
    Intersections(IntersectionsArgs),

    /// Analyze artifacts in sequencing data
    ///
    /// Reads SAM from stdin, analyzes mismatches against reference,
    /// and outputs CSV table to stdout.
    Analyze(AnalyzeArgs),
}

#[derive(Parser, Debug)]
pub struct IntersectionsArgs {
    /// Path to BED/BEDGraph file with regions of interest
    #[arg(short = 'b', long = "bed-file", required = true)]
    pub bed_file: PathBuf,

    /// Output file prefix
    #[arg(short = 'o', long = "output")]
    pub output_prefix: Option<String>,

    /// Enable debug logging
    #[arg(short = 'd', long = "debug")]
    pub debug: bool,
}

#[derive(Parser, Debug)]
pub struct AnalyzeArgs {
    /// Path to reference genome FASTA file
    #[arg(short = 'r', long = "reference-fasta", required = true)]
    pub reference_fasta: PathBuf,

    /// Convert masked bases to uppercase and include them
    #[arg(short = 'u', long = "unmask-genome")]
    pub unmask_genome: bool,

    /// Minimum base quality (default: 20)
    #[arg(short = 'q', long = "base-quality", default_value = "20")]
    pub base_quality: u8,

    /// Exclude reads with indels
    #[arg(short = 'f', long = "filter-indel")]
    pub filter_indel: bool,

    /// Read length range (min,max) (default: 0,350)
    #[arg(short = 'l', long = "filter-length", value_parser = parse_range, default_value = "0,350")]
    pub filter_length: (u32, u32),

    /// Softclip handling: 0=check, 1=skip, 2=force use (default: 0)
    #[arg(short = 's', long = "soft-clip-bypass", default_value = "0")]
    pub soft_clip_bypass: u8,

    /// Minimum mapping quality (default: 20)
    #[arg(short = 'm', long = "mapping-quality", default_value = "20")]
    pub mapping_quality: u8,

    /// Fragment length range (min,max) (default: 0,10000)
    #[arg(short = 'g', long = "fragment-length", value_parser = parse_range, default_value = "0,10000")]
    pub fragment_length: (u32, u32),

    /// Output file prefix for HTML report
    #[arg(short = 'o', long = "output-prefix")]
    pub output_prefix: Option<String>,

    /// Confidence threshold for base counting (default: 20)
    #[arg(short = 'c', long = "confidence", default_value = "20")]
    pub confidence: u32,

    /// Enable debug logging
    #[arg(short = 'd', long = "debug")]
    pub debug: bool,

    /// Oxford Nanopore mode (long reads)
    #[arg(short = 'O', long = "ont")]
    pub ont: bool,

    /// Normalize tables using Picard logic
    #[arg(short = 'p', long = "picard-logic")]
    pub picard_logic: bool,

    /// Mask C->T (read 1) and G->A (read 2) as methylation artifacts
    #[arg(long = "mask-methyl-c")]
    pub mask_methyl_c: bool,

    /// Mask methylation only at CpG sites
    #[arg(long = "mask-methyl-cpg")]
    pub mask_methyl_cpg: bool,
}

/// Parse a comma-separated range like "0,350" into (min, max)
fn parse_range(s: &str) -> Result<(u32, u32), String> {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 2 {
        return Err(format!(
            "Expected format 'min,max', got '{}' ({} parts)",
            s,
            parts.len()
        ));
    }

    let min = parts[0]
        .trim()
        .parse::<u32>()
        .map_err(|e| format!("Invalid min value '{}': {}", parts[0], e))?;
    let max = parts[1]
        .trim()
        .parse::<u32>()
        .map_err(|e| format!("Invalid max value '{}': {}", parts[1], e))?;

    Ok((min, max))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_range() {
        assert_eq!(parse_range("0,350").unwrap(), (0, 350));
        assert_eq!(parse_range("100,500").unwrap(), (100, 500));
        assert!(parse_range("invalid").is_err());
        assert!(parse_range("1,2,3").is_err());
    }
}
