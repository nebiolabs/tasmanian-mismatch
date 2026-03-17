# BED File Filtering Guide

## Overview

This tool now supports filtering and masking reads based on BED file regions, leveraging both custom implementation and htslib-compatible approaches.

## Usage

### Basic Command with BED Filtering

```bash
# Mask individual bases that overlap with BED regions (default)
./target/release/rustmanian-mismatch \
  input.bam \
  reference.fa \
  -b regions.bed \
  --bed-filter-mode mask

# Filter entire reads that overlap with BED regions
./target/release/rustmanian-mismatch \
  input.bam \
  reference.fa \
  -b regions.bed \
  --bed-filter-mode filter
```

## Filter Modes

### 1. Mask Mode (default: `--bed-filter-mode mask`)
- Skips individual bases within reads that overlap BED regions
- Keeps the read but excludes specific positions from analysis
- Useful for masking known problematic regions while retaining read-level information
- Best for: Excluding repetitive regions, low-complexity regions, or known artifacts

### 2. Filter Mode (`--bed-filter-mode filter`)
- Skips entire reads if they overlap with any BED region
- More aggressive filtering approach
- Useful for completely excluding reads from specific genomic regions
- Best for: Removing reads from entire genes, chromosomes, or large problematic regions

## BED File Format

The tool accepts standard BED format files:

```
chr1    1000    2000    region1
chr1    5000    6000    region2
chr2    100     500     region3
```

Minimum required columns:
1. Chromosome name
2. Start position (0-based)
3. End position (0-based, exclusive)

Additional columns are ignored but preserved compatibility with standard BED files.

## Examples

### Example 1: Mask repetitive regions

```bash
# Create BED file with repetitive regions to mask
cat > repetitive_regions.bed <<EOF
chr1    1000000    1001000
chr1    5000000    5002000
chr2    3000000    3000500
EOF

# Run analysis with masking
rustmanian-mismatch input.bam reference.fa -b repetitive_regions.bed
```

### Example 2: Filter reads overlapping with specific genes

```bash
# Use UCSC or Ensembl BED file for genes to exclude
rustmanian-mismatch input.bam reference.fa \
  -b genes_to_exclude.bed \
  --bed-filter-mode filter
```

### Example 3: Combined with other filters

```bash
# Combine BED filtering with quality filters
rustmanian-mismatch input.bam reference.fa \
  -b problematic_regions.bed \
  -q 30 \
  --min-map-quality 20 \
  -t 8
```

## Using with samtools-style filtering

For compatibility with samtools-style workflows, you can pre-filter your BAM file:

```bash
# htslib approach: exclude regions (creates a new BAM)
samtools view -L ^exclude_regions.bed -b input.bam > filtered.bam

# Then run rustmanian-mismatch on filtered BAM
rustmanian-mismatch filtered.bam reference.fa
```

Or use the built-in BED filtering:

```bash
# Direct approach: use built-in BED filtering
rustmanian-mismatch input.bam reference.fa -b exclude_regions.bed
```

## Performance Notes

- BED regions are loaded once at startup and indexed by chromosome
- Intervals are sorted for efficient binary search
- Minimal overhead for mask mode (~5-10%)
- Filter mode is slightly faster as it skips entire reads

## Integration with htslib

This implementation is compatible with htslib's BED handling conventions:
- 0-based coordinates (same as BAM)
- Half-open intervals [start, end)
- Standard BED format support
- Compatible with bedtools and samtools workflows

## Common Use Cases

1. **Exclude centromeric regions**: Filter out high-error regions
2. **Mask low-complexity regions**: Skip repetitive DNA
3. **Remove blacklisted regions**: Exclude ENCODE blacklist regions
4. **Gene-specific analysis**: Filter reads from specific genes
5. **Sex chromosome filtering**: Exclude X/Y chromosomes if needed

## Notes

- Empty BED files or files with no matching chromosomes are handled gracefully
- Chromosome name mismatches (chr1 vs 1) are NOT automatically handled - ensure consistency
- Overlapping BED intervals are supported (any overlap triggers filtering)
- BED file is validated at load time with helpful error messages
