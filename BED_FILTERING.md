# BED File Filtering Guide

## Overview

This tool supports filtering and masking reads based on BED file regions, leveraging both custom implementation and htslib-compatible approaches.

## Usage

### Basic Command with BED Filtering

```bash
# Mask individual bases that overlap with BED regions (default)
./target/release/tasmanian-mismatch \
  input.bam \
  reference.fa \
  -b regions.bed \
  --bed-filter-mode mask

# Filter entire reads that overlap with BED regions
./target/release/tasmanian-mismatch \
  input.bam \
  reference.fa \
  -b regions.bed \
  --bed-filter-mode filter

# Keep ONLY reads that overlap with BED regions (in-silico exome/panel restriction)
./target/release/tasmanian-mismatch \
  input.bam \
  reference.fa \
  -b exome_targets.bed \
  --bed-filter-mode include
```

BED filtering is available in both `tasmanian-mismatch` and `tasmanian-diagnostics`, with the
same three modes. When `tasmanian-mismatch` runs with `--discount-table`, run
`tasmanian-diagnostics` with the same BED and mode so the discounts come from the same reads.

`--bed-filter-mode` requires `--bed-file`; giving a mode without a BED is rejected at argument
parsing, as is an unknown mode name.

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

### 3. Include Mode (`--bed-filter-mode include`)
- The exact inverse of Filter mode: skips a read UNLESS it overlaps a BED region
- Best for: in-silico exome or targeted-panel restriction -- keeping only reads over a
  specific set of regions (e.g. protein-coding CDS) rather than excluding regions
- `mask`/`filter` can only exclude BED regions; there is no way to express "include only"
  by inverting them (masking or filtering the *complement* of a BED file doesn't give the
  same result for reads straddling a region boundary), so this is a genuinely separate mode
- Like Filter mode, this check is whole-read (not per-base like Mask mode)

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

### Parsing tolerance

- Header/comment lines (or any line with a non-numeric start/end column) are skipped with a warning instead of aborting the run.
- Intervals are not required to be pre-sorted: if a chromosome's intervals arrive unsorted by start position, they are sorted automatically (with a warning).
- Overlapping, nested, or adjacent intervals on the same chromosome are merged into non-overlapping intervals at load time. This is what makes the binary-search overlap check both correct and fast.

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
tasmanian-mismatch input.bam reference.fa \
  -b genes_to_exclude.bed \
  --bed-filter-mode filter
```

### Example 3: Combined with other filters

```bash
# Combine BED filtering with quality and fragment-length filters
tasmanian-mismatch input.bam reference.fa \
  -b problematic_regions.bed \
  -q 30 \
  --min-map-quality 20 \
  --min-fragment-length 25 \
  --max-fragment-length 10000 \
  -t 8
```

Fragment-length filtering (`--min-fragment-length` / `--max-fragment-length`, defaults 25/10000) excludes reads whose estimated fragment length falls outside the given range. It is independent of BED filtering but commonly used alongside it to further restrict which reads contribute counts.

### Example 4: In-silico exome or panel restriction

```bash
# Keep only reads over a set of protein-coding (CDS) or capture-panel regions
tasmanian-mismatch input.bam reference.fa \
  -b exome_targets.bed \
  --bed-filter-mode include
```

## Using with samtools-style filtering

For compatibility with samtools-style workflows, you can pre-filter your BAM file instead of
using `--bed-filter-mode`:

```bash
# htslib approach: exclude regions (creates a new BAM)
samtools view -L ^exclude_regions.bed -b input.bam > filtered.bam

# htslib approach: keep only regions (samtools -L is inclusion-only by default)
samtools view -L exome_targets.bed -b input.bam > exome_only.bam

# Then run tasmanian-mismatch on the pre-filtered BAM, with no -b needed
tasmanian-mismatch filtered.bam reference.fa
```

Or use the built-in BED filtering, which avoids writing an intermediate BAM:

```bash
# Direct approach: use built-in BED filtering
tasmanian-mismatch input.bam reference.fa -b exclude_regions.bed --bed-filter-mode filter
tasmanian-mismatch input.bam reference.fa -b exome_targets.bed --bed-filter-mode include
```

Both the BED intervals and read alignment spans are half-open `[start, end)`, the same
convention `samtools view -L` uses: a read overlaps an interval only if it shares at least one
reference base with it, and a read that merely abuts one does not. The built-in filter therefore
selects the same reads as a `samtools view -L` pre-filter. `include` and `filter` are exact
complements: on the same BAM and BED, every read that passes the other filters is kept by
exactly one of them.

## Performance Notes

- BED regions are loaded once at startup and indexed by chromosome
- Intervals are sorted for efficient binary search
- Minimal overhead for mask mode (~5-10%)
- Filter mode is slightly faster as it skips entire reads
- Include mode shares Filter mode's per-read overlap check (just inverted), so its
  performance characteristics are the same

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
6. **In-silico exome/panel restriction**: Keep only reads over a set of CDS or capture-panel
   regions (`--bed-filter-mode include`)

## Notes

- Empty BED files or files with no matching chromosomes are handled gracefully
- Chromosome name mismatches (chr1 vs 1) are NOT automatically handled - ensure consistency
- Overlapping, nested, and adjacent BED intervals are merged automatically at load time (any overlap with the merged interval triggers filtering)
- Header lines and unsorted input are tolerated (see "Parsing tolerance" above)
- BED file is validated at load time with helpful error messages
