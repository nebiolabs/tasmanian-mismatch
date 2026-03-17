# Rustmanian-Mismatch

A high-performance tool for the analysis of mismatches in high-throughput sequencing data from a reference genome.

## Features

- **Fast parallel processing** for multi-threaded BAM analysis
- **BED file filtering** - Filter or mask reads based on genomic regions
- **Methylation-aware analysis** for bisulfite/EM-seq data
- **Paired-end overlap detection** and inconsistency analysis
- **Soft-clip handling** with configurable thresholds
- **Genomic variant calling** with depth filtering
- **SAM flag filtering** compatible with samtools

## Installation

### Build from source

```bash
# Clone the repository
cd /path/to/tasmanian-mismatch

# Build release version
cargo build --release

# Binary will be at: ./target/release/rustmanian-mismatch
```

### Dependencies

- Rust 1.70+
- rust-htslib 0.47
- rayon 1.10
- bio 1.6
- clap 4.5

## Usage

### Basic Usage

```bash
rustmanian-mismatch <BAM_FILE> <REFERENCE_FASTA> [OPTIONS]
```

### Common Options

```bash
-t, --threads <N>              Number of threads (0 = all cores)
-q <QUAL>                      Minimum base quality (default: 20)
--min-map-quality <MAPQ>       Minimum mapping quality (default: 10)l
-r, --region-size <SIZE>       Region size for parallel processing (default: 10000000)
--softclip-threshold <FRAC>    Fraction of matching bases in softclips (default: 0.66)
```

### BED File Filtering

Filter or mask reads based on genomic regions defined in a BED file:

```bash
# Mask individual bases in BED regions (default)
rustmanian-mismatch input.bam reference.fa \
  -b regions_to_mask.bed \
  --bed-filter-mode mask

# Filter entire reads overlapping BED regions
rustmanian-mismatch input.bam reference.fa \
  -b regions_to_exclude.bed \
  --bed-filter-mode filter
```

**Filter Modes:**
- `mask`: Skip individual bases within reads that overlap BED regions (keeps read)
- `filter`: Skip entire reads that overlap any BED region (removes read)

See [BED_FILTERING.md](BED_FILTERING.md) for detailed documentation.

### Methylation Analysis

```bash
# Basic methylation mode (all C>T in read 1)
rustmanian-mismatch input.bam reference.fa -m

# CpG-only methylation mode
rustmanian-mismatch input.bam reference.fa --cpg-only
```

### Flag Filtering (samtools compatible)

```bash
# Require specific flags
rustmanian-mismatch input.bam reference.fa -f 3  # Require paired and properly paired

# Filter out flags
rustmanian-mismatch input.bam reference.fa -F 1804  # Exclude unmapped, secondary, supplementary, duplicates

# Exclude if ALL flags are set
rustmanian-mismatch input.bam reference.fa -G 256  # Exclude secondary alignments
```

### Complete Example

```bash
rustmanian-mismatch \
  sample.bam \
  reference.fa \
  -t 16 \
  -q 30 \
  --min-map-quality 20 \
  -b problematic_regions.bed \
  --bed-filter-mode mask \
  -m \
  --genomic-threshold 7 \
  --genomic-depth-threshold 10 \
  > sample_mismatches.csv
```

## Output Files

### Standard Output (CSV)
Main mismatch counts table with columns:
- `Read`: Read number (1 or 2)
- `Position`: Position in read
- `REF>ALT`: Count of each mismatch type

### potential_variants.tsv
Genomic positions with high mismatch counts:
- `chromosome`: Chromosome name
- `position`: Genomic position
- `reference_base`: Reference base
- `mismatch_base`: Observed base
- `count`: Number of mismatches
- `depth`: Total depth at position

### Overlap Analysis (if paired reads)
Additional sections in output:
- Overlap region mismatches
- Read pair inconsistencies

## Performance

- **Speed**: ~10-50x faster than Python version depending on dataset
- **Memory**: Efficient memory usage with streaming BAM processing
- **Scalability**: Linear scaling with number of threads
- **BED filtering**: Minimal overhead (~5-10% with mask mode)

## Example Workflow

```bash
# 1. Index your BAM file
samtools index input.bam

# 2. Create BED file with regions to exclude (optional)
cat > exclude.bed <<EOF
chr1    1000000    1001000
chr2    5000000    5002000
EOF

# 3. Run analysis
rustmanian-mismatch input.bam reference.fa \
  -b exclude.bed \
  -t 8 \
  -q 30 \
  > results.csv 2> analysis.log

# 4. Check potential variants
cat potential_variants.tsv
```

## Comparison with Python Version

| Feature | Rust Version | Python Version |
|---------|-------------|----------------|
| Speed | ✓✓✓ | ✓ |
| Memory | ✓✓✓ | ✓✓ |
| Parallelization | Native | Limited |
| BED filtering | ✓ (built-in) | ✓ (external) |
| Methylation | ✓ | ✓ |
| Dependencies | Minimal | Many |

## License

AGPL v3

## Citation

If you use this tool, please cite:
[Original Tasmanian paper/repository]

## Related Tools

- [Python version](python_version/) - Original implementation
- [samtools](http://www.htslib.org/) - SAM/BAM manipulation
- [bedtools](https://bedtools.readthedocs.io/) - Genomic interval operations


## ToDo

  - Add option to report by read or by insert/fragment (main.rs:109), where the insert (in between the reads) can be 1 dash or many, depending on  the length of the reads, so that the end of the    reads can be stack and he result is the sum of observations at specific positions FROM the ends.  There is already an unused variable **_insert_position_mode.**
  - We can already filter through the flag (which can identify read number and strand to which it aligns to). It would be great if we can also filter through RF/FR (read that aligns to + has lower grenomic coordinate = FR (fwd/rev), and opposite for RF). 
  - In processing.rs:140, __adjust_methylation_base__ is applied if tracking genomic positions (for cheap variant calling). However, that function has been already applied within __create_mismatch_key__ (line 112). Perhaps we can output that info and avoid running it twice (for each C>T found!)
  - In each chunk, there will be a few reads without pair, because the bam is probably sorted by genomic position rather than qname. 2 options: 1. Brad proposed exceeding the loop boundary - e.g. if we loop over 10k, keep going beyond, until we find the mate or we hit a limit...  say 11K, as any fragment should not be that long. 2. Keep these reads in the heap until the end of the program and loop through these few reads.  
