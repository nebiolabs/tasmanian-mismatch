# Rustmanian-Mismatch

[![CI](https://github.com/nebiolabs/tasmanian-mismatch/actions/workflows/ci.yml/badge.svg)](https://github.com/nebiolabs/tasmanian-mismatch/actions/workflows/ci.yml)

## Abstract

DNA preservation, extraction, fragmentation, and library preparation can introduce systematic sequencing errors that are distinct from true biological variation. These artifacts often follow characteristic patterns tied to sequence context, strand orientation, or position within a read — for example, FFPE samples are known to produce read-position-dependent mismatch signatures. Tasmanian-Mismatch was built to characterize and quantify such error patterns across sequencing datasets, helping researchers distinguish technical artifacts from genuine genomic variation. The resulting models can also be used to recalibrate base-quality scores, reducing the technical noise passed to downstream variant callers and improving the accuracy of variant discovery.

Tasmanian-Mismatch is a toolkit for mismatch analysis on indexed BAM files against a reference FASTA. The repository currently builds three binaries:

- `tasmanian-mismatch`: count mismatches by read position or fragment position
- `tasmanian-diagnostics`: report genomic mismatch sites, overlap inconsistencies, and read-position discount tables
- `tasmanian-rescale-quality`: rescale BAM quality scores from a tab-delimited matrix

## Example Visualization
![Visualization](test_data_visualization.png)
[Interactive version](test_data_visualization.html)
Create this using
```bash
pixi run python scripts/visualize.py tests/fixtures/test_data.tsv
```

## Build

```bash
cargo build --release
```

Release binaries:

```text
./target/release/tasmanian-mismatch
./target/release/tasmanian-diagnostics
./target/release/tasmanian-rescale-quality
```

## Inputs

- Coordinate-sorted BAM with index (`.bai`)
- Reference FASTA
- Optional BED file for masking or whole-read filtering
- Optional discount table from `tasmanian-diagnostics`
- Optional rescaling matrix for `tasmanian-rescale-quality`

## Binary Overview

### `tasmanian-mismatch`

Count mismatch classes such as `C>T`, `G>A`, `A>C` by either read position or fragment position.

Basic usage:

```bash
tasmanian-mismatch <BAM> <REFERENCE_FASTA> [OPTIONS]
```

Common options:

```bash
-t, --threads <N>                 Number of threads (0 keeps rayon default)
-r, --region-size <BP>            Region size for indexed chunking
-q, --min-base-quality <QUAL>     Minimum base quality
-m, --min-map-quality <MAPQ>      Minimum mapping quality
--position-mode <read|insert>     Position mode (default: insert)
--overlap-mode <cut|stretch>      Overlap handling mode
--discount-table <TSV>            Discount table from tasmanian-diagnostics
-b, --bed-file <BED>              BED file for masking/filtering
--bed-filter-mode <mask|filter>   BED handling mode
-f <FLAGS>                        SAM flags that must be present
-F <FLAGS>                        SAM flags that, if present, skip a read
-G <FLAGS>                        SAM flags that, if all present, skip a read
--min-fragment-length <LEN>       Minimum fragment length for insert mode
--max-fragment-length <LEN>       Maximum fragment length for insert mode
--methylation-mode                Collapse methylation-driven mismatch classes
--normalize                       Write normalized frequencies instead of raw counts
--emit-rescaling-matrix           Emit matrix rows for tasmanian-rescale-quality
--plot                            Launch the optional Bokeh visualization
-o, --output-file <TSV>           Output path
```

Examples:

```bash
# Raw mismatch counts by fragment position
tasmanian-mismatch sample.bam reference.fa \
  --position-mode insert \
  -o mismatch_counts.tsv

# Read-position counts with a diagnostics-derived discount table
tasmanian-mismatch sample.bam reference.fa \
  --position-mode read \
  --discount-table variant_discounts.tsv \
  -o mismatch_counts.tsv

# Pipe discount table via stdin (use '-' to read --discount-table from stdin)
cat variant_discounts.tsv | tasmanian-mismatch sample.bam reference.fa \
  --position-mode read \
  --discount-table - \
  -o mismatch_counts.tsv

# Normalized frequencies instead of integer counts
tasmanian-mismatch sample.bam reference.fa \
  --position-mode read \
  --normalize \
  -o mismatch_normalized.tsv

# Emit rescaling matrix rows to stdout (or -o file.tsv)
tasmanian-mismatch sample.bam reference.fa \
  --position-mode read \
  --emit-rescaling-matrix
```

### `tasmanian-diagnostics`

Produce three diagnostic outputs:

- genomic mismatch sites (`potential_variants.tsv`)
- paired-read overlap inconsistencies (`read_pair_inconsistencies.tsv`)
- mismatch discount table (`variant_discounts.tsv`)

Basic usage:

```bash
tasmanian-diagnostics <BAM> <REFERENCE_FASTA> [OPTIONS]
```

Common options:

```bash
-t, --threads <N>                    Number of threads (0 keeps rayon default)
-r, --region-size <BP>                Region size for parallel processing
-q, --min-base-quality <QUAL>         Minimum base quality
--min-map-quality <MAPQ>              Minimum mapping quality
-m, --methylation                     Convert C/T in read 1 back to C (bisulfite/EM-seq)
--use-insert-mode                     Use fragment-level instead of read-level positions
--min-fragment-length <LEN>           Minimum estimated fragment length for a read to be counted
--max-fragment-length <LEN>           Maximum estimated fragment length for a read to be counted
--genomic-threshold <N>               Minimum mismatch count for reporting a genomic site
--genomic-depth-threshold <N>         Minimum depth for reporting a genomic site
-b, --bed-file <BED>                  BED file for masking/filtering
--bed-filter-mode <mask|filter>       BED handling mode
-f <FLAGS>                            SAM flags that must be present
-F <FLAGS>                            SAM flags that, if present, skip a read
-G <FLAGS>                            SAM flags that, if all present, skip a read
--variants-output <TSV>               Genomic variants output (default: potential_variants.tsv)
--inconsistencies-output <TSV>        Overlap inconsistencies output (default: read_pair_inconsistencies.tsv)
--discount-output <TSV>               Discount table output (default: variant_discounts.tsv), use '-' for stdout
```

Example:

```bash
tasmanian-diagnostics sample.bam reference.fa \
  --variants-output potential_variants.tsv \
  --inconsistencies-output read_pair_inconsistencies.tsv \
  --discount-output variant_discounts.tsv
```

Direct piping to `tasmanian-mismatch` is supported by writing discounts to stdout:

```bash
tasmanian-diagnostics sample.bam reference.fa \
  --variants-output potential_variants.tsv \
  --inconsistencies-output read_pair_inconsistencies.tsv \
  --discount-output - | \
tasmanian-mismatch sample.bam reference.fa \
  --position-mode read \
  --discount-table - \
  -o mismatch_counts.tsv
```

The `variant_discounts.tsv` output can be consumed by `tasmanian-mismatch` using `--discount-table`.

### `tasmanian-rescale-quality`

Rewrite BAM base qualities from a rescaling matrix keyed by read number, read position, reference base, and observed read base.

Basic usage:

```bash
tasmanian-rescale-quality <BAM> <REFERENCE_FASTA> <MATRIX_TSV> [OPTIONS]
```

Use `-` as `<MATRIX_TSV>` to read matrix rows from stdin.

Example:

```bash
tasmanian-rescale-quality sample.bam reference.fa quality_matrix.tsv \
  -o rescaled.bam

# Consume matrix from stdin
cat quality_matrix.tsv | tasmanian-rescale-quality sample.bam reference.fa - \
  -o rescaled.bam
```

## Workflow

The intended workflow is:

1. Run `tasmanian-mismatch` to get raw mismatch counts.
2. Optionally run `tasmanian-mismatch --normalize` to get within-group mismatch frequencies.
3. Optionally convert normalized frequencies into a rescaling matrix for `tasmanian-rescale-quality`.
4. Run `tasmanian-diagnostics` when you need genomic-site summaries or a discount table.
5. Feed `variant_discounts.tsv` back into `tasmanian-mismatch` via `--discount-table` if desired.

Important distinction:

- `tasmanian-mismatch` output is not the direct input to `tasmanian-rescale-quality`
- `tasmanian-rescale-quality` expects a matrix with columns `read_num`, `position`, `ref_base`, `read_base`, `scaling_factor`
- `variant_discounts.tsv` is for `tasmanian-mismatch`, not for `tasmanian-rescale-quality`

Direct pipe from mismatch into rescale-quality:

```bash
tasmanian-mismatch sample.bam reference.fa \
  --position-mode read \
  --emit-rescaling-matrix | \
tasmanian-rescale-quality sample.bam reference.fa - \
  -o rescaled.bam
```

## Output Formats

### `tasmanian-mismatch` raw output

Read-position mode:

```tsv
base_change	read_num	reference_order	read_position	count
C>T	1	1	42	18
C>A	1	1	42	3
G>A	2	2	17	11
```

Fragment-position mode uses the same columns except `read_position` becomes `fragment_position`.

### `tasmanian-mismatch --normalize`

```tsv
base_change	read_num	reference_order	read_position	normalized_frequency
C>T	1	1	42	0.782609
C>A	1	1	42	0.130435
C>G	1	1	42	0.043478
```

Normalization is performed within each `(read_num, position, ref_base)` group. For example:

```text
C>T / (C>A + C>C + C>G + C>T)
```

If a group's total count is zero (e.g. after `--discount-table` has reduced every count in the group to zero), its frequency is reported as `0.0` rather than `NaN`.

### `potential_variants.tsv`

```tsv
chromosome	position	reference_base	mismatch_base	count	depth
chr1	10452	C	T	12	38
chr1	20891	G	A	9	27
chr2	450103	A	C	15	41
```

### `read_pair_inconsistencies.tsv`

```tsv
read1_position	read2_position	discordance_type	count
18	83	R1:A_R2:G	5
19	82	R1:C_R2:T	3
20	81	R1:G_R2:A	7
```

### `variant_discounts.tsv`

```tsv
mismatch_type	read_num	read_position	discount_count
C>T	1	42	6
G>A	2	17	4
A>C	1	88	9
```

This table is read-position keyed and can be supplied to `tasmanian-mismatch` through `--discount-table`.

### Rescaling matrix for `tasmanian-rescale-quality`

```tsv
1	42	C	T	0.85
1	42	C	A	1.10
2	17	G	A	0.65
```

Columns:

```text
read_num    position    ref_base    read_base    scaling_factor
```

## BED Filtering

BED-aware handling is available in both `tasmanian-mismatch` and `tasmanian-diagnostics`.

```bash
# Mask bases overlapping BED intervals
tasmanian-mismatch sample.bam reference.fa \
  -b regions.bed \
  --bed-filter-mode mask

# Drop whole reads overlapping BED intervals
tasmanian-mismatch sample.bam reference.fa \
  -b regions.bed \
  --bed-filter-mode filter
```

Modes:

- `mask`: skip only aligned positions overlapping BED intervals
- `filter`: skip the whole read if any aligned portion overlaps a BED interval

See [BED_FILTERING.md](BED_FILTERING.md) for details.

## Integration Tests

Current integration coverage:

- [tests/mismatch_integration.rs](tests/mismatch_integration.rs): fixture-driven test for `tasmanian-mismatch`
- [tests/diagnostics_integration.rs](tests/diagnostics_integration.rs): fixture-driven test for `tasmanian-diagnostics`
- [tests/rescale_integration.rs](tests/rescale_integration.rs): end-to-end test for `tasmanian-rescale-quality`
- [tests/all_binaries_and_pipes_integration.rs](tests/all_binaries_and_pipes_integration.rs): end-to-end test piping all three binaries together

Each integration test synthesizes its own minimal BAM fixture at runtime (into a temp directory); no BAM files are checked into the repo. Reusable helpers live in [tests/test_utils.rs](tests/test_utils.rs).

Run the full test suite with:

```bash
cargo test --all-targets
```

## Notes

- The `-m` short flag differs by binary: it's `--min-map-quality` on `tasmanian-mismatch` but `--methylation` on `tasmanian-diagnostics`. Double-check which binary you're invoking when using short flags.
- `tasmanian-mismatch` and `tasmanian-diagnostics` currently use slightly different emitted position conventions in some outputs; tests document the current behavior.
- `frequencies_to_rescaling_matrix` exists as a placeholder conversion step in the library and currently emits `1.0` scaling factors.
- `reference_order` in mismatch output indicates which read in a pair appears first in reference coordinates.

## License
AGPL v3
