# Tasmanian Rust Rewrite Plan

## Overview

Rewrite the Tasmanian mismatch analysis tool in Rust using the noodles library for SAM/BAM I/O. The tool will maintain Unix philosophy: reading from stdin and outputting to stdout.

## Architecture

### Project Structure

```
tasmanian-rs/
├── Cargo.toml
├── src/
│   ├── main.rs              # Entry point, CLI handling
│   ├── lib.rs               # Library root
│   ├── cli.rs               # CLI argument definitions (clap)
│   ├── bed.rs               # BED file parsing
│   ├── reference.rs         # FASTA reference loading
│   ├── sam/
│   │   ├── mod.rs           # SAM module root
│   │   ├── record.rs        # Read record wrapper
│   │   └── flags.rs         # SAM flag handling
│   ├── intersection.rs      # BED-read intersection logic
│   ├── analysis.rs          # Artifact analysis engine
│   ├── error_table.rs       # Mismatch counting tables
│   └── output.rs            # CSV output formatting
└── tests/
    └── integration_test.rs
```

### Dependencies

```toml
[dependencies]
noodles = { version = "0.85", features = ["sam", "fasta"] }
clap = { version = "4", features = ["derive"] }
anyhow = "1.0"
thiserror = "2.0"
csv = "1.3"
log = "0.4"
env_logger = "0.11"
```

## Core Components

### 1. CLI Interface (`cli.rs`)

Two subcommands matching the Python tools:

```
tasmanian intersections -b <bed_file>    # Preprocessor
tasmanian analyze -r <reference.fa>      # Main analyzer
```

#### `intersections` Options:
- `-b, --bed-file` (required): Path to BED/BEDGraph file
- `-o, --output-prefix`: Output file prefix
- `-d, --debug`: Enable debug logging

#### `analyze` Options:
- `-r, --reference-fasta` (required): Path to reference genome
- `-u, --unmask-genome`: Include masked (lowercase) bases
- `-q, --base-quality <N>`: Minimum base quality (default: 20)
- `-f, --filter-indel`: Exclude reads with indels
- `-l, --filter-length <MIN,MAX>`: Read length range (default: 0,350)
- `-s, --soft-clip-bypass <0|1|2>`: Softclip handling mode
- `-m, --mapping-quality <N>`: Minimum mapping quality (default: 20)
- `-g, --fragment-length <MIN,MAX>`: Fragment length filter
- `-o, --output-prefix`: Output file prefix
- `-c, --confidence <N>`: Confidence threshold for base counting
- `-d, --debug`: Enable debug logging
- `-O, --ont`: Oxford Nanopore mode (long reads)
- `--mask-methyl-c`: Mask C->T (R1) and G->A (R2) artifacts
- `--mask-methyl-cpg`: Mask methylation only at CpG sites

### 2. BED File Parsing (`bed.rs`)

```rust
pub struct BedRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Option<char>,
    pub name: Option<String>,
    pub class: Option<String>,
    pub family: Option<String>,
}

pub struct BedIndex {
    regions: HashMap<String, Vec<BedRegion>>,
    // Sorted by start position per chromosome
}
```

### 3. Reference Loading (`reference.rs`)

```rust
pub struct Reference {
    sequences: HashMap<String, Vec<u8>>,
}

impl Reference {
    pub fn load(path: &Path) -> Result<Self>;
    pub fn get(&self, chrom: &str, start: u64, end: u64) -> Option<&[u8]>;
}
```

### 4. SAM Record Processing (`sam/record.rs`)

```rust
pub struct TasmanianRecord {
    pub id: String,
    pub flag: u16,
    pub chrom: String,
    pub start: u64,
    pub mapq: u8,
    pub cigar: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub tlen: i64,

    // Intersection metadata
    pub category: Option<u8>,
    pub intersection_bounds: Option<(usize, usize)>,
    pub complement_bases: usize,
    pub junction: Option<String>,
}

impl TasmanianRecord {
    pub fn from_noodles(record: &Record, header: &Header) -> Result<Self>;
    pub fn is_proper_pair(&self) -> bool;
    pub fn read_number(&self) -> Option<u8>;  // 1 or 2
    pub fn is_reverse(&self) -> bool;
}
```

### 5. Intersection Detection (`intersection.rs`)

Categories:
- **1**: Unrelated (no intersection)
- **2a-d**: Only read 1 intersects (a=complete, b=left, c=right, d=both sides)
- **3a-d**: Only read 2 intersects
- **4xy**: Both reads intersect same region
- **5xy**: Both reads intersect different regions

```rust
pub fn detect_intersections(
    record: &mut TasmanianRecord,
    bed_index: &BedIndex,
    current_bed_idx: &mut usize,
) -> Option<IntersectionResult>;
```

### 6. Error Tables (`error_table.rs`)

```rust
pub struct ErrorTable {
    // [read_number][position][ref_base][alt_base] = count
    counts: Vec<Vec<[[u64; 4]; 4]>>,  // [2][read_length][4][4]
}

pub struct ErrorTables {
    intersection: ErrorTable,      // Bases overlapping BED regions
    complement: ErrorTable,        // Non-overlapping bases from intersecting reads
    unrelated: ErrorTable,         // All bases from non-intersecting reads
    intersection_conf: ErrorTable, // High-confidence intersection
    complement_conf: ErrorTable,   // High-confidence complement
}
```

### 7. Analysis Engine (`analysis.rs`)

```rust
pub struct Analyzer {
    reference: Reference,
    config: AnalysisConfig,
    tables: ErrorTables,
    read_lengths: Vec<usize>,  // For mode calculation
}

impl Analyzer {
    pub fn process_record(&mut self, record: &TasmanianRecord) -> Result<()>;
    pub fn finalize(&mut self) -> Result<ErrorTables>;
}
```

Key processing steps:
1. Parse CIGAR to handle indels/softclips
2. Align sequence to reference
3. For each base: quality check, reverse complement if needed
4. Categorize by intersection status
5. Increment appropriate error table

### 8. CSV Output (`output.rs`)

Output format matches Python:
```
read,position,Ia_a,Ia_t,Ia_c,Ia_g,...,Ca_a,...,Na_a,...,cIa_a,...,cCa_a,...
```

80 data columns (5 tables × 16 mutation types) plus read and position.

## Data Flow

### Pipeline Mode (Full)
```
samtools view data.bam | tasmanian intersections -b regions.bed | tasmanian analyze -r ref.fa > output.csv
```

### Simple Mode (No BED)
```
samtools view data.bam | tasmanian analyze -r ref.fa > output.csv
```

## Implementation Notes

### Reading from stdin with noodles

```rust
use noodles::sam;
use std::io::{self, BufReader};

let stdin = io::stdin().lock();
let mut reader = sam::io::Reader::new(BufReader::new(stdin));
let header = reader.read_header()?;

for result in reader.records() {
    let record = result?;
    // Process record
}
```

### SAM Tags for Intersection Data

Custom tags added by `intersections`:
- `tm:Z:a.b;start.end` - Junction coordinates
- `tc:i:N` - Number of confident bases in complement

### CIGAR Parsing

noodles provides CIGAR parsing, but we need custom logic for:
- Softclip detection and handling
- Tracking reference vs sequence position offsets
- Detecting "garbage" softclips

### Reverse Complement Handling

When flag indicates reverse strand:
- Transform reference base: A↔T, C↔G
- Transform called base: A↔T, C↔G
- Read position counts from 3' end

### Methylation Masking

When enabled:
- Read 1: C→T artifacts masked (becomes C)
- Read 2: G→A artifacts masked (becomes G)
- CpG mode: only at CpG dinucleotides

## Performance Considerations

1. **Reference Loading**: Load entire reference into memory (matches Python behavior)
2. **BED Index**: Pre-sort by chromosome and position for efficient lookup
3. **Streaming**: Process records one at a time, don't buffer entire file
4. **Parallel Processing**: Future enhancement - process chromosomes in parallel

## Testing Strategy

1. Use existing test data in `tasmanian/tests/`
2. Compare output with Python implementation
3. Test edge cases:
   - Empty BED files
   - No intersections
   - All intersections
   - Various CIGAR patterns (indels, softclips)
   - Both strands

## Migration Path

1. Phase 1: Core functionality (MVP)
   - Basic SAM reading from stdin
   - Reference loading
   - Simple mismatch counting (no BED)
   - CSV output

2. Phase 2: BED intersection support
   - BED file parsing
   - Intersection detection
   - Category assignment
   - Multiple error tables

3. Phase 3: Advanced features
   - Methylation masking
   - Oxford Nanopore mode
   - HTML report generation (optional)
   - PWM analysis (optional)

## Success Criteria

- [x] Processes SAM from stdin
- [x] Outputs CSV to stdout
- [x] Results match Python implementation on test data
- [x] Handles all CLI options
- [x] Proper error handling and logging
- [ ] Documentation and examples (in progress)

## Implementation Notes

### Completed (MVP - December 2024)

1. **Core Analysis Engine** (`analyze` subcommand)
   - Reads SAM from stdin using noodles library
   - Loads FASTA reference genome
   - Processes CIGAR alignments including indels and softclips
   - Counts mismatches in 5 error tables (intersection, complement, unrelated, confident intersection, confident complement)
   - Outputs CSV to stdout in same format as Python implementation
   - Supports all CLI options (quality filters, fragment length, methylation masking, ONT mode)

2. **Intersection Processing** (`intersections` subcommand)
   - Loads BED files with optional metadata
   - Detects read-region intersections
   - Adds `tm:Z:` and `tc:i:` tags to SAM output
   - Applies lowercase masking to intersecting bases
   - Uses efficient line-based processing (avoids noodles RecordBuf complexity)

3. **SAM/noodles API Usage**
   - Uses noodles 0.85 for SAM parsing
   - Handles proper paired reads (flags 99, 163, 83, 147)
   - Parses custom Tasmanian tags from pre-processed data
   - Note: SAM input requires header lines (@HD, @SQ) for proper parsing

### Future Enhancements

- HTML report generation (Plotly integration)
- PWM analysis output
- Performance optimization for large files
- Async I/O support
