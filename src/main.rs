use rust_htslib::bam::{Read, Reader, IndexedReader, Record, FetchDefinition};
use rayon::prelude::*;
use std::env;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::collections::HashMap;
use bio::io::fasta;
use std::sync::Arc;


// We are still not filtering out based on TLEN (python version does this)
// We are not skipping the whole read when INDELS (python version does this)
// Still needs to check if it's handling ONT data properly.


#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct MismatchKey {
    mismatch_type: String,  // e.g., "A>G", "C>T"
    read_position: usize,   // Position in the read (0-based)
    read_num: u8,           // 1 or 2 for paired-end reads
}

// Track inconsistencies between read1 and read2 in overlap regions
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct InconsistencyKey {
    discordance_type: String,  // e.g., "R1:A_R2:G"
    read1_position: usize,
    read2_position: usize,
}

// Store read information for overlap detection
#[derive(Debug, Clone)]
struct ReadInfo {
    qname: Vec<u8>,
    tid: i32,
    pos: i64,
    seq_len: usize,
    cigar: Vec<u8>,  // Store cigar as bytes for later use
    is_reverse: bool,
    read_num: u8,
}

// Store reference sequences in memory
type ReferenceGenome = HashMap<String, Vec<u8>>;

fn load_reference_genome(fasta_path: &str) -> ReferenceGenome {
    eprintln!("Loading reference genome from: {}", fasta_path);
    let reader = fasta::Reader::from_file(fasta_path)
        .expect("Failed to open reference FASTA file");
    
    let mut genome: ReferenceGenome = HashMap::new();
    for result in reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let chr_name = record.id().to_string();
        let sequence = record.seq().to_vec(); //Vec<u8> = byte, not UTF-8 char (overhead) 
        genome.insert(chr_name, sequence);
    }
    
    eprintln!("Loaded {} chromosome(s)", genome.len());
    genome
}

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => base,
    }
}

fn base_to_char(byte: u8) -> Option<char> {
    match byte {
        b'A' => Some('A'),
        b'C' => Some('C'),
        b'G' => Some('G'),
        b'T' => Some('T'),
        _ => None,
    }
}

fn parse_md_tag(md_string: &str) -> (Vec<(usize, char)>, Vec<(usize, usize)>) {
    // Parse MD tag and return:
    // 1. (position_offset, ref_base) for mismatches
    // 2. (start_pos, length) for matches
    let mut mismatches = Vec::new();
    let mut matches = Vec::new();
    let mut pos = 0;
    let mut num_str = String::new();
    
    for ch in md_string.chars() {
        if ch.is_numeric() {
            num_str.push(ch);
        } else if ch == '^' {
            // Deletion - skip the deleted bases
            if !num_str.is_empty() {
                let match_len = num_str.parse::<usize>().unwrap_or(0);
                if match_len > 0 {
                    matches.push((pos, match_len));
                    pos += match_len;
                }
                num_str.clear();
            }
        } else if ch.is_alphabetic() {
            // This is a reference base at a mismatch position
            if !num_str.is_empty() {
                let match_len = num_str.parse::<usize>().unwrap_or(0);
                if match_len > 0 {
                    matches.push((pos, match_len));
                    pos += match_len;
                }
                num_str.clear();
            }
            mismatches.push((pos, ch));
            pos += 1;
        }
    }
    
    // Handle trailing matches
    if !num_str.is_empty() {
        let match_len = num_str.parse::<usize>().unwrap_or(0);
        if match_len > 0 {
            matches.push((pos, match_len));
        }
    }
    
    (mismatches, matches)
}

fn compute_read_len_mode_from_sample_bam(bam_path: &str, sample_size: usize) -> usize {
    let mut bam = Reader::from_path(bam_path)
        .expect("Failed to open BAM file for read length sampling");
    
    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    
    for (_i, result) in bam.records().enumerate().take(sample_size) {
        if let Ok(record) = result {
            let read_len = record.seq().len();
            *length_counts.entry(read_len).or_insert(0) += 1;
        }
    }
    
    // Find mode of read lengths
    let mode_length = length_counts.iter()
                                    .max_by_key(|&(_, count)| count)
                                    .map(|(&key, _)| key)
                                    .unwrap_or(0);

    eprintln!("Computed mode read length from sample: {}", mode_length);
    mode_length
}

fn adjust_methylation_base(
    read_base: char,
    ref_base: char,
    read_num: u8,
    is_methylation: bool,
    cpg_only: bool,
    ref_seq: &[u8],
    genome_pos: usize
) -> char {
    if !is_methylation {
        return read_base;
    }
    
    if cpg_only {
        // Only convert in CpG context
        match (read_num, ref_base, read_base) {
            // Read 1: C in ref, next base is G, T in read -> treat as C (unmethylated CpG)
            (1, 'C', 'T') => {
                if genome_pos + 1 < ref_seq.len() {
                    let next_base = ref_seq[genome_pos + 1];
                    if next_base == b'G' || next_base == b'g' {
                        'C'
                    } else {
                        read_base
                    }
                } else {
                    read_base
                }
            }
            // Read 2: G in ref, next base is C, A in read -> treat as G (unmethylated CpG on reverse)
            (2, 'G', 'A') => {
                if genome_pos + 1 < ref_seq.len() {
                    let next_base = ref_seq[genome_pos + 1];
                    if next_base == b'C' || next_base == b'c' {
                        'G'
                    } else {
                        read_base
                    }
                } else {
                    read_base
                }
            }
            _ => read_base,
        }
    } else {
        // Convert all C>T and G>A
        match (read_num, ref_base, read_base) {
            // Read 1: C in ref, T in read -> treat as C (unmethylated C)
            (1, 'C', 'T') => 'C',
            // Read 2: G in ref, A in read -> treat as G (unmethylated C on reverse strand)
            (2, 'G', 'A') => 'G',
            _ => read_base,
        }
    }
}

fn create_mismatch_key(
    read_base: char,
    ref_base: char,
    r_pos: usize,
    seq_len: usize,
    is_reverse: bool,
    read_num: u8,
    is_methylation: bool,
    cpg_only: bool,
    ref_seq: &[u8],
    genome_pos: usize,
    mode_len: usize
) -> MismatchKey {
    // Apply reverse complement if read is mapped to reverse strand
    let strand_adjusted_read_base = if is_reverse { complement(read_base) } else { read_base };
    let strand_adjusted_ref_base = if is_reverse { complement(ref_base) } else { ref_base };

    // Apply methylation-aware base conversion if enabled
    let meth_and_strand_adjusted_read_base = adjust_methylation_base(
        strand_adjusted_read_base, 
        strand_adjusted_ref_base, 
        read_num, 
        is_methylation, 
        cpg_only, 
        ref_seq, 
        genome_pos
    );
    let adjusted_r_pos = if is_reverse {seq_len - r_pos -1} else {r_pos};
    let mode_adjusted_r_pos = correct_read_len_with_mode(adjusted_r_pos, seq_len, mode_len);

    MismatchKey {
        mismatch_type: format!("{}>{}", strand_adjusted_ref_base, meth_and_strand_adjusted_read_base),
        read_position: mode_adjusted_r_pos,
        read_num,
    }
}

fn compare_and_count(
    seq: &rust_htslib::bam::record::Seq,
    qual: &[u8],
    ref_seq: &[u8],
    r_pos: usize,
    genome_pos: usize,
    is_reverse: bool,
    read_num: u8,
    min_base_quality: u8,
    is_methylation: bool,
    cpg_only: bool,
    local_counts: &mut HashMap<MismatchKey, usize>,
    mode_len: usize
) {
    let seq_len = seq.len();
    let ref_len = ref_seq.len();
    
    if r_pos >= seq_len || genome_pos >= ref_len { return; }
    
    // Check base quality score (skip if below threshold)
    if r_pos >= qual.len() || qual[r_pos] < min_base_quality { return; }
    
    let Some(read_base) = base_to_char(seq[r_pos]) else { return; };    
    let Some(ref_base) = base_to_char(ref_seq[genome_pos]) else { return; };
    
    // Count both matches and mismatches
    let key = create_mismatch_key(read_base, ref_base, r_pos, seq_len, is_reverse, read_num, is_methylation, cpg_only, ref_seq, genome_pos, mode_len);
    *local_counts.entry(key).or_insert(0) += 1;
}

fn correct_read_len_with_mode(
    read_pos: usize,
    seq_len: usize,
    mode_len: usize
) -> usize {
    if mode_len > 0 && read_pos > seq_len/2 {  
        return read_pos + (mode_len - seq_len);
    } else {
        return read_pos;
    }
}

// Calculate genomic end position from CIGAR
fn calculate_end_pos(start_pos: i64, cigar: &rust_htslib::bam::record::CigarStringView) -> i64 {
    use rust_htslib::bam::record::Cigar::*;
    let mut end = start_pos;
    for op in cigar.iter() {
        match op {
            Match(len) | Equal(len) | Diff(len) | Del(len) | RefSkip(len) => {
                end += *len as i64;
            }
            _ => {}
        }
    }
    end
}

// Check if two reads overlap genomically
fn get_overlap_region(read1: &ReadInfo, read2: &ReadInfo, record1: &Record, record2: &Record) -> Option<(i64, i64)> {
    if read1.tid != read2.tid {
        return None;
    }
    
    let end1 = calculate_end_pos(read1.pos, &record1.cigar());
    let end2 = calculate_end_pos(read2.pos, &record2.cigar());
    
    let overlap_start = read1.pos.max(read2.pos);
    let overlap_end = end1.min(end2);
    
    if overlap_start < overlap_end {
        Some((overlap_start, overlap_end))
    } else {
        None
    }
}

// Process overlap region between two reads
// Process both reads - MismatchKey has read_num so they're counted separately
// Also detect inconsistencies between read1 and read2
fn process_overlap_region(
    record1: &Record,
    record2: &Record,
    overlap_start: i64,
    overlap_end: i64,
    local_counts: &mut HashMap<MismatchKey, usize>,
    inconsistency_counts: &mut HashMap<InconsistencyKey, usize>,
    reference: &ReferenceGenome,
    tid_to_name: &HashMap<i32, String>,
    min_base_quality: u8,
    is_methylation: bool,
    cpg_only: bool,
    mode_len: usize,
) {
    // Build genome_pos -> (read_pos, base, qual) mapping for both reads
    let mut read1_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();
    let mut read2_map: HashMap<usize, (usize, u8, u8)> = HashMap::new();
    
    for (record, map) in [(record1, &mut read1_map), (record2, &mut read2_map)] {
        if record.is_unmapped() || record.mapq() <= min_map_quality {
            continue;
        }
        
        let ref_start = record.pos() as usize;
        let seq = record.seq();
        let qual = record.qual();
        let cigar = record.cigar();
        
        let mut read_pos = 0;
        let mut ref_pos = ref_start;
        
        use rust_htslib::bam::record::Cigar::*;
        for cigar_op in cigar.iter() {
            match cigar_op {
                Match(len) | Equal(len) | Diff(len) => {
                    for i in 0..*len {
                        let genome_pos = ref_pos + i as usize;
                        if genome_pos >= overlap_start as usize && genome_pos < overlap_end as usize {
                            let r_pos = read_pos + i as usize;
                            if r_pos < seq.len() && r_pos < qual.len() {
                                map.insert(genome_pos, (r_pos, seq[r_pos], qual[r_pos]));
                            }
                        }
                    }
                    read_pos += *len as usize;
                    ref_pos += *len as usize;
                }
                Ins(len) => { read_pos += *len as usize; }
                Del(len) | RefSkip(len) => { ref_pos += *len as usize; }
                _ => {}
            }
        }
    }
    
    // Compare bases at overlapping genomic positions
    for (genome_pos, (r1_pos, r1_base, r1_qual)) in read1_map.iter() {
        if let Some((r2_pos, r2_base, r2_qual)) = read2_map.get(genome_pos) {
            // Both reads cover this position - check for inconsistency
            if r1_qual >= &min_base_quality && r2_qual >= &min_base_quality {
                if r1_base != r2_base {
                    // Bases disagree - record inconsistency
                    let r1_char = base_to_char(*r1_base).unwrap_or('N');
                    let r2_char = base_to_char(*r2_base).unwrap_or('N');
                    
                    let key = InconsistencyKey {
                        discordance_type: format!("R1:{}_R2:{}", r1_char, r2_char),
                        read1_position: *r1_pos,
                        read2_position: *r2_pos,
                    };
                    *inconsistency_counts.entry(key).or_insert(0) += 1;
                }
            }
        }
    }
    
    // Process mismatch counts for both reads
    for record in [record1, record2] {
        if record.is_unmapped() || record.mapq() <= min_map_quality {
            continue;
        }
        
        let tid = record.tid();
        let Some(chr_name) = tid_to_name.get(&tid) else { continue };
        let Some(ref_seq) = reference.get(chr_name.as_str()) else { continue };
        let ref_start = record.pos() as usize;
        let seq = record.seq();
        let qual = record.qual();
        let cigar = record.cigar();
        let read_num = if record.is_first_in_template() { 1 } else if record.is_last_in_template() { 2 } else { 1 };
        
        let mut read_pos = 0;
        let mut ref_pos = ref_start;
        
        use rust_htslib::bam::record::Cigar::*;
        for cigar_op in cigar.iter() {
            match cigar_op {
                Match(len) | Equal(len) | Diff(len) => {
                    for i in 0..*len {
                        let genome_pos = ref_pos + i as usize;
                        
                        // Only count if in overlap region
                        if genome_pos >= overlap_start as usize && genome_pos < overlap_end as usize {
                            let r_pos = read_pos + i as usize;
                            compare_and_count(
                                &seq, &qual, ref_seq, r_pos, genome_pos,
                                record.is_reverse(), read_num, min_base_quality,
                                is_methylation, cpg_only, local_counts, mode_len
                            );
                        }
                    }
                    read_pos += *len as usize;
                    ref_pos += *len as usize;
                }
                Ins(len) => { read_pos += *len as usize; }
                Del(len) | RefSkip(len) => { ref_pos += *len as usize; }
                _ => {}
            }
        }
    }
}

// Calculate genomic end position from CIGAR
fn process_record(
    record: &Record, 
    local_counts: &mut HashMap<MismatchKey, usize>, 
    reference: &ReferenceGenome, 
    tid_to_name: &HashMap<i32, String>,
    softclip_threshold: f64,
    min_base_quality: u8,
    is_methylation: bool,
    cpg_only: bool,
    mode_len: usize
) {
    
    // Skip read
    if record.is_unmapped() 
        || record.is_secondary() 
        || record.is_supplementary() 
        || record.mapq() <= min_map_quality
    { return; }
    
    // Get chromosome name (target_id - defined from bam header) and other read info
    let tid = record.tid();
    if tid < 0 { return; } // unmapped    
    let Some(chr_name) = tid_to_name.get(&tid) else { return }; // This directly binds chr_name (extract from Option or return)
    let Some(ref_seq) = reference.get(chr_name.as_str()) else { return; };
    let ref_start = record.pos() as usize; // (0-based)
    let seq = record.seq();
    let qual = record.qual();
    let cigar = record.cigar();
    let read_num = if record.is_first_in_template() { 1 } else if record.is_last_in_template() { 2 } else { 1 };
    
    // Walk through CIGAR and compare read to reference
    let mut read_pos = 0;
    let mut ref_pos = ref_start;
    
    for cigar_op in cigar.iter() {
        use rust_htslib::bam::record::Cigar::*;
        match cigar_op {
            Match(len) | Equal(len) | Diff(len) => {
                // Compare read and reference bases
                for i in 0..*len {
                    let r_pos = read_pos + i as usize; //if mode_len > 0 { mode_len - read_pos + i as usize } else { read_pos + i as usize };
                    let genome_pos = ref_pos + i as usize;
                    compare_and_count(&seq, &qual, ref_seq, r_pos, genome_pos, record.is_reverse(), read_num, min_base_quality, is_methylation, cpg_only, local_counts, mode_len);
                }
                read_pos += *len as usize;
                ref_pos += *len as usize;
            }
            SoftClip(len) => {
                // Soft-clipped bases - only count mismatches if at least 66% of bases match reference
                let mut total_bases = 0;
                let mut matching_bases = 0;
                let mut temp_mismatch_keys = Vec::<MismatchKey>::new();
                let seq_len = seq.len();
                
                for i in 0..*len {
                    let r_pos = read_pos + i as usize;
                    let genome_pos = if read_pos == 0 { // beginning of read
                        if ref_pos >= (*len - i) as usize {
                            ref_pos - (*len - i) as usize
                        } else {
                            continue;
                        }
                    } else { // end of read
                        ref_pos + i as usize
                    };
                    
                    if r_pos >= seq_len || genome_pos >= ref_seq.len() {
                        continue;
                    }
                    
                    // Skip bases with low quality
                    if r_pos >= qual.len() || qual[r_pos] < min_base_quality {
                        continue;
                    }
                    
                    let read_base = match seq[r_pos] {
                        b'A' => 'A',
                        b'C' => 'C',
                        b'G' => 'G',
                        b'T' => 'T',
                        _ => continue,
                    };
                    
                    let ref_base = match ref_seq[genome_pos] {
                        b'A' | b'a' => 'A',
                        b'C' | b'c' => 'C',
                        b'G' | b'g' => 'G',
                        b'T' | b't' => 'T',
                        _ => continue,
                    };
                    
                    total_bases += 1;
                    if read_base == ref_base {
                        matching_bases += 1;
                    } else {
                        // Create MismatchKey using helper function
                        let key = create_mismatch_key(read_base, ref_base, r_pos, seq_len, record.is_reverse(), read_num, is_methylation, cpg_only, ref_seq, genome_pos, mode_len);
                        temp_mismatch_keys.push(key);
                    }
                }
                
                // Only add mismatches to counts if at least 66% (default) match
                if total_bases > 0 && (matching_bases as f64 / total_bases as f64) >= softclip_threshold {
                    for key in temp_mismatch_keys {
                        *local_counts.entry(key).or_insert(0) += 1;
                    }
                }
                read_pos += *len as usize;
            }
            Ins(len) => { read_pos += *len as usize; }
            Del(len) | RefSkip(len) => { ref_pos += *len as usize; }
            HardClip(_) | Pad(_) => {}
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    
    if args.len() < 3 {
        eprintln!("Usage: {} <bam_file> <reference_fasta> [--threads N] [--region-size N] [--softclip-threshold F] [--min-base-quality Q] [--methylation] [--cpg-only] [----use-read-len-mode] [--insert-position-mode]", args[0]);
        std::process::exit(1);
    }
    
    let bam_path = &args[1];
    let fasta_path = &args[2];
    
    // Parse optional arguments
    let mut region_size: usize = 10_000_000; // 10MB regions by default
    let mut num_threads: usize = 0; // 0 means use all available cores
    let mut softclip_threshold: f64 = 0.66; // 66% threshold for soft-clipped bases
    let mut min_base_quality: u8 = 20; // Minimum phred quality score
    let mut min_map_quality: u8 = 20; // Minimum mapping quality
    let mut is_methylation: bool = false; // Methylation-aware mode (for bisulfite/EM-seq)
    let mut cpg_only: bool = false; // Only apply methylation conversion in CpG context
    let mut mode_len = 0;
    let mut insert_position_mode = false; // Not implemented yet (Not sure about it yet...)

    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--threads" | "-t" => {
                if i + 1 < args.len() {
                    num_threads = args[i + 1].parse().unwrap_or(0);
                    i += 2;
                } else {
                    eprintln!("Error: --threads requires a value");
                    std::process::exit(1);
                }
            }
            "--region-size" | "-r" => {
                if i + 1 < args.len() {
                    region_size = args[i + 1].parse().unwrap_or(10_000_000);
                    i += 2;
                } else {
                    eprintln!("Error: --region-size requires a value");
                    std::process::exit(1);
                }
            }
            "--softclip-threshold" => {
                if i + 1 < args.len() {
                    softclip_threshold = args[i + 1].parse().unwrap_or(0.66);
                    i += 2;
                }
            }
            "--min-base-quality" | "-q" => {
                if i + 1 < args.len() {
                    min_base_quality = args[i + 1].parse().unwrap_or(min_base_quality);
                    i += 2;
                } else {
                    eprintln!("Error: --min-base-quality requires a value");
                    std::process::exit(1);
                }
            }
            "--min-map-quality" => {
                if i + 1 < args.len() {
                    min_map_quality = args[i + 1].parse().unwrap_or(min_map_quality);
                    i += 2;
                } else {
                    eprintln!("Error: --min-map-quality requires a value");
                    std::process::exit(1);
                }
            }
            "--methylation" | "-m" => {
                is_methylation = true;
                i += 1;
            }
            "--cpg-only" => {
                cpg_only = true;
                i += 1;
            }
            "--use-read-len-mode" => {
                mode_len = compute_read_len_mode_from_sample_bam(bam_path, 100000);
                eprintln!("Using mode read length: {} to renumber read positions", mode_len);
                i += 1;
            }
            "--insert-position-mode" => {
                eprintln!("Using --insert-position-mode. Not read position");
                insert_position_mode = true;
                i += 1;
            }
            _ => {
                eprintln!("Unknown argument: {}", args[i]);
                std::process::exit(1);
            }
        }
    }
    
    // Set thread pool size
    if num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .expect("Failed to set thread pool size");
    }
    
    // Load reference genome
    let reference = Arc::new(load_reference_genome(fasta_path));
    
    // Print first and last 100 characters of each contig for verification
    eprintln!("\n{}", "=".repeat(80));
    eprintln!("Checking contig sequences:");
    eprintln!("{}", "=".repeat(80));
    
    let mut contig_names: Vec<&String> = reference.keys().collect();
    contig_names.sort();
    
    eprintln!("\n{}", "=".repeat(80));
    eprintln!("Reading BAM file: {}", bam_path);
    eprintln!("Region size: {} bp", region_size);
    eprintln!("Using {} threads", rayon::current_num_threads());
    
    // Open the BAM file to read header
    let bam = Reader::from_path(bam_path)
        .expect("Failed to open BAM file");
    
    // Create thread-safe chromosome name mapping (tid -> target_id=chr_name)
    let tid_to_name: Arc<HashMap<i32, String>> = Arc::new(
        (0..bam.header().target_count())
            .map(|i| (i as i32, String::from_utf8_lossy(bam.header().tid2name(i as u32)).to_string()))
            .collect()
    );
    
    // Create regions to process in parallel
    let mut regions = Vec::new();
    for (&tid, chr_name) in tid_to_name.iter() {
        let chr_len = bam.header().target_len(tid as u32).unwrap() as i64;
        
        // Split chromosome into regions
        let mut start = 0i64;
        while start < chr_len {
            let end = std::cmp::min(start + region_size as i64, chr_len);
            regions.push((tid as i32, start, end, chr_name.clone()));
            start = end;
        }
    }
    
    eprintln!("Created {} regions to process", regions.len());
    drop(bam); // Close the initial BAM reader
    
    let total_count = AtomicUsize::new(0);
    let mismatch_counts: Arc<Mutex<HashMap<MismatchKey, usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let overlap_mismatch_counts: Arc<Mutex<HashMap<MismatchKey, usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let overlap_pairs_count = AtomicUsize::new(0);
    let inconsistency_counts: Arc<Mutex<HashMap<InconsistencyKey, usize>>> = Arc::new(Mutex::new(HashMap::new()));
    
    // Process regions in parallel
    let bam_path_arc = Arc::new(bam_path.to_string());
    regions.par_iter().for_each(|(tid, start, end, chr_name)| {
        let mut bam = IndexedReader::from_path(bam_path_arc.as_str())
            .expect("Failed to open BAM file");
        
        // Fetch records in this region
        if let Err(e) = bam.fetch(FetchDefinition::Region(*tid, *start, *end)) {
            eprintln!("Warning: Failed to fetch region {}:{}-{}: {}", chr_name, start, end, e);
            return;
        }
        
        let ref_clone = Arc::clone(&reference);
        let tid_clone = Arc::clone(&tid_to_name); // This because in some cases tid could be different within chunck (e.g. chimera)
        let overlap_counts_clone = Arc::clone(&overlap_mismatch_counts);
        let inconsistency_clone = Arc::clone(&inconsistency_counts);
        
        // Local counts for this region (no locking during processing)
        let mut region_counts = HashMap::new();
        let mut overlap_region_counts = HashMap::new();
        let mut inconsistency_region_counts = HashMap::new();
        let mut local_count = 0;
        let mut local_overlap_count = 0;
        
        // Cache for read pairs
        let mut read_cache: HashMap<Vec<u8>, (ReadInfo, Record)> = HashMap::new();
        
        for result in bam.records() {
            if let Ok(record) = result {
                // Only process reads that START in this region to avoid duplicates
                let read_start = record.pos();
                if read_start >= *start && read_start < *end {
                    // Process normally
                    process_record(
                        &record, &mut region_counts, &ref_clone, &tid_clone, softclip_threshold, 
                        min_base_quality, is_methylation, cpg_only, mode_len
                    );
                    local_count += 1;
                    
                    // Cache for overlap detection if paired
                    if record.is_paired() && !record.is_unmapped() && !record.is_mate_unmapped() {
                        let qname = record.qname().to_vec();
                        let read_info = ReadInfo {
                            qname: qname.clone(),
                            tid: record.tid(),
                            pos: record.pos(),
                            seq_len: record.seq().len(),
                            cigar: Vec::new(),
                            is_reverse: record.is_reverse(),
                            read_num: if record.is_first_in_template() { 1 } else if record.is_last_in_template() { 2 } else { 1 },
                        };
                        
                        // Check if mate already cached
                        if let Some((mate_info, mate_record)) = read_cache.remove(&qname) {
                            // Found mate - check for overlap
                            if let Some((overlap_start, overlap_end)) = get_overlap_region(&read_info, &mate_info, &record, &mate_record) {
                                local_overlap_count += 1;
                                // Process overlap region for both reads
                                process_overlap_region(
                                    &record, &mate_record, overlap_start, overlap_end,
                                    &mut overlap_region_counts, &mut inconsistency_region_counts,
                                    &ref_clone, &tid_clone,
                                    min_base_quality, is_methylation, cpg_only, mode_len
                                );
                            }
                        } else {
                            // Store for when mate arrives
                            read_cache.insert(qname, (read_info, record.clone()));
                        }
                    }
                }
            }
        }
        
        // Merge region counts into global counts (single lock per region)
        // global_counts: MutexGuard<HashMap<MismatchKey, usize>> 
        // *global_counts: HashMap<MismatchKey, usize> 
        if !region_counts.is_empty() {
            let mut global_counts = mismatch_counts.lock().unwrap(); 
            for (key, count) in region_counts {
                *global_counts.entry(key).or_insert(0) += count;
            }
        }
        
        // Merge overlap counts
        if !overlap_region_counts.is_empty() {
            let mut global_overlap_counts = overlap_counts_clone.lock().unwrap();
            for (key, count) in overlap_region_counts {
                *global_overlap_counts.entry(key).or_insert(0) += count;
            }
        }
        
        // Merge inconsistency counts
        if !inconsistency_region_counts.is_empty() {
            let mut global_inconsistency_counts = inconsistency_clone.lock().unwrap();
            for (key, count) in inconsistency_region_counts {
                *global_inconsistency_counts.entry(key).or_insert(0) += count;
            }
        }
        
        if local_overlap_count > 0 {
            overlap_pairs_count.fetch_add(local_overlap_count, Ordering::Relaxed);
        }
        
        let prev_total = total_count.fetch_add(local_count, Ordering::Relaxed);
        if (prev_total + local_count) / 100000 > prev_total / 100000 { // every 100k records (arbitrary)
            eprintln!("Processed {} records...", prev_total + local_count);
        }
    });
    
    eprintln!("\nTotal records processed: {}", total_count.load(Ordering::Relaxed));
    let overlap_pairs = overlap_pairs_count.load(Ordering::Relaxed);
    eprintln!("Read pairs with overlaps: {}", overlap_pairs);
    
    // Restructure data: (read_num, position) -> mismatch_type -> count
    let counts = mismatch_counts.lock().unwrap();
    let mut position_map: HashMap<(u8, usize), HashMap<String, usize>> = HashMap::new();
    
    for (key, count) in counts.iter() {
        position_map
            .entry((key.read_num, key.read_position))
            .or_insert_with(HashMap::new)
            .insert(key.mismatch_type.clone(), *count);
    }
    
    // Get all unique mismatch types and sort them
    let mut all_mismatch_types: Vec<String> = counts
        .keys()
        .map(|k| k.mismatch_type.clone())
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    all_mismatch_types.sort();
    
    // Get all (read_num, position) pairs and sort them
    let mut positions: Vec<(u8, usize)> = position_map.keys().copied().collect();
    positions.sort();
    
    // Print header
    // print!("{:<10} {:<10}", "Read", "Position");
    print!("Read,Position");
    for mismatch_type in &all_mismatch_types {
        // print!(" {:<10}", mismatch_type);
        print!(",{}", mismatch_type);
    }
    println!();
    
    // Print data rows
    for &(read_num, pos) in positions.iter() {
        // print!("{:<10} {:<10}", read_num, pos);
        print!("{},{}", read_num, pos);
        
        if let Some(mismatch_counts) = position_map.get(&(read_num, pos)) {
            for mismatch_type in &all_mismatch_types {
                let count = mismatch_counts.get(mismatch_type).unwrap_or(&0);
                // print!(" {:<10}", count);
                print!(",{}", count);
            }
        }
        println!();
    }
    
    eprintln!("\nTotal unique (read, position) combinations: {}", positions.len());
    eprintln!("Total unique mismatch types: {}", all_mismatch_types.len());
    
    // Print overlap table if there are overlaps
    if overlap_pairs > 0 {
        let overlap_counts = overlap_mismatch_counts.lock().unwrap();
        if !overlap_counts.is_empty() {
            eprintln!("\n{}", "=".repeat(80));
            eprintln!("OVERLAP REGION MISMATCHES");
            eprintln!("{}", "=".repeat(80));
            
            let mut overlap_position_map: HashMap<(u8, usize), HashMap<String, usize>> = HashMap::new();
            for (key, count) in overlap_counts.iter() {
                overlap_position_map
                    .entry((key.read_num, key.read_position))
                    .or_insert_with(HashMap::new)
                    .insert(key.mismatch_type.clone(), *count);
            }
            
            let mut overlap_mismatch_types: Vec<String> = overlap_counts
                .keys()
                .map(|k| k.mismatch_type.clone())
                .collect::<std::collections::HashSet<_>>()
                .into_iter()
                .collect();
            overlap_mismatch_types.sort();
            
            let mut overlap_positions: Vec<(u8, usize)> = overlap_position_map.keys().copied().collect();
            overlap_positions.sort();
            
            println!("\n# Overlap Region Data");
            print!("Read,Position");
            for mismatch_type in &overlap_mismatch_types {
                print!(",{}", mismatch_type);
            }
            println!();
            
            for &(read_num, pos) in overlap_positions.iter() {
                print!("{},{}", read_num, pos);
                if let Some(mismatch_counts_map) = overlap_position_map.get(&(read_num, pos)) {
                    for mismatch_type in &overlap_mismatch_types {
                        let count = mismatch_counts_map.get(mismatch_type).unwrap_or(&0);
                        print!(",{}", count);
                    }
                }
                println!();
            }
            
            eprintln!("\nTotal unique overlap (read, position) combinations: {}", overlap_positions.len());
            eprintln!("{}", "=".repeat(80));
        }
        
        // Print inconsistencies
        let incons_counts = inconsistency_counts.lock().unwrap();
        if !incons_counts.is_empty() {
            eprintln!("\n{}", "=".repeat(80));
            eprintln!("READ PAIR INCONSISTENCIES IN OVERLAP REGIONS");
            eprintln!("{}", "=".repeat(80));
            
            // Group by read positions
            let mut incons_position_map: HashMap<(usize, usize), HashMap<String, usize>> = HashMap::new();
            for (key, count) in incons_counts.iter() {
                incons_position_map
                    .entry((key.read1_position, key.read2_position))
                    .or_insert_with(HashMap::new)
                    .insert(key.discordance_type.clone(), *count);
            }
            
            let mut incons_types: Vec<String> = incons_counts
                .keys()
                .map(|k| k.discordance_type.clone())
                .collect::<std::collections::HashSet<_>>()
                .into_iter()
                .collect();
            incons_types.sort();
            
            let mut incons_positions: Vec<(usize, usize)> = incons_position_map.keys().copied().collect();
            incons_positions.sort();
            
            println!("\n# Read Pair Inconsistencies");
            print!("Read1_Pos,Read2_Pos");
            for incons_type in &incons_types {
                print!(",{}", incons_type);
            }
            println!();
            
            for &(r1_pos, r2_pos) in incons_positions.iter() {
                print!("{},{}", r1_pos, r2_pos);
                if let Some(incons_map) = incons_position_map.get(&(r1_pos, r2_pos)) {
                    for incons_type in &incons_types {
                        let count = incons_map.get(incons_type).unwrap_or(&0);
                        print!(",{}", count);
                    }
                }
                println!();
            }
            
            let total_inconsistencies: usize = incons_counts.values().sum();
            eprintln!("\nTotal inconsistencies: {}", total_inconsistencies);
            eprintln!("Unique position pairs with inconsistencies: {}", incons_positions.len());
            eprintln!("{}", "=".repeat(80));
        }
    }
}
