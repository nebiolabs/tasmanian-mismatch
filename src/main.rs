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
    genome_pos: usize
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

    MismatchKey {
        mismatch_type: format!("{}>{}", strand_adjusted_ref_base, meth_and_strand_adjusted_read_base),
        read_position: adjusted_r_pos,
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
    local_counts: &mut HashMap<MismatchKey, usize>
) {
    let seq_len = seq.len();
    let ref_len = ref_seq.len();
    
    if r_pos >= seq_len || genome_pos >= ref_len { return; }
    
    // Check base quality score (skip if below threshold)
    if r_pos >= qual.len() || qual[r_pos] < min_base_quality { return; }
    
    let Some(read_base) = base_to_char(seq[r_pos]) else { return; };    
    let Some(ref_base) = base_to_char(ref_seq[genome_pos]) else { return; };
    
    // Count both matches and mismatches
    let key = create_mismatch_key(read_base, ref_base, r_pos, seq_len, is_reverse, read_num, is_methylation, cpg_only, ref_seq, genome_pos);
    *local_counts.entry(key).or_insert(0) += 1;
}

fn process_record(
    record: &Record, 
    local_counts: &mut HashMap<MismatchKey, usize>, 
    reference: &ReferenceGenome, 
    tid_to_name: &HashMap<i32, String>,
    softclip_threshold: f64,
    min_base_quality: u8,
    is_methylation: bool,
    cpg_only: bool
) {
    
    // Skip read
    if record.is_unmapped() 
        || record.is_secondary() 
        || record.is_supplementary() 
        || record.mapq() <= 20 
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
                    let r_pos = read_pos + i as usize;
                    let genome_pos = ref_pos + i as usize;
                    compare_and_count(&seq, &qual, ref_seq, r_pos, genome_pos, record.is_reverse(), read_num, min_base_quality, is_methylation, cpg_only, local_counts);
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
                        let key = create_mismatch_key(read_base, ref_base, r_pos, seq_len, record.is_reverse(), read_num, is_methylation, cpg_only, ref_seq, genome_pos);
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
        eprintln!("Usage: {} <bam_file> <reference_fasta> [--threads N] [--region-size N] [--softclip-threshold F] [--min-base-quality Q] [--methylation] [--cpg-only]", args[0]);
        std::process::exit(1);
    }
    
    let bam_path = &args[1];
    let fasta_path = &args[2];
    
    // Parse optional arguments
    let mut region_size: usize = 10_000_000; // 10MB regions by default
    let mut num_threads: usize = 0; // 0 means use all available cores
    let mut softclip_threshold: f64 = 0.66; // 66% threshold for soft-clipped bases
    let mut min_base_quality: u8 = 20; // Minimum phred quality score
    let mut is_methylation: bool = false; // Methylation-aware mode (for bisulfite/EM-seq)
    let mut cpg_only: bool = false; // Only apply methylation conversion in CpG context
    
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
                    min_base_quality = args[i + 1].parse().unwrap_or(20);
                    i += 2;
                } else {
                    eprintln!("Error: --min-base-quality requires a value");
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
    
    // Create thread-safe chromosome name mapping
    let tid_to_name: Arc<HashMap<i32, String>> = Arc::new(
        (0..bam.header().target_count())
            .map(|i| (i as i32, String::from_utf8_lossy(bam.header().tid2name(i as u32)).to_string()))
            .collect()
    );
    
    // Create regions to process in parallel
    let mut regions = Vec::new();
    for tid in 0..bam.header().target_count() {
        let chr_len = bam.header().target_len(tid).unwrap() as i64;
        let chr_name = String::from_utf8_lossy(bam.header().tid2name(tid)).to_string();
        
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
        let tid_clone = Arc::clone(&tid_to_name);
        
        // Local counts for this region (no locking during processing)
        let mut region_counts = HashMap::new();
        let mut local_count = 0;
        
        for result in bam.records() {
            if let Ok(record) = result {
                // Only process reads that START in this region to avoid duplicates
                let read_start = record.pos();
                if read_start >= *start && read_start < *end {
                    process_record(&record, &mut region_counts, &ref_clone, &tid_clone, softclip_threshold, min_base_quality, is_methylation, cpg_only);
                    local_count += 1;
                }
            }
        }
        
        // Merge region counts into global counts (single lock per region)
        if !region_counts.is_empty() {
            let mut global_counts = mismatch_counts.lock().unwrap();
            for (key, count) in region_counts {
                *global_counts.entry(key).or_insert(0) += count;
            }
        }
        
        let prev_total = total_count.fetch_add(local_count, Ordering::Relaxed);
        if (prev_total + local_count) / 100000 > prev_total / 100000 {
            eprintln!("Processed {} records...", prev_total + local_count);
        }
    });
    
    eprintln!("\nTotal records processed: {}", total_count.load(Ordering::Relaxed));
    
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
}
