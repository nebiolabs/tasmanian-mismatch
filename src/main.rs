use rust_htslib::bam::{Read, Reader, IndexedReader, Record, FetchDefinition};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::collections::HashMap;
use std::sync::Arc;
use std::fs::File;
use std::io::{Write, BufWriter};
use clap::Parser;
// Import from library
use rustmanian_mismatch::*;

/// Tasmanian Mismatch - A tool for analyzing mismatches in BAM files
#[derive(Parser, Debug)]
#[command(name = "rustmanian-mismatch")]
#[command(version, about, long_about = None)]
struct Args {
    bam_file: String, //path::PathBuf
    reference_fasta: String, //path::PathBuf
    
    /// Number of threads to use (0 = all available cores)
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
    
    /// Region size in base pairs for parallel processing
    #[arg(short = 'r', long, default_value_t = 10_000_000)]
    region_size: usize,
    
    /// At least that fraction of bases in softclip correspond to reference to be considered mapped
    #[arg(long, default_value_t = 0.66)]
    softclip_threshold: f64,
    
    #[arg(short = 'q', long, default_value_t = 20)]
    min_base_quality: u8,
    
    #[arg(long, default_value_t = 10)]
    min_map_quality: u8,
    
    /// All C/T in read 1 are converted back to C (for bisulfite/EM-seq)
    #[arg(short = 'm', long)]
    methylation: bool,
    
    /// C/T in CpG context in read1 are converted back to C
    #[arg(long)]
    cpg_only: bool,
    
    /// Use read length mode to renumber read positions (can stack start and end of reads)
    #[arg(long)]
    use_read_len_mode: bool,
    
    /// Use insert position mode (not read position)
    #[arg(long)]
    insert_position_mode: bool,
    
    /// Genomic thresholds for calling a potential variant
    #[arg(long, default_value_t = 7)]
    genomic_threshold: usize,
    #[arg(long, default_value_t = 10)]
    genomic_depth_threshold: usize,
    
    /// Next 3 are SAM flag filters (same as samtools)
    #[arg(short = 'f', default_value_t = 0)]
    required_flags: u16,
    #[arg(short = 'F', default_value_t = 0)]
    filter_flags: u16,
    #[arg(short = 'G', default_value_t = 0)]
    excl_flags: u16,
}


fn main() {
    let args = Args::parse();
    
    let bam_path = &args.bam_file;
    let fasta_path = &args.reference_fasta;
    let region_size = args.region_size;
    let num_threads = args.threads;
    let softclip_threshold = args.softclip_threshold;
    let min_base_quality = args.min_base_quality;
    let min_map_quality = args.min_map_quality;
    let is_methylation = args.methylation;
    let cpg_only = args.cpg_only;
    let genomic_threshold = args.genomic_threshold;
    let required_flags = args.required_flags;
    let filter_flags = args.filter_flags;
    let excl_flags = args.excl_flags;
    let _insert_position_mode = args.insert_position_mode;
    let genomic_depth_threshold = args.genomic_depth_threshold;
    let mode_len = if args.use_read_len_mode {
        let mode = compute_read_len_mode_from_sample_bam(bam_path, 100000);
        eprintln!("Using mode read length: {} to renumber read positions", mode);
        mode
    } else {
        0
    };
    
    if _insert_position_mode {
        eprintln!("Using --insert-position-mode. Not read position");
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
    let genomic_mismatch_counts: Arc<Mutex<HashMap<GenomicMismatchKey, (usize, usize)>>> = Arc::new(Mutex::new(HashMap::new()));

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

        // for variant "cheap" identification
        let mut genomic_region_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> = HashMap::new();
        let mut genomic_position_depth: HashMap<i64, usize> = HashMap::new(); // genomic_position -> total count of reads covering position (i64 not u64 for compatibility later. same usize)

        // Cache for read pairs
        let mut read_cache: HashMap<Vec<u8>, (ReadInfo, Record)> = HashMap::new();
        
        for result in bam.records() {
            if let Ok(record) = result {
                // Only process reads that START in this region to avoid duplicates
                let read_start = record.pos();
                if read_start >= *start && read_start < *end {
                    // Process normally - update both region_counts and genomic_region_counts
                    process_record(
                        &record, &mut region_counts, Some(&mut genomic_region_counts), &ref_clone, &tid_clone, softclip_threshold, 
                        min_base_quality, is_methylation, cpg_only, mode_len, min_map_quality, Some(&mut genomic_position_depth),
                        required_flags, filter_flags, excl_flags
                    );
                    local_count += 1;
        

                    // Cache for overlap detection if paired
                    if record.is_paired() && !record.is_unmapped() && !record.is_mate_unmapped() {
                        let qname = record.qname().to_vec();
                        let read_info = ReadInfo {
                            tid: record.tid(),
                            pos: record.pos(),
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
                                    min_base_quality, is_methylation, cpg_only, mode_len, min_map_quality,
                                    required_flags, filter_flags, excl_flags
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
        
        // Filter region_counts based on genomic_region_counts threshold 
        // and merge values to global genomic counts (only if >= threshold).
        if !genomic_region_counts.is_empty() {
            let mut global_genomic_counts = genomic_mismatch_counts.lock().unwrap();

            for (genomic_key, genomic_value) in genomic_region_counts.iter_mut() {
                let genomic_depth = genomic_position_depth.get(&genomic_key.genomic_position).cloned().unwrap_or(0);

                if genomic_value.count >= genomic_threshold && genomic_depth >= genomic_depth_threshold{
                    // Update global genomic counts only for positions above threshold
                    *global_genomic_counts.entry(genomic_key.clone()).or_insert((0, 0))  = (
                        genomic_value.count, 
                        genomic_depth // single thread -> single contig.
                    );
                    
                    // Decrement count for each mismatch_key that contributed to this genomic position
                    for mismatch_key in &genomic_value.mismatch_keys {
                        match region_counts.get_mut(mismatch_key) {
                            Some(count) => *count -= 1,  // Will underflow/panic in debug if count is 0
                            None => panic!("MismatchKey {:?} not found in region_counts", mismatch_key),
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
    
    // Write genomic mismatches to potential_variants.tsv and report statistics
    {
        let genomic_counts = genomic_mismatch_counts.lock().unwrap();
        eprintln!("Genomic positions with {} or more mismatches: {}", genomic_threshold, genomic_counts.len());
        
        // Write to file
        let file = File::create("potential_variants.tsv").expect("Failed to create potential_variants.tsv");
        let mut writer = BufWriter::new(file);

        writeln!(writer, "chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth").expect("Failed to write header");
        
        //let mut sorted_keys: Vec<_> = genomic_counts.iter().collect();
        // sorted_keys.sort_by_key(|(k, _)| (&k.chromosome, k.genomic_position)); // This takes forever for large datasets
        let sorted_keys: Vec<_> = genomic_counts.iter().collect(); // unsorted for speed -> can change later to the above lines.

        for (key, count) in sorted_keys {
            // Parse mismatch_type "REF>ALT" to extract ref and alt bases
            let parts: Vec<&str> = key.mismatch_type.split('>').collect();
            if parts.len() == 2 {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    key.chromosome,
                    key.genomic_position,
                    parts[0],  // reference base
                    parts[1],  // mismatch base
                    count.0,
                    count.1
                ).expect("Failed to write to file");
            }
        }
        
        eprintln!("Wrote {} genomic positions to potential_variants.tsv", genomic_counts.len());
    }
    
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