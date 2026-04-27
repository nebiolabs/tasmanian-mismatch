use clap::Parser;
use rayon::prelude::*;
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Reader, Record};
use rustmanian_mismatch::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

#[derive(Parser, Debug)]
#[command(name = "tasmanian-diagnostics")]
#[command(version, about = "Collect genomic variants and read-pair overlap inconsistencies")]
struct Args {
    bam_file: String,
    reference_fasta: String,

    /// Number of threads to use (0 = all available cores)
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,

    /// Region size in base pairs for parallel processing
    #[arg(short = 'r', long, default_value_t = 10_000_000)]
    region_size: usize,

    /// At least that fraction of soft-clipped bases must match reference
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

    /// Optional read-length normalization argument preserved for processing compatibility
    #[arg(long, num_args = 0..=1, value_name = "LENGTH")]
    use_read_len_max: Option<Option<u32>>,

    /// Keep processing behavior compatible with split-read vs insert-position mode
    #[arg(long, default_value_t = false)]
    use_insert_mode: bool,

    /// Minimum mismatch count for reporting a genomic site
    #[arg(long, default_value_t = 7)]
    genomic_threshold: usize,

    /// Minimum depth for reporting a genomic site
    #[arg(long, default_value_t = 10)]
    genomic_depth_threshold: usize,

    /// SAM flag filters
    #[arg(short = 'f', default_value_t = 0)]
    required_flags: u16,

    #[arg(short = 'F', default_value_t = 0)]
    filter_flags: u16,

    #[arg(short = 'G', default_value_t = 0)]
    excl_flags: u16,

    /// BED file with regions to filter/mask (optional)
    #[arg(short = 'b', long)]
    bed_file: Option<String>,

    /// Filter mode: 'mask' (skip individual bases) or 'filter' (skip whole reads)
    #[arg(long, default_value = "mask")]
    bed_filter_mode: String,

    /// Output path for genomic potential variants table
    #[arg(long, default_value = "potential_variants.tsv")]
    variants_output: String,

    /// Output path for overlap inconsistencies table
    #[arg(long, default_value = "read_pair_inconsistencies.tsv")]
    inconsistencies_output: String,

    /// Output path for mismatch-key discounts to be consumed by main.rs
    #[arg(long, default_value = "variant_discounts.tsv")]
    discount_output: String,
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("Failed to set thread pool size");
    }

    let mode_len = match args.use_read_len_max {
        None => 0,
        Some(None) => {
            let max_l = compute_read_len_max_from_sample_bam(&args.bam_file, 100_000);
            eprintln!("Using max read length: {} for processing compatibility", max_l);
            max_l
        }
        Some(Some(val)) => {
            eprintln!("Using max read length: {} for processing compatibility", val);
            val as usize
        }
    };

    let mut reference = load_reference_genome(&args.reference_fasta);

    let bed_for_filtering = if let Some(bed_path) = args.bed_file.as_ref() {
        let regions = parse_bed_file(bed_path)
            .expect("Failed to load BED file. Please check the file path and format.");

        match args.bed_filter_mode.as_str() {
            "filter" => Some(Arc::new(regions)),
            "mask" => {
                let masked_bases = mask_reference_with_bed(&mut reference, &regions);
                eprintln!("Masked {} bases in reference genome", masked_bases);
                None
            }
            _ => {
                panic!(
                    "Invalid --bed-filter-mode '{}'. Expected 'mask' or 'filter'.",
                    args.bed_filter_mode
                );
            }
        }
    } else {
        None
    };

    let reference = Arc::new(reference);
    let bam = Reader::from_path(&args.bam_file).expect("Failed to open BAM file");

    let tid_to_name: Arc<HashMap<i32, String>> = Arc::new(
        (0..bam.header().target_count())
            .map(|i| {
                (
                    i as i32,
                    String::from_utf8_lossy(bam.header().tid2name(i as u32)).to_string(),
                )
            })
            .collect(),
    );

    let mut regions = Vec::new();
    for (&tid, chr_name) in tid_to_name.iter() {
        let chr_len = bam.header().target_len(tid as u32).unwrap_or(0) as i64;
        let mut start = 0i64;
        while start < chr_len {
            let end = std::cmp::min(start + args.region_size as i64, chr_len);
            regions.push((tid, start, end, chr_name.clone()));
            start = end;
        }
    }
    drop(bam);

    eprintln!("Created {} regions to process", regions.len());

    let total_count = AtomicUsize::new(0);
    let overlap_pairs_count = AtomicUsize::new(0);
    let genomic_mismatch_counts: Arc<Mutex<HashMap<GenomicMismatchKey, (usize, usize)>>> =
        Arc::new(Mutex::new(HashMap::new()));
    let inconsistency_counts: Arc<Mutex<HashMap<InconsistencyKey, usize>>> =
        Arc::new(Mutex::new(HashMap::new()));
    let mismatch_discounts: Arc<Mutex<HashMap<MismatchKey, usize>>> =
        Arc::new(Mutex::new(HashMap::new()));

    let processing_config = ProcessingConfig {
        softclip_threshold: args.softclip_threshold,
        min_base_quality: args.min_base_quality,
        is_methylation: args.methylation,
        cpg_only: args.cpg_only,
        mode_len,
        min_map_quality: args.min_map_quality,
        required_flags: args.required_flags,
        filter_flags: args.filter_flags,
        excl_flags: args.excl_flags,
        use_insert_mode: args.use_insert_mode,
    };

    let bam_path_arc = Arc::new(args.bam_file.clone());

    regions.par_iter().for_each(|(tid, start, end, chr_name)| {
        let mut bam =
            IndexedReader::from_path(bam_path_arc.as_str()).expect("Failed to open BAM file");

        if let Err(e) = bam.fetch(FetchDefinition::Region(*tid, *start, *end)) {
            eprintln!(
                "Warning: Failed to fetch region {}:{}-{}: {}",
                chr_name, start, end, e
            );
            return;
        }

        let ref_clone = Arc::clone(&reference);
        let tid_clone = Arc::clone(&tid_to_name);
        let inconsistency_clone = Arc::clone(&inconsistency_counts);
        let discount_clone = Arc::clone(&mismatch_discounts);

        let chunk_bed_intervals = if let Some(bed) = bed_for_filtering.as_ref() {
            filter_bed_for_region(bed, chr_name, *start, *end)
        } else {
            Vec::new()
        };

        let processing_context = ProcessingContext {
            reference: &ref_clone,
            tid_to_name: &tid_clone,
            bed_intervals: &chunk_bed_intervals,
        };

        let mut local_count = 0usize;
        let mut local_overlap_count = 0usize;

        let mut dummy_mismatch_counts: HashMap<MismatchKey, usize> = HashMap::new();
        let mut dummy_overlap_mismatch_counts: HashMap<MismatchKey, usize> = HashMap::new();
        let mut inconsistency_region_counts: HashMap<InconsistencyKey, usize> = HashMap::new();
        let mut genomic_region_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> =
            HashMap::new();
        let mut genomic_position_depth: HashMap<i64, usize> = HashMap::new();

        let mut read_groups: HashMap<Vec<u8>, Vec<Record>> = HashMap::new();

        for result in bam.records() {
            let Ok(record) = result else {
                continue;
            };

            let read_start = record.pos();
            let starts_in_chunk = read_start >= *start && read_start < *end;
            if !starts_in_chunk {
                continue;
            }

            let should_skip_whole_read = args.bed_filter_mode == "filter"
                && !chunk_bed_intervals.is_empty()
                && {
                    let read_end = calculate_end_pos(record.pos(), &record.cigar());
                    chunk_bed_intervals
                        .iter()
                        .any(|interval| record.pos() <= interval.end && read_end >= interval.start)
                };
            if should_skip_whole_read {
                continue;
            }

            let has_mapped_pair =
                record.is_paired() && !record.is_unmapped() && !record.is_mate_unmapped();
            if has_mapped_pair {
                let qname = record.qname().to_vec();
                read_groups.entry(qname).or_insert_with(Vec::new).push(record);
                continue;
            }

            process_single_record(
                &record,
                &mut dummy_mismatch_counts,
                &mut genomic_region_counts,
                &mut genomic_position_depth,
                &processing_context,
                processing_config,
            );
            local_count += 1;
        }

        for (_qname, records) in read_groups {
            if records.len() == 2 {
                let record1 = &records[0];
                let record2 = &records[1];

                let info1 = ReadInfo {
                    tid: record1.tid(),
                    pos: record1.pos(),
                };
                let info2 = ReadInfo {
                    tid: record2.tid(),
                    pos: record2.pos(),
                };

                if let Some((overlap_start, overlap_end)) =
                    get_overlap_region(&info1, &info2, record1, record2)
                {
                    local_overlap_count += 1;
                    process_paired_reads_with_overlap(
                        record1,
                        record2,
                        overlap_start,
                        overlap_end,
                        &mut dummy_mismatch_counts,
                        &mut dummy_overlap_mismatch_counts,
                        &mut inconsistency_region_counts,
                        Some(&mut genomic_region_counts),
                        Some(&mut genomic_position_depth),
                        &processing_context,
                        processing_config,
                    );
                    local_count += 2;
                } else {
                    for record in &records {
                        process_single_record(
                            record,
                            &mut dummy_mismatch_counts,
                            &mut genomic_region_counts,
                            &mut genomic_position_depth,
                            &processing_context,
                            processing_config,
                        );
                        local_count += 1;
                    }
                }
            } else {
                for record in &records {
                    process_single_record(
                        record,
                        &mut dummy_mismatch_counts,
                        &mut genomic_region_counts,
                        &mut genomic_position_depth,
                        &processing_context,
                        processing_config,
                    );
                    local_count += 1;
                }
            }
        }

        if !genomic_region_counts.is_empty() {
            let mut global_genomic_counts = genomic_mismatch_counts.lock().expect("Lock poisoned");
            let mut global_discount_counts = discount_clone.lock().expect("Lock poisoned");

            for (genomic_key, genomic_value) in genomic_region_counts {
                let genomic_depth = genomic_position_depth
                    .get(&genomic_key.genomic_position)
                    .copied()
                    .unwrap_or(0);

                if genomic_value.count >= args.genomic_threshold
                    && genomic_depth >= args.genomic_depth_threshold
                {
                    *global_genomic_counts.entry(genomic_key).or_insert((0, 0)) =
                        (genomic_value.count, genomic_depth);

                    for mismatch_key in genomic_value.mismatch_keys {
                        *global_discount_counts.entry(mismatch_key).or_insert(0) += 1;
                    }
                }
            }
        }

        if !inconsistency_region_counts.is_empty() {
            let mut global_inconsistency_counts = inconsistency_clone.lock().expect("Lock poisoned");
            for (key, count) in inconsistency_region_counts {
                *global_inconsistency_counts.entry(key).or_insert(0) += count;
            }
        }

        if local_overlap_count > 0 {
            overlap_pairs_count.fetch_add(local_overlap_count, Ordering::Relaxed);
        }

        let prev_total = total_count.fetch_add(local_count, Ordering::Relaxed);
        if (prev_total + local_count) / 100_000 > prev_total / 100_000 {
            eprintln!("Processed {} records...", prev_total + local_count);
        }
    });

    eprintln!(
        "Total records processed: {}",
        total_count.load(Ordering::Relaxed)
    );
    eprintln!(
        "Read pairs with overlaps: {}",
        overlap_pairs_count.load(Ordering::Relaxed)
    );

    {
        let genomic_counts = genomic_mismatch_counts.lock().expect("Lock poisoned");
        write_potential_variants_tsv(&genomic_counts, &args.variants_output)
            .expect("Failed to write variants output");
        eprintln!("Wrote variants table to {}", args.variants_output);
    }

    {
        let incons_counts = inconsistency_counts.lock().expect("Lock poisoned");
        write_inconsistencies_tsv(&incons_counts, &args.inconsistencies_output)
            .expect("Failed to write inconsistencies output");
        eprintln!(
            "Wrote {} inconsistency keys to {}",
            incons_counts.len(),
            args.inconsistencies_output
        );
    }

    {
        let discounts = mismatch_discounts.lock().expect("Lock poisoned");
        write_mismatch_discounts_tsv(&discounts, &args.discount_output)
            .expect("Failed to write mismatch discount output");
        eprintln!(
            "Wrote {} mismatch discount keys to {}",
            discounts.len(),
            args.discount_output
        );
    }
}
