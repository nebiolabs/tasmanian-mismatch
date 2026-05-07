use clap::Parser;
use rayon::prelude::*;
use rustmanian_mismatch::{
    apply_external_discounts, build_tid_map_and_regions, compute_read_len_max_from_sample_bam,
    configure_thread_pool, load_discount_table, load_reference_genome, mask_reference_with_bed,
    maybe_parse_bed_file, process_region, write_normalized_output, write_output,
    write_rescaling_matrix_output, Args, InsertKey, PositionMode,
};
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

fn main() {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();
    log::info!(
        "position mode {:?}, overlap mode {:?}",
        args.position_mode,
        args.overlap_mode
    );
    let bed_filter_whole_reads = args.bed_filter_mode == "filter";

    configure_thread_pool(args.threads);

    let (max_read_len, mut reference, bed_regions) = std::thread::scope(|s| {
        let t1 = s.spawn(|| compute_read_len_max_from_sample_bam(&args.bam_path, 10_000));
        let t2 = s.spawn(|| load_reference_genome(&args.reference_path));
        let t3 = s.spawn(|| maybe_parse_bed_file(args.bed_file.as_deref()));
        (t1.join().unwrap(), t2.join().unwrap(), t3.join().unwrap())
    });
    log::info!("Sampled max read length: {}", max_read_len);

    let bed_for_filtering = if let Some(regions) = bed_regions {
        log::info!("BED filter mode: {}", args.bed_filter_mode);
        match args.bed_filter_mode.as_str() {
            "filter" => {
                log::info!("Keeping BED regions for whole-read overlap filtering...");
                Some(Arc::new(regions))
            }
            "mask" => {
                log::info!("Masking reference genome at BED regions...");
                let masked_bases = mask_reference_with_bed(&mut reference, &regions);
                log::info!("Masked {} bases in reference genome", masked_bases);
                None
            }
            _ => panic!(
                "Invalid --bed-filter-mode '{}'. Expected 'mask' or 'filter'.",
                args.bed_filter_mode
            ),
        }
    } else {
        None
    };

    let reference = Arc::new(reference);
    let (tid_to_name, regions) = build_tid_map_and_regions(&args.bam_path, args.region_size);
    log::info!("Processing {} indexed regions", regions.len());

    let total_reads = AtomicUsize::new(0);
    let global_counts: Arc<Mutex<HashMap<InsertKey, usize>>> = Arc::new(Mutex::new(HashMap::new()));

    regions.par_iter().for_each(|(tid, start, end)| {
        let Some(chr_name) = tid_to_name.get(tid) else {
            return;
        };
        let (region_counts, region_reads) = process_region(
            &args.bam_path,
            *tid,
            chr_name,
            *start,
            *end,
            &reference,
            &tid_to_name,
            &args,
            max_read_len,
            bed_for_filtering.as_deref(),
            bed_filter_whole_reads,
        );

        if !region_counts.is_empty() {
            let mut counts = global_counts.lock().expect("Lock poisoned");
            for (key, count) in region_counts {
                *counts.entry(key).or_insert(0) += count;
            }
        }

        let prev = total_reads.fetch_add(region_reads, Ordering::Relaxed);
        if (prev + region_reads) / 100_000 > prev / 100_000 {
            log::info!("Processed {} reads", prev + region_reads);
        }
    });

    log::info!(
        "Total reads processed: {}",
        total_reads.load(Ordering::Relaxed)
    );

    let mut counts = global_counts.lock().expect("Lock poisoned");

    if let Some(path) = args.discount_table.as_deref() {
        if args.position_mode != PositionMode::Read {
            log::warn!(
                "--discount-table was provided in {:?} mode. Discount rows are read-position keyed; applying anyway by matching base_position.",
                args.position_mode
            );
        }

        let discounts = load_discount_table(path)
            .unwrap_or_else(|e| panic!("Failed to load discount table from {}: {}", path, e));
        let touched = apply_external_discounts(&mut counts, discounts);
        log::info!("Applied external discounts to {} key groups", touched);
    }

    if args.emit_rescaling_matrix {
        if args.normalize {
            log::warn!("--normalize is ignored when --emit-rescaling-matrix is enabled.");
        }
        log::info!("Writing rescaling matrix rows...");
        write_rescaling_matrix_output(&counts, args.output_file.as_deref())
            .expect("Failed to write rescaling matrix output");
    } else if args.normalize {
        log::info!("Writing normalized frequencies...");
        write_normalized_output(&counts, args.output_file.as_deref(), args.position_mode)
            .expect("Failed to write normalized output");
    } else {
        write_output(&counts, args.output_file.as_deref(), args.position_mode)
            .expect("Failed to write output");
    }
}
