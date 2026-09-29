use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use tasmanian_mismatch::{
    Args, BedFilter, BedFilterMode, BlockCounts, InsertKey, PositionMode, ProcessingConfig,
    ProcessingContext, WindowSummary, apply_external_discounts, block_residual_diagnostics,
    bootstrap_intervals, build_tid_map_and_regions, compute_read_len_max_from_sample_bam,
    configure_thread_pool, launch_visualization, load_discount_table, load_reference_genome,
    mask_reference_with_bed, maybe_parse_bed_file, position_label, process_region,
    process_region_windows, summarize_window_counts, write_block_diagnostics,
    write_bootstrap_output, write_normalized_output, write_output, write_rate_intervals,
    write_rescaling_matrix_output, write_window_output,
};

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

    configure_thread_pool(args.threads);

    let (max_read_len, mut reference, bed_regions) = std::thread::scope(|s| {
        let t1 = s.spawn(|| compute_read_len_max_from_sample_bam(&args.bam_path, 10_000));
        let t2 = s.spawn(|| load_reference_genome(&args.reference_path));
        let t3 = s.spawn(|| maybe_parse_bed_file(args.bed_file.as_deref()));
        (t1.join().unwrap(), t2.join().unwrap(), t3.join().unwrap())
    });
    log::info!("Sampled max read length: {}", max_read_len);

    let bed_for_filtering = if let Some(regions) = bed_regions {
        log::info!("BED filter mode: {:?}", args.bed_filter_mode);
        match args.bed_filter_mode {
            BedFilterMode::Filter => {
                log::info!("Keeping BED regions for whole-read overlap filtering...");
                Some(Arc::new(regions))
            }
            BedFilterMode::Include => {
                log::info!(
                    "Keeping only reads overlapping BED regions (in-silico exome/panel restriction)..."
                );
                Some(Arc::new(regions))
            }
            BedFilterMode::Mask => {
                log::info!("Masking reference genome at BED regions...");
                let masked_bases = mask_reference_with_bed(&mut reference, &regions);
                log::info!("Masked {} bases in reference genome", masked_bases);
                None
            }
        }
    } else {
        None
    };

    let reference = Arc::new(reference);
    // In window mode, round chunks up to whole windows so no window straddles two chunks.
    let chunk_size = args.window_size.map_or(args.region_size, |size| {
        (args.region_size as u64).div_ceil(size).max(1) as usize * size as usize
    });
    let (tid_to_name, regions) = build_tid_map_and_regions(&args.bam_path, chunk_size);
    log::info!("Processing {} indexed regions", regions.len());

    let config = ProcessingConfig {
        softclip_threshold: 0.66,
        min_base_quality: args.min_base_quality,
        is_methylation: args.methylation_mode,
        mode_len: max_read_len,
        min_map_quality: args.min_map_quality,
        required_flags: args.required_flags,
        filter_flags: args.filter_flags,
        excl_flags: args.excl_flags,
        use_insert_mode: args.position_mode == PositionMode::Insert,
        position_mode: args.position_mode,
        overlap_mode: args.overlap_mode,
        min_fragment_length: args.min_fragment_length,
        max_fragment_length: args.max_fragment_length,
        min_read_position: args.min_read_position.unwrap_or(0),
        max_read_position: args.max_read_position.unwrap_or(usize::MAX),
    };

    let context = ProcessingContext {
        reference: &reference,
        tid_to_name: &tid_to_name,
        bed_intervals: &[],
    };
    let bed_filter = BedFilter {
        regions: bed_for_filtering.as_deref(),
        filter_whole_reads: args.bed_filter_mode.filters_whole_reads(),
        include_only: args.bed_filter_mode.include_only(),
    };
    let total_reads = AtomicUsize::new(0);
    let record_reads = |region_reads: usize| {
        let prev = total_reads.fetch_add(region_reads, Ordering::Relaxed);
        if (prev + region_reads) / 100_000 > prev / 100_000 {
            log::info!("Processed {} reads", prev + region_reads);
        }
    };

    if let Some(window_size) = args.window_size {
        let window_size = window_size as i64;
        // Regions are in genomic order and rayon's collect keeps it, so windows come out sorted.
        let windows: Vec<WindowSummary> = regions
            .par_iter()
            .flat_map_iter(|region| {
                let mut windows = Vec::new();
                let region_reads = process_region_windows(
                    &args.bam_path,
                    region,
                    &context,
                    config,
                    &bed_filter,
                    window_size,
                    |start, counts| {
                        let rows = summarize_window_counts(counts);
                        if !rows.is_empty() {
                            windows.push(WindowSummary {
                                tid: region.tid,
                                start,
                                end: (start + window_size).min(region.end),
                                rows,
                            });
                        }
                    },
                );
                record_reads(region_reads);
                windows
            })
            .collect();
        log::info!(
            "Total reads processed: {}",
            total_reads.load(Ordering::Relaxed)
        );
        log::info!("Writing {} window summaries...", windows.len());
        write_window_output(&windows, &tid_to_name, args.output_file.as_deref())
            .expect("Failed to write window output");
        return;
    }

    let global_counts: Mutex<HashMap<InsertKey, usize>> = Mutex::new(HashMap::new());
    let bootstrap_blocks: Option<Mutex<BlockCounts>> =
        args.bootstrap.map(|_| Mutex::new(BlockCounts::new()));

    regions.par_iter().for_each(|region| {
        let (region_counts, region_reads) =
            process_region(&args.bam_path, region, &context, config, &bed_filter);

        if !region_counts.is_empty() {
            if let Some(blocks) = &bootstrap_blocks {
                blocks
                    .lock()
                    .expect("Lock poisoned")
                    .add_block((region.tid, region.start), &region_counts);
            }
            let mut counts = global_counts.lock().expect("Lock poisoned");
            for (key, count) in region_counts {
                *counts.entry(key).or_insert(0) += count;
            }
        }
        record_reads(region_reads);
    });

    log::info!(
        "Total reads processed: {}",
        total_reads.load(Ordering::Relaxed)
    );

    let mut counts = global_counts.into_inner().expect("Lock poisoned");
    let mut discounted_amounts: HashMap<InsertKey, usize> = HashMap::new();

    if let Some(path) = args.discount_table.as_deref() {
        if args.position_mode != PositionMode::Read {
            log::warn!(
                "--discount-table was provided in {:?} mode. Discount rows are read-position keyed; applying anyway by matching base_position.",
                args.position_mode
            );
        }

        let discounts = load_discount_table(path)
            .unwrap_or_else(|e| panic!("Failed to load discount table from {}: {}", path, e));
        discounted_amounts = apply_external_discounts(&mut counts, discounts);
        log::info!(
            "Applied external discounts: removed {} observations from {} keys",
            discounted_amounts.values().sum::<usize>(),
            discounted_amounts.len()
        );
    }

    let intervals = bootstrap_blocks.map(|blocks| {
        run_bootstrap(
            &blocks.into_inner().expect("Lock poisoned"),
            &args,
            &tid_to_name,
            &discounted_amounts,
        )
    });

    if args.emit_rescaling_matrix {
        if args.normalize {
            log::warn!("--normalize is ignored when --emit-rescaling-matrix is enabled.");
        }
        log::info!("Writing rescaling matrix rows...");
        write_rescaling_matrix_output(&counts, args.output_file.as_deref())
            .expect("Failed to write rescaling matrix output");
    } else if let Some(intervals) = &intervals {
        write_bootstrap_output(
            &counts,
            intervals,
            args.output_file.as_deref(),
            args.position_mode,
            args.normalize,
        )
        .expect("Failed to write bootstrap output");
    } else if args.normalize {
        log::info!("Writing normalized frequencies...");
        write_normalized_output(&counts, args.output_file.as_deref(), args.position_mode)
            .expect("Failed to write normalized output");
    } else {
        write_output(&counts, args.output_file.as_deref(), args.position_mode)
            .expect("Failed to write output");
    }

    if args.plot {
        if args.emit_rescaling_matrix {
            log::warn!("--plot is not supported with --emit-rescaling-matrix; skipping.");
        } else {
            // If -o was given, use that file; otherwise write a temp TSV for the visualizer
            let tsv_path = if let Some(ref p) = args.output_file {
                p.clone()
            } else {
                let tmp = std::env::temp_dir().join("tasmanian_output.tsv");
                let tmp_str = tmp.to_str().expect("temp path not valid UTF-8").to_string();
                if args.normalize {
                    write_normalized_output(&counts, Some(&tmp_str), args.position_mode)
                        .expect("Failed to write temp TSV for visualization");
                } else {
                    write_output(&counts, Some(&tmp_str), args.position_mode)
                        .expect("Failed to write temp TSV for visualization");
                }
                tmp_str
            };
            launch_visualization(&tsv_path).expect("Visualization failed");
        }
    }
}

/// Report block diagnostics and per-position rate intervals as requested, and return the
/// frequency interval of each key for the main table.
fn run_bootstrap(
    blocks: &BlockCounts,
    args: &Args,
    tid_to_name: &HashMap<i32, String>,
    discounted_amounts: &HashMap<InsertKey, usize>,
) -> HashMap<InsertKey, (f64, f64)> {
    let replicates = args.bootstrap.expect("bootstrap blocks imply --bootstrap") as usize;
    if blocks.len() < 20 {
        log::warn!(
            "Only {} non-empty blocks to resample; bootstrap intervals will be unreliable. \
             Lower --region-size for more, smaller blocks.",
            blocks.len()
        );
    }
    log::info!(
        "Bootstrapping {} replicates over {} blocks of up to {} bp...",
        replicates,
        blocks.len(),
        args.region_size
    );

    if let Some(path) = args.bootstrap_diagnostics.as_deref() {
        let diagnostics = block_residual_diagnostics(blocks);
        for diag in &diagnostics {
            let extreme = diag
                .most_extreme
                .map_or("none".to_string(), |((tid, start), z)| {
                    format!("{}:{} (z = {:+.1})", tid_to_name[&tid], start, z)
                });
            log::info!(
                "Block dispersion read {} {}: {:.2} over {} blocks ({} too sparse); most extreme {}",
                diag.read_num,
                diag.base_change,
                diag.dispersion,
                diag.blocks_used,
                diag.blocks_sparse,
                extreme
            );
        }
        write_block_diagnostics(&diagnostics, path).expect("Failed to write bootstrap diagnostics");
    }

    let intervals = bootstrap_intervals(
        blocks,
        replicates,
        args.bootstrap_seed,
        discounted_amounts,
        args.bootstrap_rates.is_some(),
    );
    if let (Some(path), Some(rates)) = (args.bootstrap_rates.as_deref(), &intervals.rates) {
        write_rate_intervals(rates, position_label(args.position_mode), path)
            .expect("Failed to write bootstrap rates");
    }
    intervals.frequencies
}
