use clap::Parser;
use rayon::prelude::*;
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Reader, Record};
use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use tasmanian_mismatch::{OverlapMode, PositionMode, *};

#[derive(Parser, Debug)]
#[command(name = "tasmanian-diagnostics")]
#[command(
    version,
    about = "Collect genomic variants and read-pair overlap inconsistencies"
)]
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

    /// Read length for the position axis; use the same value as tasmanian-mismatch
    /// --max-read-length. Omitted or given without a value, the longest of the first
    /// 10,000 reads is used, as in tasmanian-mismatch.
    #[arg(long, num_args = 0..=1, value_name = "LENGTH")]
    use_read_len_max: Option<Option<u32>>,

    /// Keep processing behavior compatible with split-read vs insert-position mode
    #[arg(long, default_value_t = false)]
    use_insert_mode: bool,

    /// Minimum estimated fragment length required for a read to be counted
    #[arg(long, default_value_t = 25)]
    min_fragment_length: usize,

    /// Maximum estimated fragment length allowed for a read to be counted
    #[arg(long, default_value_t = 10_000)]
    max_fragment_length: usize,

    /// Minimum number of supporting reads for reporting a genomic site
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

    /// How the BED file is applied
    #[arg(long, value_enum, default_value = "mask", requires = "bed_file")]
    bed_filter_mode: BedFilterMode,

    /// Output path for genomic potential variants table
    #[arg(long, default_value = "potential_variants.tsv")]
    variants_output: String,

    /// Output path for overlap inconsistencies table
    #[arg(long, default_value = "read_pair_inconsistencies.tsv")]
    inconsistencies_output: String,

    /// Output path for overlap opportunities: compared and discordant R1/R2 base pairs per
    /// (read1_position, read2_position, read1_base, min_quality_bin). The denominator for
    /// the inconsistencies table.
    #[arg(long, default_value = "read_pair_overlap_opportunities.tsv")]
    opportunities_output: String,

    /// Also compare each optical duplicate (DT:Z:SQ from Picard, dt:Z:SQ from
    /// samtools markdup -d) with its original, cycle by cycle, and write the counts here.
    /// Differences between them come from the sequencer, not the library. Off unless given.
    #[arg(long, value_name = "TSV")]
    optical_duplicates_output: Option<String>,

    /// Output path for mismatch-key discounts to be consumed by tasmanian-mismatch.
    /// Use '-' to write the discount table to stdout.
    #[arg(long, default_value = "variant_discounts.tsv")]
    discount_output: String,
}

/// Mates whose start positions are further apart than this are not paired. A read can only
/// overlap a mate that starts inside its own reference span, so only reads spanning more
/// than this (very large deletions) lose an overlap; the window also bounds how far past
/// its end each region reads to find mates of the pairs it owns.
const PAIR_WINDOW: i64 = 1_000;
/// Reference span allowed for a read beyond PAIR_WINDOW when deciding which positions near
/// a region start may also receive bases from the previous region's reads.
const SPAN_MARGIN: i64 = 1_000;

#[derive(Clone, Copy)]
struct SiteThresholds {
    count: usize,
    depth: usize,
}

/// Raw counts for genomic sites near region edges, merged across regions.
#[derive(Default)]
struct BoundaryCounts {
    genomic: HashMap<GenomicMismatchKey, GenomicMismatchValue>,
    depth: HashMap<(i32, i64), usize>,
}

/// Per-region accumulators.
#[derive(Default)]
struct RegionCounts {
    mismatch: HashMap<MismatchKey, usize>,
    overlap_mismatch: HashMap<MismatchKey, usize>,
    inconsistency: HashMap<InconsistencyKey, usize>,
    opportunity: HashMap<OverlapOpportunityKey, OverlapOpportunityCounts>,
    genomic: HashMap<GenomicMismatchKey, GenomicMismatchValue>,
    depth: HashMap<i64, usize>,
    records: usize,
    overlap_pairs: usize,
}

impl RegionCounts {
    fn single(&mut self, record: &Record, context: &ProcessingContext, config: ProcessingConfig) {
        process_single_record(
            record,
            &mut self.mismatch,
            &mut self.genomic,
            &mut self.depth,
            context,
            config,
        );
        self.records += 1;
    }

    fn pair(
        &mut self,
        record1: &Record,
        record2: &Record,
        context: &ProcessingContext,
        config: ProcessingConfig,
    ) {
        let info1 = ReadInfo {
            tid: record1.tid(),
            pos: record1.pos(),
        };
        let info2 = ReadInfo {
            tid: record2.tid(),
            pos: record2.pos(),
        };
        let Some(overlap) = get_overlap_region(&info1, &info2, record1, record2) else {
            self.single(record1, context, config);
            self.single(record2, context, config);
            return;
        };
        self.overlap_pairs += 1;
        let mut overlap_counts = OverlapCounts {
            local_counts: &mut self.mismatch,
            overlap_counts: &mut self.overlap_mismatch,
            inconsistency_counts: &mut self.inconsistency,
            opportunity_counts: Some(&mut self.opportunity),
            genomic_counts: Some(&mut self.genomic),
            genomic_depth: Some(&mut self.depth),
        };
        process_paired_reads_with_overlap(
            record1,
            record2,
            overlap,
            &mut overlap_counts,
            context,
            config,
        );
        self.records += 2;
    }
}

fn has_mapped_pair(record: &Record) -> bool {
    record.is_paired() && !record.is_unmapped() && !record.is_mate_unmapped()
}

/// Primary alignments whose mate is close enough on the same contig to overlap them.
fn is_pairable(record: &Record) -> bool {
    has_mapped_pair(record)
        && !record.is_secondary()
        && !record.is_supplementary()
        && record.tid() == record.mtid()
        && (record.pos() - record.mpos()).abs() < PAIR_WINDOW
}

/// Keep a genomic site whose supporting reads and depth pass the thresholds, and count its
/// distinct mismatch keys as discounts.
fn report_site(
    key: GenomicMismatchKey,
    site: GenomicMismatchValue,
    depth: usize,
    thresholds: SiteThresholds,
    genomic_counts: &mut HashMap<GenomicMismatchKey, SiteCounts>,
    discounts: &mut HashMap<MismatchKey, usize>,
) {
    if site.supporting_reads < thresholds.count || depth < thresholds.depth {
        return;
    }
    genomic_counts.insert(
        key,
        SiteCounts {
            supporting_reads: site.supporting_reads,
            depth,
            distinct_read_positions: site.distinct_read_positions(),
        },
    );
    for mismatch_key in site.mismatch_keys {
        *discounts.entry(mismatch_key).or_insert(0) += 1;
    }
}

fn main() {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("Failed to set thread pool size");
    }

    // Sampled like tasmanian-mismatch so both tools use the same position axis.
    let mode_len = match args.use_read_len_max.flatten() {
        Some(val) => val as usize,
        None => compute_read_len_max_from_sample_bam(&args.bam_file, 10_000),
    };
    log::info!("Max read length: {}", mode_len);
    if args.use_insert_mode {
        log::warn!(
            "--use-insert-mode: discount table positions use diagnostics' own insert axis and do not \
             match tasmanian-mismatch --position-mode insert; use read mode (no --use-insert-mode) \
             for a table to pass to --discount-table"
        );
    }

    let mut reference = load_reference_genome(&args.reference_fasta);

    let bed_for_filtering = if let Some(regions) = maybe_parse_bed_file(args.bed_file.as_deref()) {
        match args.bed_filter_mode {
            BedFilterMode::Filter | BedFilterMode::Include => Some(regions),
            BedFilterMode::Mask => {
                let masked_bases = mask_reference_with_bed(&mut reference, &regions);
                log::info!("Masked {} bases in reference genome", masked_bases);
                None
            }
        }
    } else {
        None
    };

    let bam = Reader::from_path(&args.bam_file).expect("Failed to open BAM file");

    let tid_to_name: HashMap<i32, String> = (0..bam.header().target_count())
        .map(|i| {
            (
                i as i32,
                String::from_utf8_lossy(bam.header().tid2name(i)).to_string(),
            )
        })
        .collect();

    let mut regions = Vec::new();
    for (&tid, chr_name) in tid_to_name.iter() {
        let chr_len = bam.header().target_len(tid as u32).unwrap_or(0) as i64;
        let mut start = 0i64;
        while start < chr_len {
            let end = std::cmp::min(start + args.region_size as i64, chr_len);
            regions.push((tid, start, end, chr_len, chr_name.clone()));
            start = end;
        }
    }
    drop(bam);

    log::info!("Created {} regions to process", regions.len());

    let total_count = AtomicUsize::new(0);
    let overlap_pairs_count = AtomicUsize::new(0);
    let overreach_positions = AtomicUsize::new(0);
    let genomic_mismatch_counts: Mutex<HashMap<GenomicMismatchKey, SiteCounts>> =
        Mutex::new(HashMap::new());
    let inconsistency_counts: Mutex<HashMap<InconsistencyKey, usize>> = Mutex::new(HashMap::new());
    let mismatch_discounts: Mutex<HashMap<MismatchKey, usize>> = Mutex::new(HashMap::new());
    let opportunity_counts: Mutex<HashMap<OverlapOpportunityKey, OverlapOpportunityCounts>> =
        Mutex::new(HashMap::new());
    let boundary_counts: Mutex<BoundaryCounts> = Mutex::new(BoundaryCounts::default());
    let optical_counts: Mutex<HashMap<OpticalDuplicateKey, usize>> = Mutex::new(HashMap::new());
    let optical_stats: Mutex<OpticalDuplicateStats> = Mutex::new(OpticalDuplicateStats::default());
    let collect_optical = args.optical_duplicates_output.is_some();
    let thresholds = SiteThresholds {
        count: args.genomic_threshold,
        depth: args.genomic_depth_threshold,
    };

    let processing_config = ProcessingConfig {
        softclip_threshold: args.softclip_threshold,
        min_base_quality: args.min_base_quality,
        is_methylation: args.methylation,
        mode_len,
        min_map_quality: args.min_map_quality,
        required_flags: args.required_flags,
        filter_flags: args.filter_flags,
        excl_flags: args.excl_flags,
        use_insert_mode: args.use_insert_mode,
        position_mode: if args.use_insert_mode {
            PositionMode::Insert
        } else {
            PositionMode::Read
        },
        overlap_mode: OverlapMode::Cut,
        min_fragment_length: args.min_fragment_length,
        max_fragment_length: args.max_fragment_length,
        min_read_position: 0,
        max_read_position: usize::MAX,
    };

    regions
        .par_iter()
        .for_each(|(tid, start, end, chr_len, chr_name)| {
            let (tid, start, end) = (*tid, *start, *end);
            // Read past the region end so mates of pairs owned here (leftmost mate starts in
            // [start, end)) are found; every such mate starts before end + PAIR_WINDOW.
            let fetch_end = std::cmp::min(end + PAIR_WINDOW, *chr_len);

            let mut bam =
                IndexedReader::from_path(args.bam_file.as_str()).expect("Failed to open BAM file");
            if let Err(e) = bam.fetch(FetchDefinition::Region(tid, start, fetch_end)) {
                log::warn!(
                    "Failed to fetch region {}:{}-{}: {}",
                    chr_name,
                    start,
                    end,
                    e
                );
                return;
            }

            let (core_bed, extension_bed) = match bed_for_filtering.as_ref() {
                Some(bed) => (
                    filter_bed_for_region(bed, chr_name, start, end),
                    filter_bed_for_region(bed, chr_name, end, fetch_end),
                ),
                None => (Vec::new(), Vec::new()),
            };
            let filters_whole_reads = args.bed_filter_mode.filters_whole_reads();
            let include_only = args.bed_filter_mode.include_only();
            let mut core_cursor = 0usize;
            let mut extension_cursor = 0usize;

            let processing_context = ProcessingContext {
                reference: &reference,
                tid_to_name: &tid_to_name,
                // BED regions act per read (skip check below), never per base: mask mode masks
                // the reference instead, and per-base masking would erase every read `include` keeps.
                bed_intervals: &[],
            };
            let cfg = processing_config;

            let mut counts = RegionCounts::default();
            // First-seen mate of each pair owned by this region, until its mate arrives.
            let mut pending: HashMap<Vec<u8>, Record> = HashMap::new();
            // Optical-duplicate comparison groups whole regions, so it alone still buffers reads.
            let mut optical_records: Vec<Record> = Vec::new();

            for result in bam.records() {
                let Ok(record) = result else {
                    continue;
                };
                let pos = record.pos();
                if pos < start {
                    continue;
                }

                if pos < end {
                    if should_skip_whole_read_for_bed(
                        &record,
                        filters_whole_reads,
                        include_only,
                        &core_bed,
                        &mut core_cursor,
                    ) {
                        continue;
                    }
                    if collect_optical && has_mapped_pair(&record) {
                        optical_records.push(record.clone());
                    }
                    if !is_pairable(&record) {
                        counts.single(&record, &processing_context, cfg);
                        continue;
                    }
                    let mate_pos = record.mpos();
                    if mate_pos < start {
                        // The region holding the mate owns this pair and reads it from there.
                        continue;
                    }
                    if let Some(mate) = pending.remove(record.qname()) {
                        counts.pair(&mate, &record, &processing_context, cfg);
                    } else if mate_pos >= pos {
                        pending.insert(record.qname().to_vec(), record);
                    } else {
                        // Mate should already have been read but was filtered or is absent.
                        counts.single(&record, &processing_context, cfg);
                    }
                } else {
                    // Extension: only mates of pairs owned by this region. Their own region
                    // skipped them, so each is processed here exactly once.
                    if !is_pairable(&record) || record.mpos() < start || record.mpos() >= end {
                        continue;
                    }
                    if should_skip_whole_read_for_bed(
                        &record,
                        filters_whole_reads,
                        include_only,
                        &extension_bed,
                        &mut extension_cursor,
                    ) {
                        continue;
                    }
                    if let Some(mate) = pending.remove(record.qname()) {
                        counts.pair(&mate, &record, &processing_context, cfg);
                    } else {
                        counts.single(&record, &processing_context, cfg);
                    }
                }
            }
            for (_qname, record) in pending.drain() {
                counts.single(&record, &processing_context, cfg);
            }

            // Duplicates and their originals start at the same 5' end, so both are in this region.
            if collect_optical {
                let mut region_counts: HashMap<OpticalDuplicateKey, usize> = HashMap::new();
                let mut region_stats = OpticalDuplicateStats::default();
                collect_optical_duplicate_comparisons(
                    optical_records.iter(),
                    &processing_config,
                    &mut region_counts,
                    &mut region_stats,
                );
                *optical_stats.lock().expect("Lock poisoned") += region_stats;
                if !region_counts.is_empty() {
                    let mut global = optical_counts.lock().expect("Lock poisoned");
                    for (key, count) in region_counts {
                        *global.entry(key).or_insert(0) += count;
                    }
                }
            }

            // Genomic sites near either edge also get reads from the neighbouring region, so
            // they are merged across regions and thresholded only once all regions are done.
            // Reads from the previous region reach at most PAIR_WINDOW + a read span past its end.
            let core_start = start + PAIR_WINDOW + SPAN_MARGIN;
            let in_core = |p: i64| p >= core_start && p < end;
            let overreach_limit = end + PAIR_WINDOW + SPAN_MARGIN;
            let RegionCounts {
                inconsistency: inconsistency_region_counts,
                opportunity: opportunity_region_counts,
                genomic: genomic_region_counts,
                depth: genomic_position_depth,
                records: local_count,
                overlap_pairs: local_overlap_count,
                ..
            } = counts;

            let mut boundary_sites = Vec::new();
            if !genomic_region_counts.is_empty() {
                let mut global_genomic_counts =
                    genomic_mismatch_counts.lock().expect("Lock poisoned");
                let mut global_discount_counts = mismatch_discounts.lock().expect("Lock poisoned");
                for (genomic_key, genomic_value) in genomic_region_counts {
                    let position = genomic_key.genomic_position;
                    if !in_core(position) {
                        boundary_sites.push((genomic_key, genomic_value));
                        continue;
                    }
                    let depth = genomic_position_depth.get(&position).copied().unwrap_or(0);
                    report_site(
                        genomic_key,
                        genomic_value,
                        depth,
                        thresholds,
                        &mut global_genomic_counts,
                        &mut global_discount_counts,
                    );
                }
            }
            {
                let mut boundary = boundary_counts.lock().expect("Lock poisoned");
                for (key, site) in boundary_sites {
                    boundary.genomic.entry(key).or_default().merge(site);
                }
                let mut overreach = 0usize;
                for (&position, &depth) in genomic_position_depth.iter() {
                    if !in_core(position) {
                        *boundary.depth.entry((tid, position)).or_insert(0) += depth;
                        if position >= overreach_limit {
                            overreach += 1;
                        }
                    }
                }
                if overreach > 0 {
                    overreach_positions.fetch_add(overreach, Ordering::Relaxed);
                }
            }

            if !inconsistency_region_counts.is_empty() {
                let mut global = inconsistency_counts.lock().expect("Lock poisoned");
                for (key, count) in inconsistency_region_counts {
                    *global.entry(key).or_insert(0) += count;
                }
            }

            if !opportunity_region_counts.is_empty() {
                let mut global_opportunities = opportunity_counts.lock().expect("Lock poisoned");
                for (key, value) in opportunity_region_counts {
                    let entry = global_opportunities.entry(key).or_default();
                    entry.compared += value.compared;
                    entry.discordant += value.discordant;
                }
            }

            if local_overlap_count > 0 {
                overlap_pairs_count.fetch_add(local_overlap_count, Ordering::Relaxed);
            }

            let prev_total = total_count.fetch_add(local_count, Ordering::Relaxed);
            if (prev_total + local_count) / 100_000 > prev_total / 100_000 {
                log::info!("Processed {} records...", prev_total + local_count);
            }
        });

    // Sites near region edges, now complete across regions.
    {
        let boundary = boundary_counts.into_inner().expect("Lock poisoned");
        let name_to_tid: HashMap<&str, i32> = tid_to_name
            .iter()
            .map(|(&tid, name)| (name.as_str(), tid))
            .collect();
        let mut global_genomic_counts = genomic_mismatch_counts.lock().expect("Lock poisoned");
        let mut global_discount_counts = mismatch_discounts.lock().expect("Lock poisoned");
        for (genomic_key, site) in boundary.genomic {
            let tid = name_to_tid[genomic_key.chromosome.as_str()];
            let depth = boundary
                .depth
                .get(&(tid, genomic_key.genomic_position))
                .copied()
                .unwrap_or(0);
            report_site(
                genomic_key,
                site,
                depth,
                thresholds,
                &mut global_genomic_counts,
                &mut global_discount_counts,
            );
        }
    }
    let overreach = overreach_positions.load(Ordering::Relaxed);
    if overreach > 0 {
        log::warn!(
            "{} genomic positions were covered by reads reaching more than {} bp past their \
             region end; their depth near region edges may be split between regions",
            overreach,
            PAIR_WINDOW + SPAN_MARGIN
        );
    }

    log::info!(
        "Total records processed: {}",
        total_count.load(Ordering::Relaxed)
    );
    log::info!(
        "Read pairs with overlaps: {}",
        overlap_pairs_count.load(Ordering::Relaxed)
    );

    {
        let genomic_counts = genomic_mismatch_counts.lock().expect("Lock poisoned");
        write_potential_variants_tsv(&genomic_counts, &args.variants_output)
            .expect("Failed to write variants output");
        log::info!("Wrote variants table to {}", args.variants_output);
    }

    {
        let incons_counts = inconsistency_counts.lock().expect("Lock poisoned");
        write_inconsistencies_tsv(&incons_counts, &args.inconsistencies_output)
            .expect("Failed to write inconsistencies output");
        log::info!(
            "Wrote {} inconsistency keys to {}",
            incons_counts.len(),
            args.inconsistencies_output
        );
    }

    {
        let opportunities = opportunity_counts.lock().expect("Lock poisoned");
        let (compared, discordant) = opportunities.values().fold((0usize, 0usize), |(c, d), v| {
            (c + v.compared, d + v.discordant)
        });
        log::info!(
            "Overlap base pairs compared (both Q>={}): {}, discordant: {} ({:.3e})",
            args.min_base_quality,
            compared,
            discordant,
            if compared > 0 {
                discordant as f64 / compared as f64
            } else {
                0.0
            }
        );
        write_overlap_opportunities_tsv(&opportunities, &args.opportunities_output)
            .expect("Failed to write overlap opportunities output");
        log::info!(
            "Wrote {} overlap opportunity keys to {}",
            opportunities.len(),
            args.opportunities_output
        );
    }

    if let Some(path) = &args.optical_duplicates_output {
        let stats = *optical_stats.lock().expect("Lock poisoned");
        let counts = optical_counts.lock().expect("Lock poisoned");
        let (compared, discordant) = counts.iter().fold((0usize, 0usize), |(c, d), (k, n)| {
            (
                c + n,
                d + if k.original_base != k.duplicate_base {
                    *n
                } else {
                    0
                },
            )
        });
        log::info!(
            "Optical duplicates: {} SQ-tagged reads; {} sets compared ({} duplicates); skipped {} \
             sets with library/untagged duplicates, {} without an original, {} with several",
            stats.optical_reads,
            stats.sets_compared,
            stats.duplicates_compared,
            stats.sets_skipped_mixed,
            stats.sets_skipped_no_original,
            stats.sets_skipped_multiple_originals
        );
        log::info!(
            "Optical duplicate bases compared (both Q>={}): {}, discordant: {} ({:.3e})",
            args.min_base_quality,
            compared,
            discordant,
            if compared > 0 {
                discordant as f64 / compared as f64
            } else {
                0.0
            }
        );
        if stats.optical_reads == 0 {
            log::warn!(
                "--optical-duplicates-output: no reads tagged as optical duplicates. Mark them with \
                 `samtools markdup -d <distance>` or Picard MarkDuplicates TAGGING_POLICY=All \
                 (with an optical pixel distance suited to the flowcell)"
            );
        }
        write_optical_duplicates_tsv(&counts, path)
            .expect("Failed to write optical duplicates output");
        log::info!("Wrote {} optical duplicate keys to {}", counts.len(), path);
    }

    {
        let discounts = mismatch_discounts.lock().expect("Lock poisoned");
        if args.discount_output == "-" {
            let stdout = std::io::stdout();
            let mut handle = stdout.lock();
            write_mismatch_discounts_to_writer(&discounts, &mut handle)
                .expect("Failed to write mismatch discount output to stdout");
            log::info!("Wrote {} mismatch discount keys to stdout", discounts.len());
        } else {
            write_mismatch_discounts_tsv(&discounts, &args.discount_output)
                .expect("Failed to write mismatch discount output");
            log::info!(
                "Wrote {} mismatch discount keys to {}",
                discounts.len(),
                args.discount_output
            );
        }
    }
}
