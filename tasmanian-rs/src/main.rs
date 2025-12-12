//! Tasmanian CLI entry point
//!
//! Usage:
//!   samtools view data.bam | tasmanian intersections -b regions.bed | tasmanian analyze -r ref.fa > output.csv
//!   samtools view data.bam | tasmanian analyze -r ref.fa > output.csv

use anyhow::Result;
use clap::Parser;
use log::info;
use std::io::{self, BufReader, BufWriter, Write};

use tasmanian::cli::{Cli, Commands};
use tasmanian::{analyze, write_csv, AnalysisConfig, BedIndex, IntersectionProcessor, Reference};

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Intersections(args) => {
            run_intersections(args)?;
        }
        Commands::Analyze(args) => {
            run_analyze(args)?;
        }
    }

    Ok(())
}

fn run_intersections(args: tasmanian::cli::IntersectionsArgs) -> Result<()> {
    if args.debug {
        log::set_max_level(log::LevelFilter::Debug);
    }

    eprintln!("Loading BED file: {}", args.bed_file.display());
    let bed_index = BedIndex::load(&args.bed_file)?;
    eprintln!(
        "Loaded {} regions across {} chromosomes",
        bed_index.total_regions(),
        bed_index.num_chromosomes()
    );

    let stdin = io::stdin().lock();
    let reader = BufReader::new(stdin);

    let stdout = io::stdout().lock();
    let writer = BufWriter::new(stdout);

    let mut processor = IntersectionProcessor::new(bed_index);
    processor.process_stream(reader, writer)?;

    eprintln!("Intersection processing complete");
    Ok(())
}

fn run_analyze(args: tasmanian::cli::AnalyzeArgs) -> Result<()> {
    use noodles::sam;

    if args.debug {
        log::set_max_level(log::LevelFilter::Debug);
    }

    eprintln!("Loading reference: {}", args.reference_fasta.display());
    let reference = Reference::load(&args.reference_fasta)?;
    eprintln!(
        "Loaded {} sequences from reference",
        reference.num_sequences()
    );

    let config = AnalysisConfig {
        min_base_quality: args.base_quality,
        min_mapping_quality: args.mapping_quality,
        filter_indels: args.filter_indel,
        min_read_length: args.filter_length.0,
        max_read_length: args.filter_length.1,
        min_fragment_length: args.fragment_length.0,
        max_fragment_length: args.fragment_length.1,
        softclip_bypass: args.soft_clip_bypass,
        unmask_genome: args.unmask_genome,
        confidence_threshold: args.confidence,
        ont_mode: args.ont,
        mask_methyl_c: args.mask_methyl_c,
        mask_methyl_cpg: args.mask_methyl_cpg,
    };

    let stdin = io::stdin().lock();
    let mut reader = sam::io::Reader::new(BufReader::new(stdin));
    let header = reader.read_header()?;

    let (tables, read_length) = analyze(&mut reader, &header, &reference, &config)?;

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);
    write_csv(&mut writer, &tables, read_length)?;
    writer.flush()?;

    // Write HTML report if output prefix specified
    if let Some(prefix) = &args.output_prefix {
        let report_path = format!("{}.html", prefix);
        eprintln!("HTML report generation not yet implemented: {}", report_path);
    }

    Ok(())
}
