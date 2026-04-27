//! Input helpers for reference FASTA and BAM-derived metadata.

use crate::types::{
    DiscountKey, GenomicMismatchKey, InconsistencyKey, InsertKey, MismatchKey, PositionMode,
};
use crate::types::ReferenceGenome;
use bio::io::fasta;
use rust_htslib::bam::{Read, Reader};
use std::error::Error;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process;

/// Load a reference genome from a FASTA file.
///
/// # Arguments
/// * `fasta_path` - Path to the reference FASTA file.
///
/// # Returns
/// * A [`ReferenceGenome`] mapping chromosome names to base sequences.
///
/// # Panics
/// * If the FASTA file cannot be opened.
/// * If any FASTA record cannot be read.
pub fn load_reference_genome(fasta_path: &str) -> ReferenceGenome {
    eprintln!("Loading reference genome from: {}", fasta_path);
    let reader = fasta::Reader::from_file(fasta_path).expect("Failed to open reference FASTA file");

    let mut genome: ReferenceGenome = HashMap::new();
    for result in reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let chr_name = record.id().to_string();
        let sequence = record.seq().to_vec(); // Vec<u8> = byte, not UTF-8 char (overhead)
        genome.insert(chr_name, sequence);
    }

    eprintln!("Loaded {} chromosome(s)", genome.len());
    genome
}

/// Compute the mode (most common) read length from a sample of a BAM file.
///
/// # Arguments
/// * `bam_path` - Path to the input BAM file.
/// * `sample_size` - Number of records to sample from the start of the BAM stream.
///
/// # Returns
/// * The maximum read length (ideally the most common) in the sampled records.
/// * Returns a `default_len` if no valid records are observed in the sample.
///
/// # Panics
/// * If the BAM file cannot be opened.
pub fn compute_read_len_max_from_sample_bam(bam_path: &str, sample_size: usize) -> usize {
    let mut bam =
        Reader::from_path(bam_path).expect("Failed to open BAM file for read length sampling");

    let mut max_length: usize = 0;

    for result in bam.records().take(sample_size) {
        match result {
            Err(_) => {
                eprintln!(
                    "ERROR: Failed to sample BAM file for read length estimation. Exiting."
                );
                process::exit(1);
            }
            Ok(record) => {
                let read_len = record.seq().len();
                if read_len > max_length {
                    max_length = read_len;
                }
            }
        }
    }
    
    eprintln!("Computed max read length from sample: {}", max_length);
    max_length
}

/// Write genomic mismatch counts to a TSV file for potential variant review.
///
/// # Arguments
/// * `genomic_counts` - Genomic mismatch counts keyed by chromosome/position/mismatch.
/// * `output_path` - Output TSV path.
///
/// # Returns
/// * `Ok(())` on success.
/// * Any I/O error encountered while creating or writing the file.
pub fn write_potential_variants_tsv(
    genomic_counts: &HashMap<GenomicMismatchKey, (usize, usize)>,
    output_path: &str,
) -> std::io::Result<()> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth"
    )?;

    // Keep unsorted iteration for speed on large datasets.
    for (key, count) in genomic_counts.iter() {
        let parts: Vec<&str> = key.mismatch_type.split('>').collect();
        if parts.len() == 2 {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                key.chromosome, key.genomic_position, parts[0], parts[1], count.0, count.1
            )?;
        }
    }

    eprintln!(
        "Wrote {} genomic positions to {}",
        genomic_counts.len(),
        output_path
    );

    Ok(())
}

/// Print read pair inconsistency counts in CSV table format.
///
/// # Arguments
/// * `incons_types` - Sorted discordance type columns.
/// * `incons_positions` - Sorted `(read1_pos, read2_pos)` rows.
/// * `incons_position_map` - Mapping from position pairs to discordance counts.
pub fn print_read_pair_inconsistency_table(
    incons_types: &[String],
    incons_positions: &[(usize, usize)],
    incons_position_map: &HashMap<(usize, usize), HashMap<String, usize>>,
) {
    println!("\n# Read Pair Inconsistencies");
    print!("Read1_Pos,Read2_Pos");
    for incons_type in incons_types {
        print!(",{}", incons_type);
    }
    println!();

    for &(r1_pos, r2_pos) in incons_positions {
        print!("{},{}", r1_pos, r2_pos);
        if let Some(incons_map) = incons_position_map.get(&(r1_pos, r2_pos)) {
            for incons_type in incons_types {
                let count = incons_map.get(incons_type).unwrap_or(&0);
                print!(",{}", count);
            }
        }
        println!();
    }
}

/// Write read-pair overlap inconsistencies to a TSV file.
pub fn write_inconsistencies_tsv(
    inconsistency_counts: &HashMap<InconsistencyKey, usize>,
    output_path: &str,
) -> std::io::Result<()> {
    let mut rows: Vec<(&InconsistencyKey, &usize)> = inconsistency_counts.iter().collect();
    rows.sort_by(|(a, _), (b, _)| {
        a.read1_position
            .cmp(&b.read1_position)
            .then(a.read2_position.cmp(&b.read2_position))
            .then(a.discordance_type.cmp(&b.discordance_type))
    });

    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "read1_position\tread2_position\tdiscordance_type\tcount")?;

    for (key, count) in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            key.read1_position, key.read2_position, key.discordance_type, count
        )?;
    }

    Ok(())
}

/// Write mismatch-key discount counts to a TSV file.
pub fn write_mismatch_discounts_tsv(
    mismatch_discounts: &HashMap<MismatchKey, usize>,
    output_path: &str,
) -> std::io::Result<()> {
    let mut rows: Vec<(&MismatchKey, &usize)> = mismatch_discounts.iter().collect();
    rows.sort_by(|(a, _), (b, _)| {
        a.read_num
            .cmp(&b.read_num)
            .then(a.read_position.cmp(&b.read_position))
            .then(a.mismatch_type.cmp(&b.mismatch_type))
    });

    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "mismatch_type\tread_num\tread_position\tdiscount_count")?;

    for (key, count) in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            key.mismatch_type, key.read_num, key.read_position, count
        )?;
    }

    Ok(())
}

/// Load mismatch-rescaling matrix rows from a TSV file.
pub fn load_rescaling_matrix(
    path: &str,
) -> Result<HashMap<(u8, u16, char, char), f32>, Box<dyn Error>> {
    let mut matrix = HashMap::new();
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.trim().split('\t').collect();

        if parts.len() < 5 {
            continue;
        }

        let read_num: u8 = parts[0].parse()?;
        let position: u16 = parts[1].parse()?;
        let ref_base: char = parts[2].chars().next().ok_or("Invalid ref_base")?;
        let read_base: char = parts[3].chars().next().ok_or("Invalid read_base")?;
        let scaling_factor: f32 = parts[4].parse()?;

        matrix.insert((read_num, position, ref_base, read_base), scaling_factor);
    }

    Ok(matrix)
}

pub fn load_discount_table(path: &str) -> std::io::Result<HashMap<DiscountKey, usize>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut discounts: HashMap<DiscountKey, usize> = HashMap::new();

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        if line_idx == 0 {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 4 {
            continue;
        }

        let read_num = match fields[1].parse::<u8>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let base_position = match fields[2].parse::<usize>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let discount_count = match fields[3].parse::<usize>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let key = DiscountKey {
            base_change: fields[0].to_string(),
            read_num,
            base_position,
        };
        *discounts.entry(key).or_insert(0) += discount_count;
    }

    Ok(discounts)
}

pub fn apply_external_discounts(
    counts: &mut HashMap<InsertKey, usize>,
    discounts: HashMap<DiscountKey, usize>,
) -> usize {
    let mut touched = 0usize;

    for (discount_key, mut remaining) in discounts {
        let k1 = InsertKey {
            base_change: discount_key.base_change.clone(),
            read_num: discount_key.read_num,
            base_position: discount_key.base_position,
            reference_order: 1,
        };
        let k2 = InsertKey {
            base_change: discount_key.base_change,
            read_num: discount_key.read_num,
            base_position: discount_key.base_position,
            reference_order: 2,
        };

        let c1 = counts.get(&k1).copied().unwrap_or(0);
        let c2 = counts.get(&k2).copied().unwrap_or(0);

        if c1 == 0 && c2 == 0 {
            continue;
        }

        let first_key = if c1 >= c2 { &k1 } else { &k2 };
        if remaining > 0 {
            if let Some(v) = counts.get_mut(first_key) {
                let take = remaining.min(*v);
                *v -= take;
                remaining -= take;
            }
        }

        if remaining > 0 {
            let second_key = if first_key.reference_order == 1 {
                &k2
            } else {
                &k1
            };
            if let Some(v) = counts.get_mut(second_key) {
                let take = remaining.min(*v);
                *v -= take;
            }
        }

        touched += 1;
    }

    touched
}

pub fn write_output(
    counts: &HashMap<InsertKey, usize>,
    output_file: Option<&str>,
    position_mode: PositionMode,
) -> std::io::Result<()> {
    let position_label = match position_mode {
        PositionMode::Read => "read_position",
        PositionMode::Insert => "fragment_position",
    };

    let mut rows: Vec<(&InsertKey, &usize)> = counts.iter().collect();
    rows.sort_by(|(a, _), (b, _)| {
        a.reference_order
            .cmp(&b.reference_order)
            .then(a.read_num.cmp(&b.read_num))
            .then(a.base_position.cmp(&b.base_position))
            .then(a.base_change.cmp(&b.base_change))
    });

    if let Some(path) = output_file {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "base_change\tread_num\treference_order\t{}\tcount",
            position_label
        )?;
        for (key, count) in rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                key.base_change, key.read_num, key.reference_order, key.base_position, count
            )?;
        }
        return Ok(());
    }

    println!(
        "base_change\tread_num\treference_order\t{}\tcount",
        position_label
    );
    for (key, count) in rows {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            key.base_change, key.read_num, key.reference_order, key.base_position, count
        );
    }
    Ok(())
}
