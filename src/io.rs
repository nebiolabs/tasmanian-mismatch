//! Input helpers for reference FASTA and BAM-derived metadata.

use crate::types::ReferenceGenome;
use crate::types::{
    DiscountKey, GenomicMismatchKey, InconsistencyKey, InsertKey, MismatchKey, PositionMode,
};
use bio::io::fasta;
use rust_htslib::bam::{Read, Reader};
use std::collections::HashMap;
use std::error::Error;
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
    log::info!("Loading reference genome from: {}", fasta_path);
    let reader = fasta::Reader::from_file(fasta_path).expect("Failed to open reference FASTA file");

    let mut genome: ReferenceGenome = HashMap::new();
    for result in reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let chr_name = record.id().to_string();
        let sequence = record.seq().to_vec(); // Vec<u8> = byte, not UTF-8 char (overhead)
        genome.insert(chr_name, sequence);
    }

    log::info!("Loaded {} chromosome(s)", genome.len());
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
                log::error!("Failed to sample BAM file for read length estimation. Exiting.");
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

    log::info!("Computed max read length from sample: {}", max_length);
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

    log::info!(
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
    writeln!(
        writer,
        "read1_position\tread2_position\tdiscordance_type\tcount"
    )?;

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
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    write_mismatch_discounts_to_writer(mismatch_discounts, &mut writer)
}

pub fn write_mismatch_discounts_to_writer<W: Write>(
    mismatch_discounts: &HashMap<MismatchKey, usize>,
    writer: &mut W,
) -> std::io::Result<()> {
    let mut rows: Vec<(&MismatchKey, &usize)> = mismatch_discounts.iter().collect();
    rows.sort_by(|(a, _), (b, _)| {
        a.read_num
            .cmp(&b.read_num)
            .then(a.read_position.cmp(&b.read_position))
            .then(a.mismatch_type.cmp(&b.mismatch_type))
    });

    writeln!(
        writer,
        "mismatch_type\tread_num\tread_position\tdiscount_count"
    )?;

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
    if path == "-" {
        let stdin = std::io::stdin();
        return load_rescaling_matrix_from_reader(stdin.lock());
    }

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    load_rescaling_matrix_from_reader(reader)
}

pub fn load_rescaling_matrix_from_reader<R: BufRead>(
    reader: R,
) -> Result<HashMap<(u8, u16, char, char), f32>, Box<dyn Error>> {
    let mut matrix = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.trim().split('\t').collect();

        if parts.len() < 5 {
            continue;
        }

        // Skip header or malformed rows instead of failing the whole matrix load.
        let read_num = match parts[0].parse::<u8>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let position = match parts[1].parse::<u16>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let ref_base = match parts[2].chars().next() {
            Some(v) => v,
            None => continue,
        };
        let read_base = match parts[3].chars().next() {
            Some(v) => v,
            None => continue,
        };
        let scaling_factor = match parts[4].parse::<f32>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        matrix.insert((read_num, position, ref_base, read_base), scaling_factor);
    }

    Ok(matrix)
}

pub fn load_discount_table(path: &str) -> std::io::Result<HashMap<DiscountKey, usize>> {
    if path == "-" {
        let stdin = std::io::stdin();
        return load_discount_table_from_reader(stdin.lock());
    }

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    load_discount_table_from_reader(reader)
}

pub fn load_discount_table_from_reader<R: BufRead>(
    reader: R,
) -> std::io::Result<HashMap<DiscountKey, usize>> {
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

/// Extract the reference base from an `InsertKey`'s `base_change` field (e.g. `"C>T"` → `'C'`).
fn ref_base_of(key: &InsertKey) -> Option<char> {
    key.base_change
        .split_once('>')
        .and_then(|(ref_part, _)| ref_part.chars().next())
}

/// Normalize mismatch counts within each `(read_num, position, ref_base)` group.
///
/// For example, for a row key `C>T`, the normalized value is:
/// `C>T / (C>A + C>C + C>G + C>T)` at the same read number
/// and base position.
pub fn normalize_mismatch_counts(counts: &HashMap<InsertKey, usize>) -> HashMap<InsertKey, f64> {
    let mut group_totals: HashMap<(u8, usize, char), usize> = HashMap::new();

    for (key, count) in counts {
        if let Some(ref_base) = ref_base_of(key) {
            *group_totals
                .entry((key.read_num, key.base_position, ref_base))
                .or_insert(0) += count;
        }
    }

    counts
        .iter()
        .filter_map(|(key, &count)| {
            let ref_base = ref_base_of(key)?;
            let total = *group_totals.get(&(key.read_num, key.base_position, ref_base))?;
            // +1.0 guards against NaN when apply_external_discounts has reduced
            // every count in a group to zero.
            Some((key.clone(), count as f64 / (total as f64 + 1.0)))
        })
        .collect()
}

/// Convert normalized mismatch frequencies to a rescaling matrix format for tasmanian-rescale-quality.
///
/// This is a placeholder function that transforms normalized frequencies into scaling factors
/// suitable for quality score rescaling. The conversion strategy (e.g., inverse probability,
/// log-odds, empirical scaling) can be customized by changing the formula.
///
/// # Arguments
/// * `normalized_counts` - HashMap of normalized frequencies from [`normalize_mismatch_counts`].
///
/// # Returns
/// * A HashMap keyed by `(read_num, position, ref_base, read_base)` with scaling factors.
/// * Each scaling factor is currently a placeholder (1.0) pending implementation strategy.
///
/// # Implementation Notes
/// * The `position` is cast to `u16`; positions > 65535 will overflow.
/// * To use different scaling strategies (e.g., inverse probability, log-odds), modify
///   the conversion formula in the loop.
pub fn frequencies_to_rescaling_matrix(
    normalized_counts: &HashMap<InsertKey, f64>,
) -> HashMap<(u8, u16, char, char), f32> {
    let mut matrix: HashMap<(u8, u16, char, char), f32> = HashMap::new();

    for key in normalized_counts.keys() {
        let Some((ref_part, alt_part)) = key.base_change.split_once('>') else {
            continue;
        };

        let Some(ref_base) = ref_part.chars().next() else {
            continue;
        };
        let Some(alt_base) = alt_part.chars().next() else {
            continue;
        };

        let position = key.base_position as u16;

        // PLACEHOLDER: Default scaling factor = 1.0 (no scaling).
        // Replace this with actual conversion strategy:
        // - Example: 1.0 / normalized_freq (inverse: rarer variants get higher scaling)
        // - Example: (1.0 - normalized_freq) / normalized_freq (odds ratio)
        // - Example: Empirical quality recalibration based on mismatch-vs-quality correlation
        let scaling_factor = 1.0f32;

        let matrix_key = (key.read_num, position, ref_base, alt_base);
        matrix.insert(matrix_key, scaling_factor);
    }

    matrix
}

pub fn write_rescaling_matrix_output(
    counts: &HashMap<InsertKey, usize>,
    output_file: Option<&str>,
) -> std::io::Result<()> {
    let normalized = normalize_mismatch_counts(counts);
    let matrix = frequencies_to_rescaling_matrix(&normalized);

    let mut rows: Vec<(&(u8, u16, char, char), &f32)> = matrix.iter().collect();
    rows.sort_by(|(a, _), (b, _)| {
        a.0.cmp(&b.0)
            .then(a.1.cmp(&b.1))
            .then(a.2.cmp(&b.2))
            .then(a.3.cmp(&b.3))
    });

    with_output_writer(output_file, |w| {
        for ((read_num, position, ref_base, read_base), scaling_factor) in &rows {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{:.6}",
                read_num, position, ref_base, read_base, scaling_factor
            )?;
        }
        Ok(())
    })
}

fn position_label(mode: PositionMode) -> &'static str {
    match mode {
        PositionMode::Read => "read_position",
        PositionMode::Insert => "fragment_position",
    }
}

fn sort_insert_rows<V>(rows: &mut Vec<(&InsertKey, V)>) {
    rows.sort_by(|(a, _), (b, _)| {
        a.reference_order
            .cmp(&b.reference_order)
            .then(a.read_num.cmp(&b.read_num))
            .then(a.base_position.cmp(&b.base_position))
            .then(a.base_change.cmp(&b.base_change))
    });
}

/// Open `output_file` for writing, or lock stdout when `None`, then call `f`.
fn with_output_writer<F>(output_file: Option<&str>, f: F) -> std::io::Result<()>
where
    F: FnOnce(&mut dyn Write) -> std::io::Result<()>,
{
    if let Some(path) = output_file {
        f(&mut BufWriter::new(File::create(path)?))
    } else {
        f(&mut std::io::stdout().lock())
    }
}

pub fn write_output(
    counts: &HashMap<InsertKey, usize>,
    output_file: Option<&str>,
    position_mode: PositionMode,
) -> std::io::Result<()> {
    let label = position_label(position_mode);
    let mut rows: Vec<(&InsertKey, &usize)> = counts.iter().collect();
    sort_insert_rows(&mut rows);

    with_output_writer(output_file, |w| {
        writeln!(
            w,
            "base_change\tread_num\treference_order\t{}\tcount",
            label
        )?;
        for (key, count) in &rows {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                key.base_change, key.read_num, key.reference_order, key.base_position, count
            )?;
        }
        Ok(())
    })
}

pub fn write_normalized_output(
    counts: &HashMap<InsertKey, usize>,
    output_file: Option<&str>,
    position_mode: PositionMode,
) -> std::io::Result<()> {
    let normalized = normalize_mismatch_counts(counts);
    let label = position_label(position_mode);
    let mut rows: Vec<(&InsertKey, &f64)> = normalized.iter().collect();
    sort_insert_rows(&mut rows);

    with_output_writer(output_file, |w| {
        writeln!(
            w,
            "base_change\tread_num\treference_order\t{}\tnormalized_frequency",
            label
        )?;
        for (key, freq) in &rows {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{:.6}",
                key.base_change, key.read_num, key.reference_order, key.base_position, freq
            )?;
        }
        Ok(())
    })
}
