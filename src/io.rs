//! Input helpers for reference FASTA and BAM-derived metadata.

use crate::types::{
    DiscountKey, GenomicMismatchKey, InconsistencyKey, InsertKey, MismatchKey, PositionMode,
    ReferenceGenome, ReferenceOrder, RescalingMatrix,
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
pub fn load_rescaling_matrix(path: &str) -> Result<RescalingMatrix, Box<dyn Error>> {
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
) -> Result<RescalingMatrix, Box<dyn Error>> {
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
            reference_order: ReferenceOrder::First,
        };
        let k2 = InsertKey {
            base_change: discount_key.base_change,
            read_num: discount_key.read_num,
            base_position: discount_key.base_position,
            reference_order: ReferenceOrder::Second,
        };

        let c1 = counts.get(&k1).copied().unwrap_or(0);
        let c2 = counts.get(&k2).copied().unwrap_or(0);

        if c1 == 0 && c2 == 0 {
            continue;
        }

        let first_key = if c1 >= c2 { &k1 } else { &k2 };
        if remaining > 0
            && let Some(v) = counts.get_mut(first_key) {
                let take = remaining.min(*v);
                *v -= take;
                remaining -= take;
            }

        if remaining > 0 {
            let second_key = if first_key.reference_order == ReferenceOrder::First {
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
            // If total is 0 (which can occur if apply_external_discounts has reduced
            // every count in a group to zero), return 0.0 to guard against NaN.
            let frequency = if total == 0 {
                0.0
            } else {
                count as f64 / total as f64
            };
            Some((key.clone(), frequency))
        })
        .collect()
}

/// Fraction of the read length excluded from each end when computing a
/// mismatch class's baseline rate (systematic damage, like FFPE, concentrates
/// near read ends, so the ends are excluded from the "no excess signal" floor).
const RESCALING_EDGE_FRACTION: f64 = 0.2;

/// Floor for scaling factors, so a hot position never fully zeroes quality.
const RESCALING_MIN_SCALING: f32 = 0.05;

/// Pseudocount used to shrink a position's scaling factor back toward 1.0
/// (no discount) when few mismatches actually back it. A position/class bin
/// with a handful of real events can look just as "elevated" as one with a
/// genuine, well-supported damage pattern once turned into a frequency --
/// the group depth doesn't help distinguish them (it pools reads from every
/// genomic locus that happens to land on that read position/class, and is
/// large even when the actual mismatch count is tiny). Weighting by the raw
/// count lets a handful of coincidental true-variant mismatches be shrunk
/// back toward "no discount" while a systematic pattern backed by many
/// mismatches (e.g. FFPE deamination near read ends) still gets discounted.
const RESCALING_SHRINKAGE_PSEUDOCOUNT: f64 = 20.0;

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

/// Build a dense per-position array over `[min_pos, max_pos]`, linearly
/// interpolating positions that fall between two observed points and
/// filling from the nearest observed value at either end. A single position
/// Tasmanian never reported (e.g. below `--min-base-quality` depth) would
/// otherwise create a smoothing discontinuity.
fn densify_and_interpolate(
    positions: &[usize],
    freqs: &[f64],
    min_pos: usize,
    max_pos: usize,
) -> Vec<f64> {
    let n = max_pos - min_pos + 1;
    let mut dense: Vec<Option<f64>> = vec![None; n];
    for (&p, &f) in positions.iter().zip(freqs.iter()) {
        dense[p - min_pos] = Some(f);
    }

    let mut result = vec![0.0; n];
    let mut last_known: Option<(usize, f64)> = None;
    let mut i = 0;
    while i < n {
        if let Some(v) = dense[i] {
            result[i] = v;
            last_known = Some((i, v));
            i += 1;
            continue;
        }

        let mut j = i;
        while j < n && dense[j].is_none() {
            j += 1;
        }
        let next_known = if j < n { dense[j].map(|v| (j, v)) } else { None };

        match (last_known, next_known) {
            (Some((li, lv)), Some((ni, nv))) => {
                for (k, slot) in result.iter_mut().enumerate().take(j).skip(i) {
                    let t = (k - li) as f64 / (ni - li) as f64;
                    *slot = lv + t * (nv - lv);
                }
            }
            (None, Some((_, nv))) => {
                for slot in result.iter_mut().take(j).skip(i) {
                    *slot = nv;
                }
            }
            (Some((_, lv)), None) => {
                for slot in result.iter_mut().take(j).skip(i) {
                    *slot = lv;
                }
            }
            (None, None) => {} // no observations at all in this class
        }
        i = j;
    }
    result
}

/// Centered rolling median with `min_periods = 1` (the window shrinks near
/// either edge instead of leaving those positions unsmoothed).
fn rolling_median(values: &[f64], window: usize) -> Vec<f64> {
    let n = values.len();
    let half = window / 2;
    (0..n)
        .map(|i| {
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            median(&values[lo..hi])
        })
        .collect()
}

/// Convert normalized mismatch frequencies to a rescaling matrix format for tasmanian-rescale-quality.
///
/// Per-position mismatch frequencies are noisy on their own (sampling noise
/// scales with local depth), so using them directly as scaling factors would
/// just re-inject noise into the BAM. Instead, each `(read_num, ref_base,
/// alt_base)` mismatch class's frequency profile across read position is:
///
///   1. Collapsed across `reference_order` (which mate maps first) since the
///      rescaling matrix and [`rescale_phred_scores`](crate::rescale_phred_scores)
///      key only on `read_num`/position/bases. Because [`normalize_mismatch_counts`]
///      divides each row by a `(read_num, position, ref_base)` total that
///      already pools both `reference_order` values, summing the two
///      normalized frequencies recovers the combined-order frequency.
///   2. Smoothed with a centered rolling median -- the "model" -- to
///      separate systematic positional bias (e.g. FFPE damage rising at
///      read ends) from sampling noise.
///   3. Compared against that class's baseline rate (the median smoothed
///      rate over the central read, away from the edges where damage
///      concentrates): positions at or below baseline get a scaling factor
///      of 1.0 (left untouched); positions with an elevated smoothed rate
///      -- more likely artifact than true variant -- get a proportionally
///      reduced factor, floored at [`RESCALING_MIN_SCALING`].
///   4. Shrunk back toward 1.0 in proportion to how few raw mismatches
///      actually back that position/class (see
///      [`RESCALING_SHRINKAGE_PSEUDOCOUNT`]), so a position that only looks
///      "elevated" because a couple of true variants happened to land there
///      isn't discounted as if it were a well-supported damage signal.
///
/// # Arguments
/// * `normalized_counts` - HashMap of normalized frequencies from [`normalize_mismatch_counts`].
/// * `raw_counts` - The same raw mismatch counts `normalized_counts` was derived from
///   (i.e. what was passed to [`normalize_mismatch_counts`]), used only to weigh how much
///   to trust each position's deviation from baseline -- never to change the profile shape.
///
/// # Returns
/// * A HashMap keyed by `(read_num, position, ref_base, read_base)` with scaling factors.
///
/// # Implementation Notes
/// * The `position` is cast to `u16`; positions > 65535 will overflow.
pub fn frequencies_to_rescaling_matrix(
    normalized_counts: &HashMap<InsertKey, f64>,
    raw_counts: &HashMap<InsertKey, usize>,
) -> RescalingMatrix {
    // Collapse reference_order into (read_num, position, ref_base, alt_base),
    // pooling both the normalized frequency and the raw mismatch count the
    // same way (summing) since each reference_order row is an independent
    // slice of the same group.
    let mut collapsed: HashMap<(u8, usize, char, char), (f64, usize)> = HashMap::new();
    for (key, &freq) in normalized_counts {
        let Some((ref_part, alt_part)) = key.base_change.split_once('>') else {
            continue;
        };
        let Some(ref_base) = ref_part.chars().next() else {
            continue;
        };
        let Some(alt_base) = alt_part.chars().next() else {
            continue;
        };
        if ref_base == alt_base {
            continue; // matches aren't scaled; rescale_phred_scores only looks up mismatches
        }

        let raw_count = raw_counts.get(key).copied().unwrap_or(0);
        let entry = collapsed
            .entry((key.read_num, key.base_position, ref_base, alt_base))
            .or_insert((0.0, 0));
        entry.0 += freq;
        entry.1 += raw_count;
    }

    // Group into per-(read_num, ref_base, alt_base) positional series.
    let mut by_class: HashMap<(u8, char, char), Vec<(usize, f64, usize)>> = HashMap::new();
    for (&(read_num, pos, ref_base, alt_base), &(freq, count)) in &collapsed {
        by_class
            .entry((read_num, ref_base, alt_base))
            .or_default()
            .push((pos, freq, count));
    }

    let mut matrix: RescalingMatrix = HashMap::new();

    for ((read_num, ref_base, alt_base), mut series) in by_class {
        series.sort_by_key(|&(pos, _, _)| pos);
        let positions: Vec<usize> = series.iter().map(|&(p, _, _)| p).collect();
        let freqs: Vec<f64> = series.iter().map(|&(_, f, _)| f).collect();
        let counts: Vec<usize> = series.iter().map(|&(_, _, c)| c).collect();

        let min_pos = *positions.first().unwrap();
        let max_pos = *positions.last().unwrap();
        let n_positions = max_pos - min_pos + 1;

        let window = (n_positions / 10).max(5);
        let window = if window % 2 == 0 { window + 1 } else { window };

        let dense = densify_and_interpolate(&positions, &freqs, min_pos, max_pos);
        let smoothed = rolling_median(&dense, window);

        let lo = min_pos as f64 + RESCALING_EDGE_FRACTION * n_positions as f64;
        let hi = max_pos as f64 - RESCALING_EDGE_FRACTION * n_positions as f64;
        let central: Vec<f64> = (min_pos..=max_pos)
            .filter(|&p| (p as f64) >= lo && (p as f64) <= hi)
            .map(|p| smoothed[p - min_pos])
            .collect();
        let baseline = if central.is_empty() {
            median(&smoothed)
        } else {
            median(&central)
        };

        const EPS: f64 = 1e-9;
        for (i, &pos) in positions.iter().enumerate() {
            let sm = smoothed[pos - min_pos];
            let raw_scaling = baseline / sm.max(EPS);

            // Trust the deviation from baseline in proportion to how many
            // actual mismatches this position/class contributed -- not the
            // (always large) pooled depth -- so low-count noise relaxes
            // toward 1.0 instead of being scaled as if it were real signal.
            let count = counts[i] as f64;
            let weight = count / (count + RESCALING_SHRINKAGE_PSEUDOCOUNT);
            let shrunk_scaling = 1.0 - weight * (1.0 - raw_scaling);

            let scaling_factor = shrunk_scaling.clamp(RESCALING_MIN_SCALING as f64, 1.0) as f32;
            matrix.insert((read_num, pos as u16, ref_base, alt_base), scaling_factor);
        }
    }

    matrix
}

pub fn write_rescaling_matrix_output(
    counts: &HashMap<InsertKey, usize>,
    output_file: Option<&str>,
) -> std::io::Result<()> {
    let normalized = normalize_mismatch_counts(counts);
    let matrix = frequencies_to_rescaling_matrix(&normalized, counts);

    let mut rows: Vec<_> = matrix.iter().collect();
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

/// Launch the embedded Bokeh visualization script on the given TSV file.
/// Writes the Python script to a temp file and invokes `python3` on it.
pub fn launch_visualization(tsv_path: &str) -> std::io::Result<()> {
    let script = if std::path::Path::new("scripts/visualize.py").exists() {
        match std::fs::read_to_string("scripts/visualize.py") {
            Ok(s) => {
                log::info!("Using live visualization script from scripts/visualize.py");
                s
            }
            Err(e) => {
                log::warn!(
                    "Could not read scripts/visualize.py ({}), using embedded fallback",
                    e
                );
                include_str!("../scripts/visualize.py").to_string()
            }
        }
    } else {
        include_str!("../scripts/visualize.py").to_string()
    };
    let tmp_dir = std::env::temp_dir();
    let script_path = tmp_dir.join("tasmanian_visualize.py");
    std::fs::write(&script_path, script)?;

    let html_path = format!(
        "{}.html",
        tsv_path.trim_end_matches(".tsv").trim_end_matches(".csv")
    );

    log::info!("Launching visualization: {} -> {}", tsv_path, html_path);

    let status = std::process::Command::new("python3")
        .arg(&script_path)
        .arg(tsv_path)
        .arg("-o")
        .arg(&html_path)
        .status();

    match status {
        Ok(s) if s.success() => {
            log::info!("Visualization saved to {}", html_path);
            Ok(())
        }
        Ok(s) => {
            log::warn!("Visualization script exited with {}", s);
            Ok(())
        }
        Err(e) => {
            log::warn!(
                "Could not launch visualization (python3 not found or failed): {}. \
                 Install bokeh and pandas: pip install bokeh pandas",
                e
            );
            Ok(())
        }
    }
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

#[cfg(test)]
mod rescaling_matrix_tests {
    use super::*;

    fn insert_key(base_change: &str, read_num: u8, base_position: usize) -> InsertKey {
        InsertKey {
            base_change: base_change.to_string(),
            read_num,
            base_position,
            reference_order: ReferenceOrder::First,
        }
    }

    /// Build a normalized-frequency map for a 76bp read where C>T is
    /// elevated near both ends (mimicking FFPE damage) and flat elsewhere,
    /// with every other mismatch class flat throughout.
    fn ffpe_like_normalized_counts() -> HashMap<InsertKey, f64> {
        let mut counts = HashMap::new();
        for pos in 0..76usize {
            let edge = pos.min(75 - pos);
            let ct_freq = if edge < 10 { 0.24 } else { 0.015 };
            counts.insert(insert_key("C>T", 1, pos), ct_freq);
            counts.insert(insert_key("C>A", 1, pos), 0.02);
            counts.insert(insert_key("C>G", 1, pos), 0.02);
            counts.insert(insert_key("C>C", 1, pos), 1.0 - ct_freq - 0.04);
        }
        counts
    }

    /// Raw mismatch counts consistent with [`ffpe_like_normalized_counts`],
    /// but well-supported (hundreds of observations per bin) so the edge
    /// signal there is trusted at close to full weight -- i.e. what a real,
    /// systematic damage pattern actually backed by data looks like.
    fn ffpe_like_raw_counts() -> HashMap<InsertKey, usize> {
        let mut counts = HashMap::new();
        for pos in 0..76usize {
            let edge = pos.min(75 - pos);
            let ct_count = if edge < 10 { 240 } else { 15 };
            counts.insert(insert_key("C>T", 1, pos), ct_count);
            counts.insert(insert_key("C>A", 1, pos), 20);
            counts.insert(insert_key("C>G", 1, pos), 20);
            counts.insert(insert_key("C>C", 1, pos), 1000);
        }
        counts
    }

    #[test]
    fn downweights_elevated_edge_positions_only() {
        let matrix = frequencies_to_rescaling_matrix(
            &ffpe_like_normalized_counts(),
            &ffpe_like_raw_counts(),
        );

        let edge_factor = matrix[&(1, 0, 'C', 'T')];
        let center_factor = matrix[&(1, 38, 'C', 'T')];

        assert!(
            edge_factor < 0.5,
            "expected the elevated read-end C>T rate to be downweighted, got {edge_factor}"
        );
        assert!(
            (center_factor - 1.0).abs() < 1e-6,
            "expected the flat central-read rate to be left at 1.0, got {center_factor}"
        );
    }

    #[test]
    fn matches_are_never_scaled() {
        let matrix = frequencies_to_rescaling_matrix(
            &ffpe_like_normalized_counts(),
            &ffpe_like_raw_counts(),
        );
        assert!(!matrix.contains_key(&(1, 0, 'C', 'C')));
    }

    #[test]
    fn low_raw_counts_are_shrunk_toward_one_even_with_elevated_frequency() {
        // The exact same elevated-frequency profile as the FFPE-like test
        // (edges "hot", center flat) -- but each bin is backed by only a
        // couple of raw mismatches instead of hundreds, i.e. the handful of
        // true variants a real "no noise" condition still carries. The rate
        // looks identical to genuine damage; only the count tells them
        // apart, and this used to get discounted just as confidently.
        let normalized = ffpe_like_normalized_counts();
        let mut sparse_raw = HashMap::new();
        for pos in 0..76usize {
            let edge = pos.min(75 - pos);
            let ct_count = if edge < 10 { 2 } else { 1 };
            sparse_raw.insert(insert_key("C>T", 1, pos), ct_count);
            sparse_raw.insert(insert_key("C>A", 1, pos), 1);
            sparse_raw.insert(insert_key("C>G", 1, pos), 1);
            sparse_raw.insert(insert_key("C>C", 1, pos), 40);
        }

        let matrix = frequencies_to_rescaling_matrix(&normalized, &sparse_raw);
        let edge_factor = matrix[&(1, 0, 'C', 'T')];

        assert!(
            edge_factor > 0.9,
            "a handful of mismatches shouldn't be discounted as confidently as a well-supported pattern, got {edge_factor}"
        );
    }

    #[test]
    fn scaling_factor_is_floored() {
        let mut normalized = HashMap::new();
        let mut raw = HashMap::new();
        // One position with a huge excess mismatch rate relative to
        // baseline, backed by enough raw mismatches to be trusted at
        // (near) full weight.
        normalized.insert(insert_key("A>G", 2, 0), 0.9);
        normalized.insert(insert_key("A>G", 2, 1), 0.9);
        raw.insert(insert_key("A>G", 2, 0), 900);
        raw.insert(insert_key("A>G", 2, 1), 900);
        for pos in 2..30usize {
            normalized.insert(insert_key("A>G", 2, pos), 0.001);
            raw.insert(insert_key("A>G", 2, pos), 50);
        }

        let matrix = frequencies_to_rescaling_matrix(&normalized, &raw);
        assert_eq!(matrix[&(2, 0, 'A', 'G')], RESCALING_MIN_SCALING);
    }

    #[test]
    fn reference_order_is_pooled() {
        // Two reference_order rows for the same (read_num, position,
        // base_change), as normalize_mismatch_counts would emit for a
        // pooled denominator: their frequencies must sum, not overwrite.
        let mut normalized = HashMap::new();
        let mut raw = HashMap::new();
        for pos in 0..20usize {
            let key_first = InsertKey {
                base_change: "A>G".to_string(),
                read_num: 1,
                base_position: pos,
                reference_order: ReferenceOrder::First,
            };
            let key_second = InsertKey {
                base_change: "A>G".to_string(),
                read_num: 1,
                base_position: pos,
                reference_order: ReferenceOrder::Second,
            };
            normalized.insert(key_first.clone(), 0.005);
            normalized.insert(key_second.clone(), 0.005);
            raw.insert(key_first, 5);
            raw.insert(key_second, 5);
        }

        let matrix = frequencies_to_rescaling_matrix(&normalized, &raw);
        // Flat profile -> every position at its own class baseline -> 1.0,
        // regardless of the per-reference_order split.
        assert!((matrix[&(1, 5, 'A', 'G')] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn densify_and_interpolate_fills_gaps_and_edges() {
        let positions = [2usize, 5, 8];
        let freqs = [0.0, 1.0, 0.5];
        let dense = densify_and_interpolate(&positions, &freqs, 0, 8);

        assert_eq!(dense.len(), 9);
        assert_eq!(dense[0], 0.0); // back-filled from position 2
        assert_eq!(dense[2], 0.0);
        assert!((dense[3] - (1.0 / 3.0)).abs() < 1e-9); // interpolated toward position 5
        assert_eq!(dense[5], 1.0);
        assert_eq!(dense[8], 0.5);
    }

    #[test]
    fn rolling_median_shrinks_window_at_edges() {
        let values = [1.0, 2.0, 3.0, 100.0, 5.0];
        let smoothed = rolling_median(&values, 3);
        // Center point sees [2,3,100] -> median 3.0, unaffected by the spike's magnitude.
        assert_eq!(smoothed[2], 3.0);
        // Left edge only has [1,2] available (min_periods=1 semantics).
        assert_eq!(smoothed[0], median(&[1.0, 2.0]));
    }
}
