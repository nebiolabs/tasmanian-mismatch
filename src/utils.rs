//! General sequence, coordinate, and output helpers.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Mismatches: list of (read_position, reference_base).
type MdMismatches = Vec<(usize, char)>;
/// Runs of matching bases: list of (start_position, length).
type MdMatches = Vec<(usize, usize)>;

/// Return the DNA complement of a nucleotide base.
///
/// # Arguments
/// * `base` - DNA base to complement.
///
/// # Returns
/// * The complemented base for `A`, `C`, `G`, or `T`.
/// * The original input for any other character.
pub fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => base,
    }
}

/// Convert an encoded base byte to an uppercase DNA character.
///
/// # Arguments
/// * `byte` - Base byte to decode.
///
/// # Returns
/// * `Some(char)` for canonical `A`, `C`, `G`, or `T` bases.
/// * `None` for any other byte.
pub fn base_to_char(byte: u8) -> Option<char> {
    match byte {
        b'A' => Some('A'),
        b'C' => Some('C'),
        b'G' => Some('G'),
        b'T' => Some('T'),
        _ => None,
    }
}

/// Normalize a read position using the modal read length.
///
/// # Arguments
/// * `read_pos` - Position within the read.
/// * `seq_len` - Observed read length.
/// * `mode_len` - Modal read length used for right-half adjustment.
///
/// # Returns
/// * The original position when no adjustment is needed.
/// * A mode-adjusted position for bases in the right half of shorter reads.
// instead of this we can use start of read + tlen and call it reference_span
pub fn correct_read_len_with_mode(
    read_pos: usize,
    seq_len: usize,
    mode_len: usize,
    use_insert_mode: bool,
    read_num: u8,
) -> usize {
    if use_insert_mode {
        match read_num {
            1 => read_pos,
            2 => (2 * mode_len + 10) - (seq_len - read_pos),
            _ => read_pos,
        }
    } else {
        if mode_len > 0 && read_pos > seq_len / 2 {
            read_pos + (mode_len - seq_len)
        } else {
            read_pos
        }
    }
}

/// Calculate the reference end position implied by a CIGAR string.
///
/// # Arguments
/// * `start_pos` - Alignment start position in 0-based coordinates.
/// * `cigar` - CIGAR string for the alignment.
///
/// # Returns
/// * The exclusive genomic end position after consuming reference bases.
pub fn calculate_end_pos(start_pos: i64, cigar: &rust_htslib::bam::record::CigarStringView) -> i64 {
    use rust_htslib::bam::record::Cigar::*;
    let mut end = start_pos;
    for op in cigar.iter() {
        match op {
            Match(len) | Equal(len) | Diff(len) | Del(len) | RefSkip(len) => {
                end += *len as i64;
            }
            _ => {}
        }
    }
    end
}

/// Print the position-based mismatch table as CSV.
///
/// # Arguments
/// * `position_map` - Nested mismatch counts keyed by read number and position.
/// * `output_path` - Optional output path. When `None`, writes to stdout.
pub fn print_position_table(
    position_map: &HashMap<(u8, usize), HashMap<String, usize>>,
    output_path: Option<&str>,
) -> std::io::Result<()> {
    fn write_position_table<W: Write>(
        writer: &mut W,
        position_map: &HashMap<(u8, usize), HashMap<String, usize>>,
    ) -> std::io::Result<()> {
        let mut all_mismatch_types: Vec<String> = position_map
            .values()
            .flat_map(|m| m.keys().cloned())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        all_mismatch_types.sort();

        let mut positions: Vec<(u8, usize)> = position_map.keys().copied().collect();
        positions.sort();

        write!(writer, "Read,Position")?;
        for mismatch_type in &all_mismatch_types {
            write!(writer, ",{}", mismatch_type)?;
        }
        writeln!(writer)?;

        for &(read_num, pos) in &positions {
            write!(writer, "{},{}", read_num, pos)?;

            if let Some(mismatch_counts) = position_map.get(&(read_num, pos)) {
                for mismatch_type in &all_mismatch_types {
                    let count = mismatch_counts.get(mismatch_type).unwrap_or(&0);
                    write!(writer, ",{}", count)?;
                }
            }
            writeln!(writer)?;
        }

        Ok(())
    }

    if let Some(path) = output_path {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        write_position_table(&mut writer, position_map)?;
    } else {
        let stdout = std::io::stdout();
        let mut writer = stdout.lock();
        write_position_table(&mut writer, position_map)?;
    }

    Ok(())
}

/// Parse an MD tag string into mismatches and match runs.
///
/// # Arguments
/// * `md` - The MD tag string from a BAM record.
/// * `read_name` - Read name used in panic messages to help locate the offending record.
///
/// # Returns
/// * A tuple of `(mismatches, matches)` where:
///   - `mismatches` is a `Vec<(usize, char)>` of (read_position, reference_base)
///   - `matches` is a `Vec<(usize, usize)>` of (start_position, length)
pub fn parse_md_tag(md: &str, read_name: &str) -> (MdMismatches, MdMatches) {
    let mut mismatches = Vec::new();
    let mut matches = Vec::new();
    let mut pos = 0usize;
    let mut chars = md.chars().peekable();

    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() {
            // digits: run of matching bases e.g. "5" in "5A4"
            let len: usize = std::iter::from_fn(|| chars.next_if(|c| c.is_ascii_digit()))
                .collect::<String>()
                .parse()
                .unwrap();
            if len > 0 {
                matches.push((pos, len));
                pos += len;
            }
        } else if c == '^' {
            // caret: deleted reference bases e.g. "^AC" — skip, no read position advance
            chars.next(); // consume '^'
            while chars.next_if(|c| c.is_ascii_alphabetic()).is_some() {}
        } else if c.is_ascii_alphabetic() {
            // letter: substitution — reference base differs from read
            chars.next();
            mismatches.push((pos, c.to_ascii_uppercase()));
            pos += 1;
        } else {
            panic!(
                "unexpected character '{}' in MD tag '{}' for read '{}'",
                c, md, read_name
            );
        }
    }

    (mismatches, matches)
}
///
/// # Arguments
/// * `position_map` - Nested mismatch counts keyed by read number and position.
/// * `output_path` - Optional output path. When `None`, writes to stdout.
pub fn print_main_output(
    position_map: &HashMap<(u8, usize), HashMap<String, usize>>,
    output_path: Option<&str>,
) -> std::io::Result<()> {
    print_position_table(position_map, output_path)?;
    log::info!(
        "Total unique (read, position) combinations: {}",
        position_map.len()
    );

    Ok(())
}
