//! General sequence, coordinate, and output helpers.

use std::collections::HashMap;

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
pub fn correct_read_len_with_mode(read_pos: usize, seq_len: usize, mode_len: usize) -> usize {
    if mode_len > 0 && read_pos > seq_len / 2 {
        read_pos + (mode_len - seq_len)
    } else {
        read_pos
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

/// Parse an MD tag into mismatch and match spans.
///
/// # Arguments
/// * `md_string` - MD tag string from a SAM/BAM alignment record.
///
/// # Returns
/// * A tuple containing mismatch positions and match spans.
#[allow(dead_code)]
pub fn parse_md_tag(md_string: &str) -> (Vec<(usize, char)>, Vec<(usize, usize)>) {
    // Parse MD tag and return:
    // 1. (position_offset, ref_base) for mismatches
    // 2. (start_pos, length) for matches
    let mut mismatches = Vec::new();
    let mut matches = Vec::new();
    let mut pos = 0;
    let mut num_str = String::new();
    let mut in_deletion = false;

    for ch in md_string.chars() {
        if ch.is_numeric() {
            if in_deletion {
                in_deletion = false;
            }
            num_str.push(ch);
        } else if ch == '^' {
            // Deletion - skip the deleted bases
            if !num_str.is_empty() {
                let match_len = num_str.parse::<usize>().unwrap_or(0);
                if match_len > 0 {
                    matches.push((pos, match_len));
                    pos += match_len;
                }
                num_str.clear();
            }
            in_deletion = true;
        } else if ch.is_alphabetic() {
            if in_deletion {
                // Deletion bases are present in reference only and should not affect read position.
                continue;
            }
            // This is a reference base at a mismatch position
            if !num_str.is_empty() {
                let match_len = num_str.parse::<usize>().unwrap_or(0);
                if match_len > 0 {
                    matches.push((pos, match_len));
                    pos += match_len;
                }
                num_str.clear();
            }
            mismatches.push((pos, ch));
            pos += 1;
        }
    }

    // Handle trailing matches
    if !num_str.is_empty() {
        let match_len = num_str.parse::<usize>().unwrap_or(0);
        if match_len > 0 {
            matches.push((pos, match_len));
        }
    }

    (mismatches, matches)
}

/// Print the position-based mismatch table as CSV.
///
/// # Arguments
/// * `position_map` - Nested mismatch counts keyed by read number and position.
pub fn print_position_table(position_map: &HashMap<(u8, usize), HashMap<String, usize>>) {
    let mut all_mismatch_types: Vec<String> = position_map
        .values()
        .flat_map(|m| m.keys().cloned())
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    all_mismatch_types.sort();

    let mut positions: Vec<(u8, usize)> = position_map.keys().copied().collect();
    positions.sort();

    print!("Read,Position");
    for mismatch_type in &all_mismatch_types {
        print!(",{}", mismatch_type);
    }
    println!();

    for &(read_num, pos) in &positions {
        print!("{},{}", read_num, pos);

        if let Some(mismatch_counts) = position_map.get(&(read_num, pos)) {
            for mismatch_type in &all_mismatch_types {
                let count = mismatch_counts.get(mismatch_type).unwrap_or(&0);
                print!(",{}", count);
            }
        }
        println!();
    }
}

/// Print the main mismatch table followed by summary statistics.
///
/// # Arguments
/// * `position_map` - Nested mismatch counts keyed by read number and position.
pub fn print_main_output(position_map: &HashMap<(u8, usize), HashMap<String, usize>>) {
    print_position_table(position_map);
    eprintln!(
        "\nTotal unique (read, position) combinations: {}",
        position_map.len()
    );
}
