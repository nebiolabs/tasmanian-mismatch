/// Return the DNA complement of a base
pub fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => base,
    }
}

/// Convert a byte to a DNA base character
pub fn base_to_char(byte: u8) -> Option<char> {
    match byte {
        b'A' => Some('A'),
        b'C' => Some('C'),
        b'G' => Some('G'),
        b'T' => Some('T'),
        _ => None,
    }
}

/// Adjust read position based on mode length for normalization --> INSTEAD OF THIS, USE THE START OF THE READ + TLEN. CALL IT REFERENCE_SPAN
pub fn correct_read_len_with_mode(
    read_pos: usize,
    seq_len: usize,
    mode_len: usize
) -> usize {
    if mode_len > 0 && read_pos > seq_len/2 {  
        read_pos + (mode_len - seq_len)
    } else {
        read_pos
    }
}

/// Calculate genomic end position from CIGAR string
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

/// Parse MD tag (currently unused but kept for potential future use)
#[allow(dead_code)]
pub fn parse_md_tag(md_string: &str) -> (Vec<(usize, char)>, Vec<(usize, usize)>) {
    // Parse MD tag and return:
    // 1. (position_offset, ref_base) for mismatches
    // 2. (start_pos, length) for matches
    let mut mismatches = Vec::new();
    let mut matches = Vec::new();
    let mut pos = 0;
    let mut num_str = String::new();
    
    for ch in md_string.chars() {
        if ch.is_numeric() {
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
        } else if ch.is_alphabetic() {
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