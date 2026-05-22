//! Methylation-aware base normalization.

/// Adjust a read base for methylation-aware mismatch analysis.
///
/// Converts bisulfite-style `C>T` events in read 1 and `G>A` events in read 2
/// back to the corresponding reference base when methylation mode is enabled.
///
/// # Arguments
/// * `read_base` - Read base after any strand normalization.
/// * `ref_base` - Reference base after any strand normalization.
/// * `read_num` - Read number within the pair, usually `1` or `2`.
/// * `is_methylation` - Whether methylation-aware conversion is enabled.
///
/// # Returns
/// * The base to use for mismatch accounting.
pub fn adjust_methylation_base(
    read_base: char,
    ref_base: char,
    read_num: u8,
    is_methylation: bool,
) -> char {
    if !is_methylation {
        return read_base;
    }

    // Convert all C>T and G>A
    match (read_num, ref_base, read_base) {
        // Read 1: C in ref, T in read -> treat as C (unmethylated C)
        (1, 'C', 'T') => 'C',
        // Read 2: G in ref, A in read -> treat as G (unmethylated C on reverse strand)
        (2, 'G', 'A') => 'G',
        _ => read_base,
    }
}
