//! Methylation-aware base normalization.

/// Adjust a read base for methylation-aware mismatch analysis.
///
/// Depending on configuration, this converts bisulfite-style `C>T` events in
/// read 1 and `G>A` events in read 2 back to the corresponding reference base,
/// optionally restricting that conversion to CpG context.
///
/// # Arguments
/// * `read_base` - Read base after any strand normalization.
/// * `ref_base` - Reference base after any strand normalization.
/// * `read_num` - Read number within the pair, usually `1` or `2`.
/// * `is_methylation` - Whether methylation-aware conversion is enabled.
/// * `cpg_only` - Whether conversion should be restricted to CpG context.
/// * `ref_seq` - Reference sequence for the current chromosome.
/// * `genome_pos` - Genomic position within `ref_seq`.
///
/// # Returns
/// * The base to use for mismatch accounting.
pub fn adjust_methylation_base(
    read_base: char,
    ref_base: char,
    read_num: u8,
    is_methylation: bool,
    cpg_only: bool,
    ref_seq: &[u8],
    genome_pos: usize,
) -> char {
    if !is_methylation {
        return read_base;
    }

    if cpg_only {
        // Only convert in CpG context
        match (read_num, ref_base, read_base) {
            // Read 1: C in ref, next base is G, T in read -> treat as C (unmethylated CpG)
            (1, 'C', 'T') => {
                if genome_pos + 1 < ref_seq.len() {
                    let next_base = ref_seq[genome_pos + 1];
                    if next_base == b'G' || next_base == b'g' {
                        'C'
                    } else {
                        read_base
                    }
                } else {
                    read_base
                }
            }
            // Read 2: G in ref, next base is C, A in read -> treat as G (unmethylated CpG on reverse)
            (2, 'G', 'A') => {
                if genome_pos + 1 < ref_seq.len() {
                    let next_base = ref_seq[genome_pos + 1];
                    if next_base == b'C' || next_base == b'c' {
                        'G'
                    } else {
                        read_base
                    }
                } else {
                    read_base
                }
            }
            _ => read_base,
        }
    } else {
        // Convert all C>T and G>A
        match (read_num, ref_base, read_base) {
            // Read 1: C in ref, T in read -> treat as C (unmethylated C)
            (1, 'C', 'T') => 'C',
            // Read 2: G in ref, A in read -> treat as G (unmethylated C on reverse strand)
            (2, 'G', 'A') => 'G',
            _ => read_base,
        }
    }
}
