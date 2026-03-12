/// Adjust base calling for methylation-aware analysis
/// Converts C>T (read1) and G>A (read2) in appropriate contexts
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
