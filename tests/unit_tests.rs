// Unit tests for tasmanian-mismatch

use rustmanian_mismatch::*;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_complement() {
        assert_eq!(complement('A'), 'T');
        assert_eq!(complement('T'), 'A');
        assert_eq!(complement('C'), 'G');
        assert_eq!(complement('G'), 'C');
        assert_eq!(complement('N'), 'N');
    }
    
    #[test]
    fn test_base_to_char() {
        assert_eq!(base_to_char(b'A'), Some('A'));
        assert_eq!(base_to_char(b'C'), Some('C'));
        assert_eq!(base_to_char(b'G'), Some('G'));
        assert_eq!(base_to_char(b'T'), Some('T'));
        assert_eq!(base_to_char(b'N'), None);
        assert_eq!(base_to_char(b'X'), None);
    }
    
    #[test]
    fn test_correct_read_len_with_mode() {
        // No mode correction
        assert_eq!(correct_read_len_with_mode(10, 100, 0), 10);
        
        // With mode correction for second half of read
        assert_eq!(correct_read_len_with_mode(60, 100, 150), 110);
        
        // No correction for first half even with mode
        assert_eq!(correct_read_len_with_mode(40, 100, 150), 40);
    }
    
    #[test]
    fn test_calculate_end_pos() {
        use rust_htslib::bam::{Header, HeaderView, Record};
        
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);
        
        let sam_line = b"read1\t0\t*\t101\t60\t10M2D5M3I4N\t*\t0\t0\tACGTACGTACGTACGTAC\t*";
        let record = Record::from_sam(&header_view, sam_line).unwrap();
        let start_pos = record.pos();
        let cigar_view = record.cigar();
        
        let end_pos = calculate_end_pos(start_pos, &cigar_view);

        // Expected end pos = 100 (0-based) + 10 + 2 + 5 + 4 = 121
        assert_eq!(end_pos, 121);
    }


    
    #[test]
    fn test_adjust_methylation_base_no_methylation() {
        let ref_seq = b"ACGT";
        
        // Without methylation mode, bases should not change
        assert_eq!(
            adjust_methylation_base('T', 'C', 1, false, false, ref_seq, 0),
            'T'
        );
    }
    
    #[test]
    fn test_adjust_methylation_base_cpg_context() {
        let ref_seq = b"CGTA"; // CG dinucleotide at positions 0-1
        
        // CpG context: C in ref, T in read (read1) -> should convert to C
        assert_eq!(
            adjust_methylation_base('T', 'C', 1, true, true, ref_seq, 0),
            'C'
        );
        
        // Non-CpG context: C in ref, T in read (read1) -> should NOT convert
        let ref_seq_non_cpg = b"CATA";
        assert_eq!(
            adjust_methylation_base('T', 'C', 1, true, true, ref_seq_non_cpg, 0),
            'T'
        );
    }
    
    #[test]
    fn test_mismatch_key_equality() {
        let key1 = MismatchKey {
            mismatch_type: "A>G".to_string(),
            read_position: 10,
            read_num: 1,
        };
        
        let key2 = MismatchKey {
            mismatch_type: "A>G".to_string(),
            read_position: 10,
            read_num: 1,
        };
        
        let key3 = MismatchKey {
            mismatch_type: "A>G".to_string(),
            read_position: 10,
            read_num: 2, // Different read number
        };
        
        assert_eq!(key1, key2);
        assert_ne!(key1, key3);
    }
}
