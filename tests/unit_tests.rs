// Unit tests for tasmanian-mismatch

use rustmanian_mismatch::*;

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{Header, HeaderView, Record, Read, Reader, Writer, Format};
    use std::collections::HashMap;

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
    fn test_compute_read_len_mode_from_sample_bam() {
        // Create a temporary BAM file with known read lengths
        let bam_path = "test_sample.bam";
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);
        
        {
            let mut bam_writer = Writer::from_path(
                bam_path,
                &header,
                Format::Bam,
            ).unwrap();
            
            let sam_lines: Vec<&[u8]> = vec![
                b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
                b"read2\t0\t*\t1\t60\t15M\t*\t0\t0\tACGTACGTACGTACG\tIIIIIIIIIIIIIII",
                b"read3\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
                b"read4\t0\t*\t1\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII",
                b"read5\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
            ];
            
            for sam_line in sam_lines {
                let record = Record::from_sam(&header_view, sam_line).unwrap();
                bam_writer.write(&record).unwrap();
            }
        }
        
        // Compute mode read length from sample BAM
        let mode_len = compute_read_len_mode_from_sample_bam(bam_path, 10);
        
        // Clean up temporary BAM file
        std::fs::remove_file(bam_path).unwrap();
        
        // The mode read length should be 10 (appears 3 times)
        assert_eq!(mode_len, 10);
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
    
    #[test]
    fn test_compare_and_count() {
        
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);
        
        // Create a test record
        let sam_line = b"read1\t0\t*\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let record = Record::from_sam(&header_view, sam_line).unwrap();
        
        let ref_seq = b"AGGAACGTAC"; // Has mismatch at position 2: read has G, ref has G (actually match at 2, mismatch at 1)
        let mut local_counts = HashMap::new();
        
        // Test compare_and_count
        compare_and_count(
            &record.seq(),
            record.qual(),
            ref_seq,
            0,      // read position
            0,      // genome position
            false,  // not reverse
            1,      // read 1
            20,     // min base quality
            false,  // not methylation
            false,  // not cpg_only
            &mut local_counts,
            0,      // mode_len
        );
        
        // Should have one entry in counts
        assert_eq!(local_counts.len(), 1);
    }
    
    #[test]
    fn test_create_mismatch_key() {
        let ref_seq = b"GCGTACGTAC";
        
        // Test creating a mismatch key for position 0 (forward strand, no methylation)
        let key = create_mismatch_key(
            'A',         // read base
            'G',         // ref base
            0,           // read position
            10,          // seq length
            false,       // not reverse
            1,           // read 1
            false,       // not methylation mode
            false,       // not cpg_only
            ref_seq,
            0,           // genome position
            0,           // mode_len
        );
        
        assert_eq!(key.mismatch_type, "G>A");
        assert_eq!(key.read_position, 0);
        assert_eq!(key.read_num, 1);
    }
    
    #[test]
    fn test_get_overlap_region() {
        
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);
        
        // Create two overlapping reads
        // Read1: pos 100 (1-based=101), 20M (covers 100-119 in 0-based)
        let sam_line1 = b"read1\t99\t*\t101\t60\t20M\t*\t111\t30\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII";
        let read1 = Record::from_sam(&header_view, sam_line1).unwrap();
        
        // Read2: pos 110 (1-based=111), 20M (covers 110-129 in 0-based)
        // This overlaps with read1 from position 110 to 119 (10 bases)
        let sam_line2 = b"read1\t147\t*\t111\t60\t20M\t*\t101\t-30\tTGCATGCATGCATGCATGCA\tIIIIIIIIIIIIIIIIIIII";
        let read2 = Record::from_sam(&header_view, sam_line2).unwrap();
        
        // Create ReadInfo structs
        let info1 = ReadInfo {
            tid: read1.tid(),
            pos: read1.pos(),
        };
        let info2 = ReadInfo {
            tid: read2.tid(),
            pos: read2.pos(),
        };
        
        let overlap = get_overlap_region(&info1, &info2, &read1, &read2);
        
        // Should have overlap from 110 to 120 (exclusive end)
        assert!(overlap.is_some(), "Expected reads to overlap");
        let (start, end) = overlap.unwrap();
        assert_eq!(start, 110);
        assert_eq!(end, 120);  // Exclusive end position
    }
    
    #[test]
    fn test_process_overlap_region() {
        
        // Create header
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);
        
        // Create overlapping paired reads with a mismatch
        // Read1: forward, pos 100
        let sam_line1 = b"read1\t99\t*\t101\t60\t10M\t*\t111\t20\tACGTACGTAC\tIIIIIIIIII\tMD:Z:10";
        let read1 = Record::from_sam(&header_view, sam_line1).unwrap();
        
        // Read2: reverse, pos 110, overlaps last base of read1
        let sam_line2 = b"read1\t147\t*\t111\t60\t10M\t*\t101\t-20\tGTACGTACGT\tIIIIIIIIII\tMD:Z:10";
        let read2 = Record::from_sam(&header_view, sam_line2).unwrap();
        
        // Create reference genome (tid -1 for unmapped)
        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"ACGTACGTACGTACGTACGT".to_vec());
        
        // Create tid_to_name mapping
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());
        
        let mut local_counts = HashMap::new();
        let mut inconsistency_counts = HashMap::new();
        
        process_overlap_region(
            &read1,
            &read2,
            109,  // overlap start
            109,  // overlap end
            &mut local_counts,
            &mut inconsistency_counts,
            &reference_genome,
            &tid_to_name,
            20,   // min base quality
            false, // not methylation
            false, // not cpg_only
            0,    // mode_len
            0,    // min_map_quality
        );
        
        // Should have processed the overlap region
        // The exact counts depend on the sequence alignment
    }
    
    #[test]
    fn test_process_record_basic() {
        
        // Create header
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);
        
        // Create a simple record - just verify it doesn't panic
        let sam_line = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tMD:Z:5A4";
        let record = Record::from_sam(&header_view, sam_line).unwrap();
        
        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"ACGTAAGCAC".to_vec());
        
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());
        
        let mut local_counts = HashMap::new();
        
        // Just verify the function runs without panicking
        process_record(
            &record,
            &mut local_counts,
            &reference_genome,
            &tid_to_name,
            0.0,   // softclip_threshold
            20,    // min base quality
            false, // not methylation
            false, // not cpg_only
            150,   // read len mode
            0,     // min map quality
        );
        
        // Function completed successfully (unmapped reads may result in no counts)
        assert!(true);
    }
}
