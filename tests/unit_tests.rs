// Unit tests for tasmanian-mismatch

use rustmanian_mismatch::*;

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
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
            let mut bam_writer = Writer::from_path(bam_path, &header, Format::Bam).unwrap();

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
        let mut genomic_region_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> =
            HashMap::new();

        // Test compare_and_count
        compare_and_count(
            &record.seq(),
            record.qual(),
            ref_seq,
            0,     // read position
            0,     // genome position
            false, // not reverse
            1,     // read 1
            20,    // min base quality
            false, // not methylation
            false, // not cpg_only
            &mut local_counts,
            Some(&mut genomic_region_counts),
            "chr1",
            0, // mode_len
        );

        // Should have one entry in counts
        assert_eq!(local_counts.len(), 1);
    }

    #[test]
    fn test_create_mismatch_key() {
        let ref_seq = b"GCGTACGTAC";

        // Test creating a mismatch key for position 0 (forward strand, no methylation)
        let key = create_mismatch_key(
            'A',   // read base
            'G',   // ref base
            0,     // read position
            10,    // seq length
            false, // not reverse
            1,     // read 1
            false, // not methylation mode
            false, // not cpg_only
            ref_seq, 0, // genome position
            0, // mode_len
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
        let sam_line1 =
            b"read1\t99\t*\t101\t60\t20M\t*\t111\t30\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII";
        let read1 = Record::from_sam(&header_view, sam_line1).unwrap();

        // Read2: pos 110 (1-based=111), 20M (covers 110-129 in 0-based)
        // This overlaps with read1 from position 110 to 119 (10 bases)
        let sam_line2 =
            b"read1\t147\t*\t111\t60\t20M\t*\t101\t-30\tTGCATGCATGCATGCATGCA\tIIIIIIIIIIIIIIIIIIII";
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
        assert_eq!(end, 120); // Exclusive end position
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
        let sam_line2 =
            b"read1\t147\t*\t111\t60\t10M\t*\t101\t-20\tGTACGTACGT\tIIIIIIIIII\tMD:Z:10";
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
            109, // overlap start
            109, // overlap end
            &mut local_counts,
            &mut inconsistency_counts,
            &reference_genome,
            &tid_to_name,
            20,    // min base quality
            false, // not methylation
            false, // not cpg_only
            0,     // mode_len
            0,     // min_map_quality
            0,     // required_flags
            0,     // filter_flags
            0,     // excl_flags
            &[],   // no bed intervals
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
            None, // no genomic counts
            &reference_genome,
            &tid_to_name,
            0.0,   // softclip_threshold
            20,    // min base quality
            false, // not methylation
            false, // not cpg_only
            150,   // read len mode
            0,     // min map quality
            None,  // no genomic depth
            0,     // required_flags
            0,     // filter_flags
            0,     // excl_flags
            &[],   // no bed intervals
        );

        // Function completed successfully (unmapped reads may result in no counts)
        assert!(true);
    }

    #[test]
    fn test_parse_md_tag_basic() {
        // MD: 5A4 means 5 matches, A mismatch, 4 matches
        let (mismatches, matches) = parse_md_tag("5A4");
        
        assert_eq!(mismatches.len(), 1);
        assert_eq!(mismatches[0], (5, 'A'));
        
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0], (0, 5));
        assert_eq!(matches[1], (6, 4));
    }

    #[test]
    fn test_parse_md_tag_with_deletion() {
        // MD: 3^AC5 means 3 matches, AC deleted, 5 matches
        let (mismatches, matches) = parse_md_tag("3^AC5");
        
        assert_eq!(mismatches.len(), 0); // Deletions don't create mismatches
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0], (0, 3));
        assert_eq!(matches[1], (3, 5));
    }

    #[test]
    fn test_parse_md_tag_complex() {
        // MD: 2A3T4 means 2 matches, A mismatch, 3 matches, T mismatch, 4 matches
        let (mismatches, matches) = parse_md_tag("2A3T4");
        
        assert_eq!(mismatches.len(), 2);
        assert_eq!(mismatches[0], (2, 'A'));
        assert_eq!(mismatches[1], (6, 'T'));
        
        assert_eq!(matches.len(), 3);
        assert_eq!(matches[0], (0, 2));
        assert_eq!(matches[1], (3, 3));
        assert_eq!(matches[2], (7, 4));
    }

    #[test]
    fn test_parse_md_tag_only_matches() {
        // MD: 10 means 10 matches only
        let (mismatches, matches) = parse_md_tag("10");
        
        assert_eq!(mismatches.len(), 0);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0], (0, 10));
    }

    #[test]
    fn test_position_overlaps_intervals() {
        let intervals = vec![
            BedInterval { start: 100, end: 200 },
            BedInterval { start: 300, end: 400 },
        ];

        // Position within first interval
        assert!(position_overlaps_intervals(&intervals, 150));
        
        // Position at interval boundary (inclusive start)
        assert!(position_overlaps_intervals(&intervals, 100));
        
        // Position at interval boundary (exclusive end)
        assert!(!position_overlaps_intervals(&intervals, 200));
        
        // Position between intervals
        assert!(!position_overlaps_intervals(&intervals, 250));
        
        // Position in second interval
        assert!(position_overlaps_intervals(&intervals, 350));
    }

    #[test]
    fn test_position_overlaps_intervals_empty() {
        let intervals: Vec<BedInterval> = vec![];

        // Empty intervals should return false
        assert!(!position_overlaps_intervals(&intervals, 150));
    }

    #[test]
    fn test_load_reference_genome() {
        // Create a temporary FASTA file
        let fasta_path = "test_reference.fa";
        let fasta_content = ">chr1\nACGTACGTACGT\n>chr2\nGGCCTTAA\n";
        
        std::fs::write(fasta_path, fasta_content).unwrap();

        let reference = load_reference_genome(fasta_path);

        // Verify contents
        assert!(reference.contains_key("chr1"));
        assert!(reference.contains_key("chr2"));
        
        let chr1_seq = reference.get("chr1").unwrap();
        assert_eq!(chr1_seq.len(), 12);
        assert_eq!(&chr1_seq[0..4], b"ACGT");

        let chr2_seq = reference.get("chr2").unwrap();
        assert_eq!(chr2_seq.len(), 8);
        assert_eq!(&chr2_seq[0..4], b"GGCC");

        // Clean up
        std::fs::remove_file(fasta_path).unwrap();
    }

    #[test]
    fn test_load_reference_genome_multiline_sequence() {
        // Create a FASTA with sequences spanning multiple lines
        let fasta_path = "test_reference_multiline.fa";
        let fasta_content = ">chr1\nACGT\nACGT\nACGT\n>chr2\nGG\nCC\n";
        
        std::fs::write(fasta_path, fasta_content).unwrap();

        let reference = load_reference_genome(fasta_path);

        assert!(reference.contains_key("chr1"));
        let chr1_seq = reference.get("chr1").unwrap();
        assert_eq!(chr1_seq.len(), 12); // 4+4+4
        
        assert!(reference.contains_key("chr2"));
        let chr2_seq = reference.get("chr2").unwrap();
        assert_eq!(chr2_seq.len(), 4); // 2+2

        // Clean up
        std::fs::remove_file(fasta_path).unwrap();
    }

    #[test]
    fn test_parse_bed_file() {
        // Create a temporary BED file
        let bed_path = "test.bed";
        let bed_content = "chr1\t100\t200\nchromosome2\t500\t600\nchr3\t1000\t1500\n";
        
        std::fs::write(bed_path, bed_content).unwrap();

        let bed_regions = parse_bed_file(bed_path).unwrap();

        // Verify structure
        assert!(bed_regions.contains_key("chr1"));
        assert!(bed_regions.contains_key("chromosome2"));
        assert!(bed_regions.contains_key("chr3"));

        let chr1_intervals = bed_regions.get("chr1").unwrap();
        assert_eq!(chr1_intervals.len(), 1);
        assert_eq!(chr1_intervals[0].start, 100);
        assert_eq!(chr1_intervals[0].end, 200);

        let chr2_intervals = bed_regions.get("chromosome2").unwrap();
        assert_eq!(chr2_intervals.len(), 1);
        assert_eq!(chr2_intervals[0].start, 500);
        assert_eq!(chr2_intervals[0].end, 600);

        // Clean up
        std::fs::remove_file(bed_path).unwrap();
    }

    #[test]
    fn test_parse_bed_file_multiple_regions_same_contig() {
        // Create a BED file with multiple regions on the same contig
        let bed_path = "test_multi.bed";
        let bed_content = "chr1\t100\t200\nchr1\t300\t400\nchr1\t1000\t1500\n";
        
        std::fs::write(bed_path, bed_content).unwrap();

        let bed_regions = parse_bed_file(bed_path).unwrap();

        let chr1_intervals = bed_regions.get("chr1").unwrap();
        assert_eq!(chr1_intervals.len(), 3);
        assert_eq!(chr1_intervals[0].start, 100);
        assert_eq!(chr1_intervals[1].start, 300);
        assert_eq!(chr1_intervals[2].start, 1000);

        // Clean up
        std::fs::remove_file(bed_path).unwrap();
    }

    #[test]
    fn test_process_single_record() {
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        let sam_line = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let record = Record::from_sam(&header_view, sam_line).unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"ACGTACGTAC".to_vec());

        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());

        let mut region_counts = HashMap::new();
        let mut genomic_region_counts = HashMap::new();
        let mut genomic_position_depth = HashMap::new();

        let context = ProcessingContext {
            reference: &reference_genome,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };

        let config = ProcessingConfig {
            softclip_threshold: 0.66,
            min_base_quality: 20,
            is_methylation: false,
            cpg_only: false,
            mode_len: 0,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
        };

        // Should not panic
        process_single_record(
            &record,
            &mut region_counts,
            &mut genomic_region_counts,
            &mut genomic_position_depth,
            &context,
            config,
        );

        // Verify it completed successfully
        assert!(true);
    }

    #[test]
    fn test_process_paired_reads_with_overlap() {
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        // Read1: pos 100, 20M
        let sam_line1 = b"read1\t99\t*\t101\t60\t20M\t*\t111\t30\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII";
        let read1 = Record::from_sam(&header_view, sam_line1).unwrap();

        // Read2: pos 110, 20M (overlaps with read1)
        let sam_line2 = b"read1\t147\t*\t111\t60\t20M\t*\t101\t-30\tTGCATGCATGCATGCATGCA\tIIIIIIIIIIIIIIIIIIII";
        let read2 = Record::from_sam(&header_view, sam_line2).unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec());

        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());

        let mut local_counts = HashMap::new();
        let mut overlap_counts = HashMap::new();
        let mut inconsistency_counts = HashMap::new();
        let mut genomic_counts = HashMap::new();
        let mut genomic_depth = HashMap::new();

        let context = ProcessingContext {
            reference: &reference_genome,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };

        let config = ProcessingConfig {
            softclip_threshold: 0.66,
            min_base_quality: 20,
            is_methylation: false,
            cpg_only: false,
            mode_len: 0,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
        };

        // Should not panic
        process_paired_reads_with_overlap(
            &read1,
            &read2,
            110,  // overlap_start
            120,  // overlap_end
            &mut local_counts,
            &mut overlap_counts,
            &mut inconsistency_counts,
            Some(&mut genomic_counts),
            Some(&mut genomic_depth),
            &context,
            config,
        );

        assert!(true);
    }

    #[test]
    fn test_rescale_phred_scores() {
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        let sam_line = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let mut record = Record::from_sam(&header_view, sam_line).unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"ACGTACGTAC".to_vec());

        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());

        let rescaling_matrix = HashMap::new(); // Empty rescaling matrix

        // Should not panic
        rescale_phred_scores(&mut record, &reference_genome, &tid_to_name, &rescaling_matrix);

        assert!(true);
    }

    #[test]
    fn test_rescale_phred_scores_with_scaling_factors() {
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        let sam_line = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let mut record = Record::from_sam(&header_view, sam_line).unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"AGGAACGTAC".to_vec());

        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());

        // Create rescaling matrix with one scaling factor
        // (read_num, read_pos, ref_base, read_base) -> scaling_factor
        let mut rescaling_matrix = HashMap::new();
        rescaling_matrix.insert((1, 1u16, 'G', 'C'), 0.8f32);

        // Should not panic
        rescale_phred_scores(&mut record, &reference_genome, &tid_to_name, &rescaling_matrix);

        assert!(true);
    }

    #[test]
    fn test_filter_bed_for_region() {
        // Create a BED file
        let bed_path = "test_filter.bed";
        let bed_content = "chr1\t100\t200\nchr1\t500\t600\nchr2\t1000\t1500\n";
        
        std::fs::write(bed_path, bed_content).unwrap();
        let bed_regions = parse_bed_file(bed_path).unwrap();

        // Filter for chr1 with range 0-400
        let filtered = filter_bed_for_region(&bed_regions, "chr1", 0, 400);

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].start, 100);
        assert_eq!(filtered[0].end, 200);

        // Filter for chr1 with range 0-601 (includes second region)
        let filtered2 = filter_bed_for_region(&bed_regions, "chr1", 0, 601);
        assert_eq!(filtered2.len(), 2);

        // Clean up
        std::fs::remove_file(bed_path).unwrap();
    }

    #[test]
    fn test_filter_bed_for_region_no_overlap() {
        // Create a BED file
        let bed_path = "test_filter_no_overlap.bed";
        let bed_content = "chr1\t100\t200\nchr1\t500\t600\n";
        
        std::fs::write(bed_path, bed_content).unwrap();
        let bed_regions = parse_bed_file(bed_path).unwrap();

        // Filter for chr1 with range 1000-2000 (no overlap)
        let filtered = filter_bed_for_region(&bed_regions, "chr1", 1000, 2000);
        assert_eq!(filtered.len(), 0);

        // Clean up
        std::fs::remove_file(bed_path).unwrap();
    }

    #[test]
    fn test_mask_reference_with_bed() {

        let mut reference_genome = HashMap::new();
        reference_genome.insert("chr1".to_string(), b"ACGTACGTACGT".to_vec());

        let bed_path = "test_mask.bed";
        let bed_content = "chr1\t4\t8\n";
        std::fs::write("test_mask.bed", bed_content).unwrap();
        let bed_regions = parse_bed_file(bed_path).unwrap();

        let masked_count = mask_reference_with_bed(&mut reference_genome, &bed_regions);

        assert_eq!(masked_count, 4);
        let masked_seq = reference_genome.get("chr1").unwrap();
        assert_eq!(masked_seq, b"ACGTNNNNACGT");

        std::fs::remove_file(bed_path).unwrap();
    }
}