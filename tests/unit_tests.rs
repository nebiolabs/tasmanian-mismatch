// Unit tests for tasmanian-mismatch

use rustmanian_mismatch::*;

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::{Format, Header, HeaderView, Record, Writer};
    use std::collections::HashMap;
    use std::io::Cursor;

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
        assert_eq!(correct_read_len_with_mode(10, 100, 0, false, 1), 10);

        // With mode correction for second half of read
        assert_eq!(correct_read_len_with_mode(60, 100, 150, false, 1), 110);

        // No correction for first half even with mode
        assert_eq!(correct_read_len_with_mode(40, 100, 150, false, 1), 40);

        // Insert mode, read 2: (2*100+10) - (90-30) = 150
        assert_eq!(correct_read_len_with_mode(30, 90, 100, true, 2), 150);
    }

    #[test]
    fn test_calculate_end_pos() {
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        let sam_line: &[u8] = b"read1\t0\t*\t101\t60\t10M2D5M3I4N\t*\t0\t0\tACGTACGTACGTACGTAC\t*";
        let record = Record::from_sam(&header_view, sam_line).unwrap();
        let start_pos = record.pos();
        let cigar_view = record.cigar();

        let end_pos = calculate_end_pos(start_pos, &cigar_view);

        // Expected end pos = 100 (0-based) + 10 + 2 + 5 + 4 = 121
        assert_eq!(end_pos, 121);
    }

    #[test]
    fn test_compute_read_len_max_from_sample_bam() {
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
        let mode_len = compute_read_len_max_from_sample_bam(bam_path, 10);

        // Clean up temporary BAM file
        std::fs::remove_file(bam_path).unwrap();

        // The max read length should be 20
        assert_eq!(mode_len, 20);
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
        let sam_line: &[u8] = b"read1\t0\t*\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let record = Record::from_sam(&header_view, sam_line).unwrap();

        let ref_seq: &[u8] = b"AGGAACGTAC"; // Has mismatch at position 2: read has G, ref has G (actually match at 2, mismatch at 1)
        let mut local_counts = HashMap::new();
        let mut genomic_region_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> =
            HashMap::new();

        // Test compare_and_count
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 20,
            is_methylation: false,
            cpg_only: false,
            mode_len: 0,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };
        let seq = record.seq();
        let read_ctx = ReadContext {
            seq: &seq,
            qual: record.qual(),
            ref_seq,
            is_reverse: false,
            read_num: 1,
        };
        compare_and_count(
            &read_ctx,
            0, // read position
            0, // genome position
            &config,
            &mut local_counts,
            Some(&mut genomic_region_counts),
            "chr1",
        );

        // Should have one entry in counts
        assert_eq!(local_counts.len(), 1);
    }

    #[test]
    fn test_create_mismatch_key() {
        let ref_seq: &[u8] = b"GCGTACGTAC";

        let header = Header::new();
        let header_view = HeaderView::from_header(&header);
        let sam_line: &[u8] = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let record = Record::from_sam(&header_view, sam_line).unwrap();

        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 0,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };
        let seq = record.seq();
        let read_ctx = ReadContext {
            seq: &seq,
            qual: record.qual(),
            ref_seq,
            is_reverse: false,
            read_num: 1,
        };
        // Test creating a mismatch key for position 0 (forward strand, no methylation)
        let key = create_mismatch_key(
            'A', // read base
            'G', // ref base
            0,   // read position
            0,   // genome position
            &read_ctx, &config,
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
        let sam_line1: &[u8] =
            b"read1\t99\t*\t101\t60\t20M\t*\t111\t30\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII";
        let read1 = Record::from_sam(&header_view, sam_line1).unwrap();

        // Read2: pos 110 (1-based=111), 20M (covers 110-129 in 0-based)
        // This overlaps with read1 from position 110 to 119 (10 bases)
        let sam_line2: &[u8] =
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
    fn test_normalize_mismatch_counts_grouped_by_ref_base() {
        let mut counts: HashMap<InsertKey, usize> = HashMap::new();

        // Group 1: read 1, reference order 1, position 10, ref base C
        counts.insert(
            InsertKey {
                base_change: "C>A".to_string(),
                read_num: 1,
                base_position: 10,
                reference_order: 1,
            },
            2,
        );
        counts.insert(
            InsertKey {
                base_change: "C>T".to_string(),
                read_num: 1,
                base_position: 10,
                reference_order: 1,
            },
            3,
        );
        counts.insert(
            InsertKey {
                base_change: "C>C".to_string(),
                read_num: 1,
                base_position: 10,
                reference_order: 1,
            },
            5,
        );

        // Group 2: same location metadata but ref base T (separate denominator)
        counts.insert(
            InsertKey {
                base_change: "T>A".to_string(),
                read_num: 1,
                base_position: 10,
                reference_order: 1,
            },
            4,
        );
        counts.insert(
            InsertKey {
                base_change: "T>T".to_string(),
                read_num: 1,
                base_position: 10,
                reference_order: 1,
            },
            6,
        );

        let normalized = normalize_mismatch_counts(&counts);

        let c_to_t = InsertKey {
            base_change: "C>T".to_string(),
            read_num: 1,
            base_position: 10,
            reference_order: 1,
        };
        let c_to_a = InsertKey {
            base_change: "C>A".to_string(),
            read_num: 1,
            base_position: 10,
            reference_order: 1,
        };
        let t_to_a = InsertKey {
            base_change: "T>A".to_string(),
            read_num: 1,
            base_position: 10,
            reference_order: 1,
        };

        // C group denominator: (2 + 3 + 5) + 1.0 = 11.0
        // T group denominator: (4 + 6) + 1.0 = 11.0
        assert!((normalized[&c_to_t] - (3.0 / 11.0)).abs() < 1e-9);
        assert!((normalized[&c_to_a] - (2.0 / 11.0)).abs() < 1e-9);
        assert!((normalized[&t_to_a] - (4.0 / 11.0)).abs() < 1e-9);
    }

    #[test]
    fn test_normalize_mismatch_counts_separates_read_number() {
        let mut counts: HashMap<InsertKey, usize> = HashMap::new();

        counts.insert(
            InsertKey {
                base_change: "C>T".to_string(),
                read_num: 1,
                base_position: 7,
                reference_order: 1,
            },
            4,
        );
        counts.insert(
            InsertKey {
                base_change: "C>C".to_string(),
                read_num: 1,
                base_position: 7,
                reference_order: 1,
            },
            6,
        );

        counts.insert(
            InsertKey {
                base_change: "C>T".to_string(),
                read_num: 2,
                base_position: 7,
                reference_order: 1,
            },
            1,
        );
        counts.insert(
            InsertKey {
                base_change: "C>C".to_string(),
                read_num: 2,
                base_position: 7,
                reference_order: 1,
            },
            3,
        );

        let normalized = normalize_mismatch_counts(&counts);

        let r1_c_to_t = InsertKey {
            base_change: "C>T".to_string(),
            read_num: 1,
            base_position: 7,
            reference_order: 1,
        };
        let r2_c_to_t = InsertKey {
            base_change: "C>T".to_string(),
            read_num: 2,
            base_position: 7,
            reference_order: 1,
        };

        // R1 C group denominator: (4 + 6) + 1.0 = 11.0
        // R2 C group denominator: (1 + 3) + 1.0 = 5.0
        assert!((normalized[&r1_c_to_t] - (4.0 / 11.0)).abs() < 1e-9);
        assert!((normalized[&r2_c_to_t] - (1.0 / 5.0)).abs() < 1e-9);
    }

    #[test]
    fn test_process_overlap_region() {
        // Use a named contig so records aren't treated as unmapped.
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // Read1: forward, positions 1-10 (0-based 0-9).
        let read1 = Record::from_sam(
            &header_view,
            b"read1\t99\tchr1\t1\t60\t10M\t=\t6\t15\tACGTACGTAC\tIIIIIIIIII",
        )
        .unwrap();
        // Read2: reverse, positions 6-15 (0-based 5-14). Overlaps read1 at 5-9.
        let read2 = Record::from_sam(
            &header_view,
            b"read1\t147\tchr1\t6\t60\t10M\t=\t1\t-15\tGTACGTACGT\tIIIIIIIIII",
        )
        .unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert("chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec());

        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        let mut local_counts = HashMap::new();
        let mut inconsistency_counts = HashMap::new();

        let context = ProcessingContext {
            reference: &reference_genome,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 10,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        // Overlap region [5, 10): 5 bases shared by both reads.
        process_overlap_region(
            &read1,
            &read2,
            (5, 10),
            &mut local_counts,
            &mut inconsistency_counts,
            &context,
            &config,
        );

        // Both reads covered positions 5-9; mismatch counts must be non-empty.
        assert!(local_counts.values().sum::<usize>() > 0);
    }

    #[test]
    fn test_process_record_basic() {
        // Create header
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        // Create a simple record - just verify it doesn't panic
        let sam_line: &[u8] = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tMD:Z:5A4";
        let record = Record::from_sam(&header_view, sam_line).unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"ACGTAAGCAC".to_vec());

        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());

        let mut local_counts = HashMap::new();

        let context = ProcessingContext {
            reference: &reference_genome,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 20,
            is_methylation: false,
            cpg_only: false,
            mode_len: 150,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        // Just verify the function runs without panicking
        process_record(
            &record,
            &mut local_counts,
            None, // no genomic counts
            &context,
            &config,
            None, // no genomic depth
        );

        // Function completed successfully (unmapped reads may result in no counts)
        assert!(true);
    }

    #[test]
    fn test_parse_md_tag_basic() {
        // MD: 5A4 means 5 matches, A mismatch, 4 matches
        let (mismatches, matches) = parse_md_tag("5A4", "test_read");

        assert_eq!(mismatches.len(), 1);
        assert_eq!(mismatches[0], (5, 'A'));

        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0], (0, 5));
        assert_eq!(matches[1], (6, 4));
    }

    #[test]
    fn test_parse_md_tag_with_deletion() {
        // MD: 3^AC5 means 3 matches, AC deleted, 5 matches
        let (mismatches, matches) = parse_md_tag("3^AC5", "test_read");

        assert_eq!(mismatches.len(), 0); // Deletions don't create mismatches
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0], (0, 3));
        assert_eq!(matches[1], (3, 5));
    }

    #[test]
    fn test_parse_md_tag_complex() {
        // MD: 2A3T4 means 2 matches, A mismatch, 3 matches, T mismatch, 4 matches
        let (mismatches, matches) = parse_md_tag("2A3T4", "test_read");

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
        let (mismatches, matches) = parse_md_tag("10", "test_read");

        assert_eq!(mismatches.len(), 0);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0], (0, 10));
    }

    #[test]
    fn test_position_overlaps_intervals() {
        let intervals = vec![
            BedInterval {
                start: 100,
                end: 200,
            },
            BedInterval {
                start: 300,
                end: 400,
            },
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

        let sam_line: &[u8] = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
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
            bed_intervals: &[] as &[BedInterval],
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
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
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
        let sam_line1: &[u8] =
            b"read1\t99\t*\t101\t60\t20M\t*\t111\t30\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII";
        let read1 = Record::from_sam(&header_view, sam_line1).unwrap();

        // Read2: pos 110, 20M (overlaps with read1)
        let sam_line2: &[u8] =
            b"read1\t147\t*\t111\t60\t20M\t*\t101\t-30\tTGCATGCATGCATGCATGCA\tIIIIIIIIIIIIIIIIIIII";
        let read2 = Record::from_sam(&header_view, sam_line2).unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert(
            "*".to_string(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        );

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
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut counts = OverlapCounts {
            local_counts: &mut local_counts,
            overlap_counts: &mut overlap_counts,
            inconsistency_counts: &mut inconsistency_counts,
            genomic_counts: Some(&mut genomic_counts),
            genomic_depth: Some(&mut genomic_depth),
        };

        // Should not panic
        process_paired_reads_with_overlap(
            &read1,
            &read2,
            (110, 120), // overlap (start, end)
            &mut counts,
            &context,
            config,
        );

        assert!(true);
    }

    #[test]
    fn test_rescale_phred_scores() {
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        let sam_line: &[u8] = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let mut record = Record::from_sam(&header_view, sam_line).unwrap();

        let mut reference_genome = HashMap::new();
        reference_genome.insert("*".to_string(), b"ACGTACGTAC".to_vec());

        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(-1i32, "*".to_string());

        let rescaling_matrix = HashMap::new(); // Empty rescaling matrix

        // Should not panic
        rescale_phred_scores(
            &mut record,
            &reference_genome,
            &tid_to_name,
            &rescaling_matrix,
        );

        assert!(true);
    }

    #[test]
    fn test_rescale_phred_scores_with_scaling_factors() {
        let header = Header::new();
        let header_view = HeaderView::from_header(&header);

        let sam_line: &[u8] = b"read1\t0\t*\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
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
        rescale_phred_scores(
            &mut record,
            &reference_genome,
            &tid_to_name,
            &rescaling_matrix,
        );

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

    #[test]
    fn test_write_potential_variants_tsv_writes_header_and_rows() {
        let mut genomic_counts = HashMap::new();
        genomic_counts.insert(
            GenomicMismatchKey {
                chromosome: "chr1".to_string(),
                mismatch_type: "A>G".to_string(),
                genomic_position: 42,
            },
            (3, 10),
        );
        genomic_counts.insert(
            GenomicMismatchKey {
                chromosome: "chr2".to_string(),
                mismatch_type: "INVALID".to_string(),
                genomic_position: 7,
            },
            (1, 2),
        );

        let output_path = "test_potential_variants.tsv";
        write_potential_variants_tsv(&genomic_counts, output_path).unwrap();

        let output = std::fs::read_to_string(output_path).unwrap();
        assert!(
            output.contains("chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth")
        );
        assert!(output.contains("chr1\t42\tA\tG\t3\t10"));
        assert!(!output.contains("chr2\t7"));

        std::fs::remove_file(output_path).unwrap();
    }

    #[test]
    fn test_write_inconsistencies_tsv_sorts_rows() {
        let mut inconsistency_counts = HashMap::new();
        inconsistency_counts.insert(
            InconsistencyKey {
                discordance_type: "R1:C_R2:T".to_string(),
                read1_position: 8,
                read2_position: 5,
            },
            2,
        );
        inconsistency_counts.insert(
            InconsistencyKey {
                discordance_type: "R1:A_R2:G".to_string(),
                read1_position: 3,
                read2_position: 9,
            },
            4,
        );

        let output_path = "test_inconsistencies.tsv";
        write_inconsistencies_tsv(&inconsistency_counts, output_path).unwrap();

        let output = std::fs::read_to_string(output_path).unwrap();
        let lines: Vec<&str> = output.lines().collect();

        assert_eq!(
            lines[0],
            "read1_position\tread2_position\tdiscordance_type\tcount"
        );
        assert_eq!(lines[1], "3\t9\tR1:A_R2:G\t4");
        assert_eq!(lines[2], "8\t5\tR1:C_R2:T\t2");

        std::fs::remove_file(output_path).unwrap();
    }

    #[test]
    fn test_write_mismatch_discounts_tsv_sorts_rows() {
        let mut mismatch_discounts = HashMap::new();
        mismatch_discounts.insert(
            MismatchKey {
                mismatch_type: "G>A".to_string(),
                read_position: 12,
                read_num: 2,
            },
            7,
        );
        mismatch_discounts.insert(
            MismatchKey {
                mismatch_type: "A>G".to_string(),
                read_position: 4,
                read_num: 1,
            },
            3,
        );

        let output_path = "test_mismatch_discounts.tsv";
        write_mismatch_discounts_tsv(&mismatch_discounts, output_path).unwrap();

        let output = std::fs::read_to_string(output_path).unwrap();
        let lines: Vec<&str> = output.lines().collect();

        assert_eq!(
            lines[0],
            "mismatch_type\tread_num\tread_position\tdiscount_count"
        );
        assert_eq!(lines[1], "A>G\t1\t4\t3");
        assert_eq!(lines[2], "G>A\t2\t12\t7");

        std::fs::remove_file(output_path).unwrap();
    }

    #[test]
    fn test_load_discount_table_from_reader_parses_and_accumulates() {
        let input = "base_change\tread_num\tread_position\tdiscount_count\nC>T\t1\t42\t3\nC>T\t1\t42\t2\nA>G\t2\t7\t1\nbad\tline\n";
        let reader = Cursor::new(input.as_bytes());

        let discounts = load_discount_table_from_reader(reader).unwrap();

        let key1 = DiscountKey {
            base_change: "C>T".to_string(),
            read_num: 1,
            base_position: 42,
        };
        let key2 = DiscountKey {
            base_change: "A>G".to_string(),
            read_num: 2,
            base_position: 7,
        };

        assert_eq!(discounts.get(&key1), Some(&5));
        assert_eq!(discounts.get(&key2), Some(&1));
        assert_eq!(discounts.len(), 2);
    }

    #[test]
    fn test_load_rescaling_matrix_from_reader_skips_header_and_invalid_rows() {
        let input = "read_num\tposition\tref_base\tread_base\tscaling_factor\n1\t42\tC\tT\t0.5\nX\tbad\tZ\tY\tnope\n2\t7\tG\tA\t1.25\n";
        let reader = Cursor::new(input.as_bytes());

        let matrix = load_rescaling_matrix_from_reader(reader).unwrap();

        assert_eq!(matrix.get(&(1, 42, 'C', 'T')), Some(&0.5f32));
        assert_eq!(matrix.get(&(2, 7, 'G', 'A')), Some(&1.25f32));
        assert_eq!(matrix.len(), 2);
    }

    #[test]
    fn test_write_rescaling_matrix_output_emits_parseable_rows() {
        let mut counts: HashMap<InsertKey, usize> = HashMap::new();
        counts.insert(
            InsertKey {
                base_change: "A>T".to_string(),
                read_num: 1,
                base_position: 4,
                reference_order: 1,
            },
            3,
        );
        counts.insert(
            InsertKey {
                base_change: "A>A".to_string(),
                read_num: 1,
                base_position: 4,
                reference_order: 1,
            },
            7,
        );

        let output_path = "test_rescaling_matrix.tsv";
        write_rescaling_matrix_output(&counts, Some(output_path)).unwrap();

        let matrix = load_rescaling_matrix(output_path).unwrap();
        assert!(matrix.contains_key(&(1, 4, 'A', 'T')));
        assert!(matrix.contains_key(&(1, 4, 'A', 'A')));
        assert_eq!(matrix.get(&(1, 4, 'A', 'T')), Some(&1.0f32));

        std::fs::remove_file(output_path).unwrap();
    }

    // ── print_position_table / print_main_output ────────────────────────────

    #[test]
    fn test_print_position_table_writes_expected_csv() {
        let mut position_map: HashMap<(u8, usize), HashMap<String, usize>> = HashMap::new();
        position_map
            .entry((1, 3))
            .or_default()
            .insert("A>G".to_string(), 2);
        position_map
            .entry((1, 3))
            .or_default()
            .insert("C>T".to_string(), 1);
        position_map
            .entry((2, 5))
            .or_default()
            .insert("C>T".to_string(), 3);

        let path = "test_print_position_table.csv";
        print_position_table(&position_map, Some(path)).unwrap();

        let contents = std::fs::read_to_string(path).unwrap();
        std::fs::remove_file(path).unwrap();

        // Columns are sorted; rows are sorted by (read, position).
        assert!(contents.starts_with("Read,Position,A>G,C>T\n"));
        assert!(contents.contains("1,3,2,1\n"));
        assert!(contents.contains("2,5,0,3\n"));
    }

    #[test]
    fn test_print_main_output_writes_expected_csv() {
        let mut position_map: HashMap<(u8, usize), HashMap<String, usize>> = HashMap::new();
        position_map
            .entry((1, 0))
            .or_default()
            .insert("G>T".to_string(), 5);

        let path = "test_print_main_output.csv";
        print_main_output(&position_map, Some(path)).unwrap();

        let contents = std::fs::read_to_string(path).unwrap();
        std::fs::remove_file(path).unwrap();

        assert!(contents.contains("G>T"));
        assert!(contents.contains("1,0,5\n"));
    }

    // ── methylation branches ────────────────────────────────────────────────

    #[test]
    fn test_adjust_methylation_base_cpg_context_read2_and_bounds() {
        // Read 2, CpG-only mode: G>A collapses when next base is C/c.
        let ref_seq = b"ACGT";
        // genome_pos=0 ('A'), next=C: not a G ref, no collapse.
        assert_eq!(
            adjust_methylation_base('A', 'A', 2, true, true, ref_seq, 0),
            'A'
        );
        // genome_pos=0, ref G, read A, next base is C → collapse to G.
        let ref_seq2 = b"GCA";
        assert_eq!(
            adjust_methylation_base('A', 'G', 2, true, true, ref_seq2, 0),
            'G'
        );
        // genome_pos=0, ref G, read A, next base is 'c' (lowercase) → collapse to G.
        let ref_seq3 = b"GcA";
        assert_eq!(
            adjust_methylation_base('A', 'G', 2, true, true, ref_seq3, 0),
            'G'
        );
        // genome_pos=0, ref G, read A, next base is 'T' (not C) → keep A.
        let ref_seq4 = b"GTA";
        assert_eq!(
            adjust_methylation_base('A', 'G', 2, true, true, ref_seq4, 0),
            'A'
        );
        // genome_pos at last position (boundary): genome_pos + 1 >= len → keep read_base.
        let ref_seq5 = b"G";
        assert_eq!(
            adjust_methylation_base('A', 'G', 2, true, true, ref_seq5, 0),
            'A'
        );
    }

    #[test]
    fn test_adjust_methylation_base_non_cpg_mode_branches() {
        let ref_seq = b"ACGT";
        // Non-CpG, read 2, G ref, A read → collapse to G.
        assert_eq!(
            adjust_methylation_base('A', 'G', 2, true, false, ref_seq, 0),
            'G'
        );
        // Non-CpG, read 1, C ref, T read → collapse to C.
        assert_eq!(
            adjust_methylation_base('T', 'C', 1, true, false, ref_seq, 0),
            'C'
        );
        // Non-CpG, default fallthrough (read 1, A ref, T read) → passthrough.
        assert_eq!(
            adjust_methylation_base('T', 'A', 1, true, false, ref_seq, 0),
            'T'
        );
    }

    // ── processing helpers ──────────────────────────────────────────────────

    #[test]
    fn test_processing_position_helpers() {
        // insert_mode_read_position: non-stretch, first read.
        assert_eq!(insert_mode_read_position(true, 0, 10, 10, false, None), 1);
        assert_eq!(insert_mode_read_position(true, 9, 10, 10, false, None), 10);
        // second read, non-stretch.
        assert_eq!(insert_mode_read_position(false, 8, 10, 12, false, None), 23);

        // read_mode_read_position: forward, first half.
        assert_eq!(read_mode_read_position(2, 10, false, 10), 3);
        // forward, second half.
        assert_eq!(read_mode_read_position(7, 10, false, 10), 8);
        // reverse.
        assert_eq!(read_mode_read_position(2, 10, true, 10), 8);

        // base_position_for_mode with Read mode.
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 10,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };
        let pos = base_position_for_mode(&config, 2, 10, false, true, false, None);
        assert_eq!(pos, read_mode_read_position(2, 10, false, 10));
    }

    #[test]
    fn test_processing_parsing_and_identity_helpers() {
        // cigar_reference_span
        assert_eq!(cigar_reference_span("10M"), Some(10));
        assert_eq!(cigar_reference_span("5M2D3M"), Some(10));
        assert_eq!(cigar_reference_span("5M3I2M"), Some(7));
        assert_eq!(cigar_reference_span(""), Some(0));
        assert_eq!(cigar_reference_span("10?"), None);

        // softclip_identity
        use rustmanian_mismatch::SoftclipComparison;
        let comps = vec![
            SoftclipComparison {
                read_pos: 0,
                ref_pos: 0,
                read_base: 'A',
                ref_base: 'A',
            },
            SoftclipComparison {
                read_pos: 1,
                ref_pos: 1,
                read_base: 'C',
                ref_base: 'T',
            },
        ];
        assert!((softclip_identity(&comps).unwrap() - 0.5).abs() < 1e-9);
        assert!(softclip_identity(&[]).is_none());
    }

    #[test]
    fn test_build_base_change_modes() {
        // Non-methylation: simple strand-adjusted change.
        assert_eq!(build_base_change(1, 'C', 'T', false, false, None), "C>T");
        assert_eq!(build_base_change(1, 'A', 'G', true, false, None), "T>C");

        // Methylation, non-CpG: C>T on read 1 collapses to C>C.
        assert_eq!(build_base_change(1, 'C', 'T', false, true, None), "C>C");
        // Methylation, G>A on read 2 collapses to G>G.
        assert_eq!(build_base_change(2, 'G', 'A', false, true, None), "G>G");
        // Methylation, no collapse: A>G stays A>G.
        assert_eq!(build_base_change(1, 'A', 'G', false, true, None), "A>G");
        // Methylation mode, CpG context passed as Some('G'): still collapses.
        assert_eq!(
            build_base_change(1, 'C', 'T', false, true, Some('G')),
            "C>C"
        );
        // The current implementation collapses this case as well.
        assert_eq!(
            build_base_change(2, 'G', 'A', false, true, Some('A')),
            "G>G"
        );
    }

    #[test]
    fn test_processing_mate_and_overlap_helpers() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        let paired = Record::from_sam(
            &header_view,
            b"read1\t99\tchr1\t1\t60\t8M\t=\t5\t12\tACGTACGT\tIIIIIIII\tMC:Z:8M",
        )
        .unwrap();
        let no_mc = Record::from_sam(
            &header_view,
            b"read1\t99\tchr1\t1\t60\t8M\t=\t5\t12\tACGTACGT\tIIIIIIII",
        )
        .unwrap();
        let unpaired = Record::from_sam(
            &header_view,
            b"read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tMC:Z:8M",
        )
        .unwrap();

        let mate_end = mc_mate_end(&paired).unwrap();
        assert_eq!(mate_end, 12);
        assert_eq!(overlap_interval(&paired, Some(mate_end)), Some((4, 8)));
        assert_eq!(estimated_fragment_length(&paired, Some(mate_end)), Some(12));

        assert_eq!(mc_mate_end(&no_mc), None);
        assert_eq!(mc_mate_end(&unpaired), None);
        assert_eq!(estimated_fragment_length(&paired, None), None);

        // record_read_num: flag 0x40 = first in template.
        let r1 = Record::from_sam(
            &header_view,
            b"r\t67\tchr1\t1\t60\t5M\t=\t6\t10\tACGTA\tIIIII",
        )
        .unwrap();
        assert_eq!(record_read_num(&r1), 1);

        // read_is_first_in_reference: unpaired returns true.
        assert!(read_is_first_in_reference(&unpaired));
    }

    #[test]
    fn test_compare_record_to_reference_counts_mismatches() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 100);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // 4M: bases ACGT vs ref ACCT → position 2 is G vs C (mismatch).
        let record =
            Record::from_sam(&header_view, b"r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII").unwrap();

        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACCTACGT".to_vec());
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        let context = ProcessingContext {
            reference: &reference,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 4,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut counts: HashMap<InsertKey, usize> = HashMap::new();
        compare_record_to_reference(&record, &context, config, None, &mut counts);

        let total: usize = counts.values().sum();
        assert!(total > 0);

        // There should be a C>G (or strand-equivalent) mismatch key.
        let has_mismatch = counts
            .keys()
            .any(|k| k.base_change.contains('>') && &k.base_change[0..1] != &k.base_change[2..3]);
        assert!(has_mismatch);
    }

    #[test]
    fn test_should_skip_record_branches() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 100);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        let make_record = |flags: u16| {
            Record::from_sam(
                &header_view,
                &format!("r\t{flags}\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII").into_bytes(),
            )
            .unwrap()
        };

        let base_config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 4,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        // No filters → not skipped.
        assert!(!should_skip_record(&make_record(0), base_config));

        // required_flags = 0x1 (paired), record flag = 0 → skipped.
        let mut cfg = base_config;
        cfg.required_flags = 0x1;
        assert!(should_skip_record(&make_record(0), cfg));

        // required_flags satisfied (flag includes 0x1) → not skipped.
        assert!(!should_skip_record(&make_record(0x1), cfg));

        // filter_flags = 0x4 (unmapped). Flag 0x4 set → skipped (is_unmapped() also true).
        let mut cfg2 = base_config;
        cfg2.filter_flags = 0x4;
        assert!(should_skip_record(&make_record(0x4), cfg2));

        // excl_flags: skipped if all bits match.
        let mut cfg3 = base_config;
        cfg3.excl_flags = 0x3;
        assert!(should_skip_record(&make_record(0x3), cfg3));
        assert!(!should_skip_record(&make_record(0x1), cfg3));
    }

    #[test]
    fn test_should_skip_whole_read_for_bed_cursor_advances() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // Record at position 50-58 (8M).
        let record = Record::from_sam(
            &header_view,
            b"r\t0\tchr1\t51\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII",
        )
        .unwrap();

        // BED interval that covers 45-65 → read should be skipped.
        let intervals = vec![
            BedInterval { start: 10, end: 20 },
            BedInterval { start: 45, end: 65 },
        ];
        let mut cursor = 0usize;
        assert!(should_skip_whole_read_for_bed(
            &record,
            true,
            &intervals,
            &mut cursor
        ));
        // Cursor should have advanced past the first interval.
        assert!(cursor >= 1);
    }

    #[test]
    fn test_process_record_with_bed_mask_and_depth() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 100);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // 4M: positions 0-3.
        let record =
            Record::from_sam(&header_view, b"r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII").unwrap();

        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGTACGT".to_vec());
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        // Mask positions 2 and 3 (0-based).
        let bed_intervals = vec![BedInterval { start: 2, end: 4 }];

        let context = ProcessingContext {
            reference: &reference,
            tid_to_name: &tid_to_name,
            bed_intervals: &bed_intervals,
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 4,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut local_counts: HashMap<MismatchKey, usize> = HashMap::new();
        let mut genomic_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> = HashMap::new();
        let mut depth: HashMap<i64, usize> = HashMap::new();

        process_record(
            &record,
            &mut local_counts,
            Some(&mut genomic_counts),
            &context,
            &config,
            Some(&mut depth),
        );

        // Only positions 0 and 1 should be in depth (2 and 3 masked).
        assert_eq!(depth.len(), 2);
        assert!(depth.contains_key(&0));
        assert!(depth.contains_key(&1));
        assert!(!depth.contains_key(&2));

        // Total count entries should be 2 (one per unmasked position).
        let total: usize = local_counts.values().sum();
        assert_eq!(total, 2);
    }

    #[test]
    fn test_process_paired_reads_with_overlap_updates_overlap_and_inconsistency() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 200);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // read1: pos 1-8 (0-based 0-7), read2: pos 5-12 (0-based 4-11). Overlap 4-7.
        let read1 = Record::from_sam(
            &header_view,
            b"r\t99\tchr1\t1\t60\t8M\t=\t5\t12\tACGTACGT\tIIIIIIII\tMC:Z:8M",
        )
        .unwrap();
        let read2 = Record::from_sam(
            &header_view,
            b"r\t147\tchr1\t5\t60\t8M\t=\t1\t-12\tACGTACGT\tIIIIIIII\tMC:Z:8M",
        )
        .unwrap();

        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec());
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        let context = ProcessingContext {
            reference: &reference,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 8,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut local_counts: HashMap<MismatchKey, usize> = HashMap::new();
        let mut overlap_counts: HashMap<MismatchKey, usize> = HashMap::new();
        let mut inconsistency_counts: HashMap<InconsistencyKey, usize> = HashMap::new();
        let mut genomic_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> = HashMap::new();
        let mut depth: HashMap<i64, usize> = HashMap::new();

        let mut counts = OverlapCounts {
            local_counts: &mut local_counts,
            overlap_counts: &mut overlap_counts,
            inconsistency_counts: &mut inconsistency_counts,
            genomic_counts: Some(&mut genomic_counts),
            genomic_depth: Some(&mut depth),
        };

        process_paired_reads_with_overlap(&read1, &read2, (4, 8), &mut counts, &context, config);

        // Overlap positions (4-7) are counted in overlap_counts.
        let overlap_total: usize = overlap_counts.values().sum();
        assert!(overlap_total > 0);
        // Non-overlap positions counted in local_counts.
        let local_total: usize = local_counts.values().sum();
        assert!(local_total > 0);
        // Depth map should have entries.
        assert!(!depth.is_empty());
    }

    #[test]
    fn test_softclip_qualification_and_side_bounds() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        let record = Record::from_sam(
            &header_view,
            b"read1\t0\tchr1\t3\t60\t2S6M2S\t*\t0\t0\tAACCCCGGGG\tIIIIIIIIII",
        )
        .unwrap();

        let ref_seq = b"AATTTTTTAA";
        let qualifying = qualifying_softclip_comparisons(&record, ref_seq, 0.66);

        // Left soft-clip (AA vs AA) qualifies; right soft-clip (GG vs AA) does not.
        assert_eq!(qualifying.len(), 2);
        assert_eq!(qualifying[0].read_pos, 0);
        assert_eq!(qualifying[0].ref_pos, 0);
        assert_eq!(qualifying[1].read_pos, 1);
        assert_eq!(qualifying[1].ref_pos, 1);

        // Out-of-bounds soft-clip extraction returns an empty comparison set.
        let out_of_bounds = softclip_side_comparisons(&record, ref_seq, 9, 4, 8);
        assert!(out_of_bounds.is_empty());
    }

    #[test]
    fn test_mate_end_overlap_and_fragment_helpers() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        let paired = Record::from_sam(
            &header_view,
            b"read1\t99\tchr1\t1\t60\t8M\t=\t5\t12\tACGTACGT\tIIIIIIII\tMC:Z:8M",
        )
        .unwrap();
        let no_mc = Record::from_sam(
            &header_view,
            b"read1\t99\tchr1\t1\t60\t8M\t=\t5\t12\tACGTACGT\tIIIIIIII",
        )
        .unwrap();
        let unpaired = Record::from_sam(
            &header_view,
            b"read1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tMC:Z:8M",
        )
        .unwrap();

        let mate_end = mc_mate_end(&paired).unwrap();
        assert_eq!(mate_end, 12);
        assert_eq!(overlap_interval(&paired, Some(mate_end)), Some((4, 8)));
        assert_eq!(estimated_fragment_length(&paired, Some(mate_end)), Some(12));

        assert_eq!(mc_mate_end(&no_mc), None);
        assert_eq!(mc_mate_end(&unpaired), None);
        assert_eq!(estimated_fragment_length(&paired, None), None);
    }

    #[test]
    fn test_process_record_skips_secondary_records() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // Flag 256 = secondary alignment.
        let secondary = Record::from_sam(
            &header_view,
            b"read1\t256\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
        )
        .unwrap();

        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGTACGT".to_vec());
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        let context = ProcessingContext {
            reference: &reference,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 4,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut local_counts: HashMap<MismatchKey, usize> = HashMap::new();
        let mut genomic_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> = HashMap::new();
        let mut depth: HashMap<i64, usize> = HashMap::new();

        process_record(
            &secondary,
            &mut local_counts,
            Some(&mut genomic_counts),
            &context,
            &config,
            Some(&mut depth),
        );

        assert!(local_counts.is_empty());
        assert!(genomic_counts.is_empty());
        assert!(depth.is_empty());
    }

    #[test]
    fn test_compare_record_cut_vs_stretch_for_read2_overlap() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // read2 (flag 0x93 = 0x80|0x10|0x2|0x1), positions 5-12 (0-based 4-11).
        // Mate is at pos 1-8. Overlap region [4, 8).
        let read2 = Record::from_sam(
            &header_view,
            b"r\t147\tchr1\t5\t60\t8M\t=\t1\t-12\tACGTACGT\tIIIIIIII\tMC:Z:8M",
        )
        .unwrap();

        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec());
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        let context = ProcessingContext {
            reference: &reference,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };

        let mate_end = mc_mate_end(&read2);

        // Cut mode: read2 skips overlap bases → fewer (or no overlap) counts.
        let config_cut = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 8,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut cut_counts: HashMap<InsertKey, usize> = HashMap::new();
        compare_record_to_reference(&read2, &context, config_cut, mate_end, &mut cut_counts);
        let cut_total: usize = cut_counts.values().sum();

        // Stretch mode: no bases are dropped.
        let config_stretch = ProcessingConfig {
            overlap_mode: OverlapMode::Stretch,
            ..config_cut
        };
        let mut stretch_counts: HashMap<InsertKey, usize> = HashMap::new();
        compare_record_to_reference(
            &read2,
            &context,
            config_stretch,
            mate_end,
            &mut stretch_counts,
        );
        let stretch_total: usize = stretch_counts.values().sum();

        // Stretch should have >= Cut (overlap bases not dropped).
        assert!(stretch_total >= cut_total);
    }

    #[test]
    fn test_compare_record_to_reference_cpg_collapse() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 100);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // Sequence: T at position 0 vs reference C at position 0.
        // Reference: CG... so next base is G → CpG context.
        // With cpg_only=true, methylation=true: C>T should collapse to C>C.
        let record =
            Record::from_sam(&header_view, b"r1\t0\tchr1\t1\t60\t2M\t*\t0\t0\tTG\tII").unwrap();

        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"CGTACGT".to_vec());
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        let context = ProcessingContext {
            reference: &reference,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.0,
            min_base_quality: 0,
            is_methylation: true,
            cpg_only: true,
            mode_len: 2,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut counts: HashMap<InsertKey, usize> = HashMap::new();
        compare_record_to_reference(&record, &context, config, None, &mut counts);

        // The C>T at a CpG site should have been collapsed to C>C.
        let has_cpg_collapse = counts.keys().any(|k| k.base_change == "C>C");
        assert!(
            has_cpg_collapse,
            "Expected C>C key for CpG collapse, got: {:?}",
            counts
        );
    }

    #[test]
    fn test_should_skip_whole_read_for_bed_early_returns() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        let record = Record::from_sam(
            &header_view,
            b"read1\t0\tchr1\t21\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII",
        )
        .unwrap();

        let mut cursor = 0usize;
        // filter_whole_reads=false → always false.
        assert!(!should_skip_whole_read_for_bed(
            &record,
            false,
            &[],
            &mut cursor
        ));
        // empty BED → always false.
        assert!(!should_skip_whole_read_for_bed(
            &record,
            true,
            &[],
            &mut cursor
        ));
    }

    // ── count_softclip_mismatches (via process_record) ──────────────────────

    #[test]
    fn test_count_softclip_mismatches_via_process_record() {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 100);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // pos=3 (1-based) → ref_pos=2 (0-based). Left clip 2 → genome 0,1 (AA).
        // 6M at ref 2-7 (CCCCGG). Right clip 2 → genome 8,9 (TT).
        let record = Record::from_sam(
            &header_view,
            b"read1\t0\tchr1\t3\t60\t2S6M2S\t*\t0\t0\tAACCCCGGTT\tIIIIIIIIII",
        )
        .unwrap();

        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AACCCCGGTT".to_vec());
        let mut tid_to_name = HashMap::new();
        tid_to_name.insert(0i32, "chr1".to_string());

        let context = ProcessingContext {
            reference: &reference,
            tid_to_name: &tid_to_name,
            bed_intervals: &[],
        };
        let config = ProcessingConfig {
            softclip_threshold: 0.5,
            min_base_quality: 0,
            is_methylation: false,
            cpg_only: false,
            mode_len: 10,
            min_map_quality: 0,
            required_flags: 0,
            filter_flags: 0,
            excl_flags: 0,
            use_insert_mode: false,
            position_mode: PositionMode::Read,
            overlap_mode: OverlapMode::Cut,
        };

        let mut local_counts: HashMap<MismatchKey, usize> = HashMap::new();
        let mut genomic_counts: HashMap<GenomicMismatchKey, GenomicMismatchValue> = HashMap::new();
        let mut depth: HashMap<i64, usize> = HashMap::new();

        process_record(
            &record,
            &mut local_counts,
            Some(&mut genomic_counts),
            &context,
            &config,
            Some(&mut depth),
        );

        // All bases match → no cross-base substitutions.
        let has_substitution = local_counts.keys().any(|k| {
            k.mismatch_type.len() == 3 && &k.mismatch_type[0..1] != &k.mismatch_type[2..3]
        });
        assert!(
            !has_substitution,
            "unexpected substitution mismatches: {:?}",
            local_counts
        );

        // Right clip CC vs TT → 0% match → below threshold → not added.
        let record2 = Record::from_sam(
            &header_view,
            b"read2\t0\tchr1\t3\t60\t2S6M2S\t*\t0\t0\tAACCCCGGCC\tIIIIIIIIII",
        )
        .unwrap();

        let mut local2: HashMap<MismatchKey, usize> = HashMap::new();
        let mut genomic2: HashMap<GenomicMismatchKey, GenomicMismatchValue> = HashMap::new();
        let mut depth2: HashMap<i64, usize> = HashMap::new();

        process_record(
            &record2,
            &mut local2,
            Some(&mut genomic2),
            &context,
            &config,
            Some(&mut depth2),
        );
        let has_sub2 = local2.keys().any(|k| {
            k.mismatch_type.len() == 3 && &k.mismatch_type[0..1] != &k.mismatch_type[2..3]
        });
        assert!(
            !has_sub2,
            "unexpected substitution mismatches: {:?}",
            local2
        );
    }
}
