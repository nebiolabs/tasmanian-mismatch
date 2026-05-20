# Tasmanian Mismatch — Control Flow

```mermaid
flowchart TD
    A([main]) --> B[Args::parse]
    B --> C[configure_thread_pool\nrayon global pool]
    C --> D[compute_read_len_max_from_sample_bam\nsample first 10 000 records → max_read_len]
    D --> E[load_reference_genome\nFASTA → HashMap chr → Vec u8]

    E --> F{bed_file\nprovided?}
    F -- No --> G[bed_for_filtering = None]
    F -- Yes → parse_bed_file --> H{bed_filter_mode}
    H -- mask --> I[mask_reference_with_bed\nN-mask ref sequence in place\nbed_for_filtering = None]
    H -- filter --> J[keep BedRegions as Arc\nfor whole-read overlap filter\nbed_for_filtering = Some]

    G & I & J --> K[build_tid_map_and_regions\ntid→name map + chunked regions\ndefault 1 Mbp chunks]

    K --> L[par_iter regions\nrayon parallel across chunks]

    subgraph process_region [process_region — one thread per chunk]
        L --> PR1[IndexedReader::fetch tid start end]
        PR1 --> PR2{for each BAM record}
        PR2 -- read_start outside chunk bounds --> PR2
        PR2 --> PR3[should_skip_whole_read_for_bed?\nBED cursor scan]
        PR3 -- overlaps BED & filter mode --> PR2
        PR3 -- passes --> PR4[should_skip_record?\nunmapped / mapq / flag filters]
        PR4 -- skip --> PR2
        PR4 -- passes --> PR5[mc_mate_end\nparse MC tag → mate end pos]
        PR5 --> PR6[compare_record_to_reference]
        PR6 --> PR2
    end

    subgraph compare_record_to_reference [compare_record_to_reference — per record]
        PR6 --> CR1[resolve chr_name & ref_seq]
        CR1 --> CR2[read metadata\nread_num 1 or 2\nreference_order First/Second by coord]
        CR2 --> CR3[overlap_interval\nfrom mpos + mc_mate_end]
        CR3 --> CR4[qualifying_softclip_comparisons\nidentity threshold 0.66]
        CR4 --> CR5{CIGAR walk\nMatch Equal Diff}

        CR5 -- Ins / SoftClip --> ADV_R[advance read_pos]
        ADV_R --> CR5
        CR5 -- Del / RefSkip --> ADV_G[advance ref_pos]
        ADV_G --> CR5

        CR5 -- aligned base --> OV{in overlap &\nCut mode & R2?}
        OV -- yes → skip base --> CR5
        OV -- no --> BC[build_base_change\nstrand adjust ref & read bases]

        BC --> METH{methylation_mode?}
        METH -- No --> BC2[ref>read as-is]
        METH -- Yes → R1 C→T or R2 G→A --> METH2{cpg_only?}
        METH2 -- No → collapse all --> BC3[C>C or G>G]
        METH2 -- Yes → check next ref base --> BC4{CpG context?}
        BC4 -- Yes --> BC3
        BC4 -- No --> BC2

        BC2 & BC3 --> POS[base_position_for_mode]
        POS --> POSM{position_mode}
        POSM -- Read --> RMP[read_mode_read_position\npos within read corrected for max_read_len]
        POSM -- Insert --> IMP[insert_mode_read_position\nR1: distance from fragment start\nR2: distance from fragment end\noptional cubic stretch for overlap]
        RMP & IMP --> KEY[InsertKey\nbase_change read_num\nbase_position reference_order\n+= 1]
        KEY --> CR5

        CR5 -- CIGAR done --> SC[process qualifying soft-clips\nsame build_base_change + pos pipeline]
        SC --> KEY
    end

    PR6 --> MERGE[merge local_counts into\nglobal_counts HashMap InsertKey usize\nMutex lock per chunk]

    MERGE --> DT{discount_table\nprovided?}
    DT -- Yes --> DTA[load_discount_table\napply_external_discounts\nscale counts by position-keyed factors]
    DTA --> OUT
    DT -- No --> OUT

    OUT{output mode}
    OUT -- emit_rescaling_matrix --> OM[write_rescaling_matrix_output\nread_num pos ref read scaling_factor TSV]
    OUT -- normalize --> ON[write_normalized_output\ncounts ÷ total per ref_base group]
    OUT -- raw --> OR[write_output\nraw InsertKey counts TSV]
```

## Key design decisions

| Decision point | Options |
|---|---|
| BED handling | `mask` mutates the reference in place; `filter` skips whole reads at query time |
| Overlap mode | `Cut` drops R2 bases in the overlap; `Stretch` applies a cubic spline to map both reads' positions into `[1, 2×max_read_len]` without dropping bases |
| Position mode | `Read` — position within the read; `Insert` — position relative to fragment start/end |
| Methylation | Global C→T / G→A collapse, or CpG-context-only |
| Parallelism | `rayon` parallel over 1 Mbp chunks; each chunk uses an independent `IndexedReader`; results merged under a `Mutex` |
| Output | Raw counts, normalized frequencies, or rescaling matrix rows for `tasmanian-rescale-quality` |
