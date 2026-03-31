%%{init: {'flowchart': {'htmlLabels': true}, 'themeVariables': {'fontSize': '30px', 'primaryTextColor': '#111827'}}}%%
flowchart LR
    A["Input Files<br>BAM + indexed BAM, reference FASTA, optional BED"] --> B["Parse CLI Options<br>threads, MAPQ, base quality, methylation, BED mode, thresholds"]
    B --> C{"use_read_len_max?"}
    C -- yes --> D["Sample first 100k BAM records<br>and compute modal read length"]
    C -- no --> E["Load reference FASTA into memory"]
    D --> E
    E --> F{"BED file provided?"}
    F -- yes --> G["Parse BED into per-chromosome intervals"]
    F -- no --> K["Open BAM header and build contigs map"]
    G --> H{"BED mode"}
    H -- mask --> I["Mask reference bases to N<br>for BED intervals"]
    H -- filter --> J["Keep BED intervals for whole-read filtering"]
    I --> K
    J --> K
    K --> L["Split each contig into regions<br>default 10 Mb chunks"]
    L --> M["Process regions in parallel<br>IndexedReader fetch per region"]
    M --> O["Prefilter BED intervals for current chunk"]
    O --> N["Keep only reads that start in this region<br>to avoid duplicate counting"]
    N --> P{"filter mode and read overlaps BED?" <div style="width:300px">}
    P -- yes --> Q["Skip whole read"]
    P -- no --> R{"paired read with mapped mate?" <div style="width:250px">}
    Q --> N
    R -- no --> T["Process single record immediately"]
    R -- yes --> S["Group reads by qname inside chunk"]
    T --> Y["Walk CIGAR operations"]
    S --> U{"two mates found in chunk?"  <div style="width:250px">}
    U -- no --> X["Non-overlap path<br>process each read separately"]
    U -- yes --> V{"mates overlap genomically?"  <div style="width:300px">}
    V -- yes --> W["Single-pass overlap path<br>count overlap mismatches and read-pair inconsistencies"]
    V -- no --> X
    W --> AG["Apply genomic thresholds<br>keep only mismatch positions meeting count and depth cutoffs"]
    X --> Y
    Y --> Z{"record fails filters?<br>unmapped, MAPQ, flags, secondary, supplementary"  <div style="width:400px">}
    Z -- yes --> AA["Skip record"]
    Z -- no --> AB["For aligned bases<br>skip BED-masked positions"]
    AA --> AG
    AB --> AC["Apply strand normalization<br>and optional methylation adjustment"]
    AC --> AD["Count match or mismatch by read number and read position"]
    AD --> AE["Track genomic depth and genomic mismatch keys"]
    AE --> AF["Handle soft clips<br>only count if enough clip bases match reference"]
    AF --> AG
    AG --> AH["Merge local region counts into global maps"]
    AH --> AI["After all regions finish"]
    AI --> AJ["Write potential_variants.tsv"]
    AJ --> AK["Print main CSV to stdout<br>Read,Position,mismatch counts"]
    AK --> AL{"overlapping pairs found?"}
    AL -- yes --> AM["Print overlap mismatch table"]
    AL -- no --> AO["Analysis complete"]
    AM --> AN["Print read-pair inconsistency table"]
    AN --> AO

    style A fill:#f94144,stroke:#7f1d1d,stroke-width:2px,color:#ffffff
    style B fill:#f3722c,stroke:#7c2d12,stroke-width:2px,color:#ffffff
    style C fill:#f8961e,stroke:#78350f,stroke-width:2px,color:#111827
    style D fill:#f9844a,stroke:#7c2d12,stroke-width:2px,color:#ffffff
    style E fill:#f9c74f,stroke:#854d0e,stroke-width:2px,color:#111827
    style F fill:#90be6d,stroke:#365314,stroke-width:2px,color:#111827
    style G fill:#43aa8b,stroke:#14532d,stroke-width:2px,color:#ffffff
    style H fill:#4d908e,stroke:#134e4a,stroke-width:2px,color:#ffffff
    style I fill:#577590,stroke:#1e3a8a,stroke-width:2px,color:#ffffff
    style J fill:#277da1,stroke:#082f49,stroke-width:2px,color:#ffffff
    style K fill:#9b5de5,stroke:#581c87,stroke-width:2px,color:#ffffff
    style L fill:#7209b7,stroke:#4c1d95,stroke-width:2px,color:#ffffff
    style M fill:#560bad,stroke:#3b0764,stroke-width:2px,color:#ffffff
    style N fill:#480ca8,stroke:#312e81,stroke-width:2px,color:#ffffff
    style O fill:#3a0ca3,stroke:#312e81,stroke-width:2px,color:#ffffff
    style P fill:#3f37c9,stroke:#1e3a8a,stroke-width:2px,color:#ffffff
    style Q fill:#4361ee,stroke:#1d4ed8,stroke-width:2px,color:#ffffff
    style R fill:#4895ef,stroke:#075985,stroke-width:2px,color:#ffffff
    style S fill:#4cc9f0,stroke:#155e75,stroke-width:2px,color:#111827
    style T fill:#06d6a0,stroke:#14532d,stroke-width:2px,color:#111827
    style U fill:#ffd166,stroke:#854d0e,stroke-width:2px,color:#111827
    style V fill:#ef476f,stroke:#881337,stroke-width:2px,color:#ffffff
    style W fill:#118ab2,stroke:#082f49,stroke-width:2px,color:#ffffff
    style X fill:#073b4c,stroke:#082f49,stroke-width:2px,color:#ffffff
    style Y fill:#8ecae6,stroke:#075985,stroke-width:2px,color:#111827
    style Z fill:#ffb703,stroke:#92400e,stroke-width:2px,color:#111827
    style AA fill:#fb8500,stroke:#9a3412,stroke-width:2px,color:#ffffff
    style AB fill:#bde0fe,stroke:#1d4ed8,stroke-width:2px,color:#111827
    style AC fill:#a2d2ff,stroke:#1e3a8a,stroke-width:2px,color:#111827
    style AD fill:#cdb4db,stroke:#6b21a8,stroke-width:2px,color:#111827
    style AE fill:#ffc8dd,stroke:#9d174d,stroke-width:2px,color:#111827
    style AF fill:#ffafcc,stroke:#9d174d,stroke-width:2px,color:#111827
    style AG fill:#d4a373,stroke:#78350f,stroke-width:2px,color:#111827
    style AH fill:#ccd5ae,stroke:#365314,stroke-width:2px,color:#111827
    style AI fill:#e9edc9,stroke:#4d7c0f,stroke-width:2px,color:#111827
    style AJ fill:#faedcd,stroke:#92400e,stroke-width:2px,color:#111827
    style AK fill:#fefae0,stroke:#a16207,stroke-width:2px,color:#111827
    style AL fill:#d8f3dc,stroke:#166534,stroke-width:2px,color:#111827
    style AM fill:#95d5b2,stroke:#166534,stroke-width:2px,color:#111827
    style AN fill:#52b788,stroke:#166534,stroke-width:2px,color:#ffffff
    style AO fill:#2d6a4f,stroke:#14532d,stroke-width:2px,color:#ffffff