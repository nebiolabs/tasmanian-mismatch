flowchart TB
    %% Input stage
    subgraph Input["📥 INPUT"]
        BAM[BAM File]
        FASTA[Reference FASTA]
        BED[BED File<br/>optional]
        ARGS[Command-line Arguments<br/>threads, quality thresholds,<br/>methylation mode, etc.]
    end
    
    %% Initialization stage
    subgraph Init["⚙️ INITIALIZATION"]
        PARSE[Parse Arguments]
        MODELEN{Use Read<br/>Length Mode?}
        CALCMODE[Compute Mode<br/>Read Length<br/>from 100k sample]
        LOADREF[Load Reference<br/>Genome into Memory]
        BEDPROC{BED File<br/>Provided?}
        BEDMODE{BED Filter<br/>Mode?}
        MASKREF[Mask Reference<br/>at BED Regions]
        KEEPBED[Keep BED Regions<br/>for Read Filtering]
    end
    
    %% Parallel processing setup
    subgraph Setup["🔧 PARALLEL SETUP"]
        HEADER[Parse BAM Header<br/>tid→chromosome mapping]
        REGIONS[Split Chromosomes<br/>into Regions<br/>default: 10Mb chunks]
        THREADPOOL[Initialize Thread Pool]
    end
    
    %% Main processing loop
    subgraph ParallelProc["⚡ PARALLEL PROCESSING<br/>one thread per region"]
        FETCH[Fetch Reads<br/>in Region]
        
        subgraph ReadLoop["For Each Read in Region"]
            CHECKSTART{Read Starts<br/>in Region?}
            BEDFILTER{BED Filter<br/>Mode?}
            SKIPREAD{Read Overlaps<br/>BED Interval?}
            
            subgraph RecordProc["🔬 RECORD PROCESSING"]
                CIGAR[Parse CIGAR String]
                ALIGN[Align Read to<br/>Reference Position]
                
                subgraph BaseLoop["For Each Base"]
                    QUALITY{Base Quality<br/>≥ threshold?}
                    BEDMASK{Position in<br/>Masked Region?}
                    GETBASES[Extract Read Base<br/>& Ref Base]
                    STRAND{Reverse<br/>Strand?}
                    COMPLEMENT[Apply Reverse<br/>Complement]
                    METH{Methylation<br/>Mode?}
                    METHADJ[Adjust C/T<br/>for Methylation]
                    COMPARE[Compare Read<br/>to Reference]
                    LOCALMM[Increment Local<br/>Mismatch Count]
                    GENOMIC{Is<br/>Mismatch?}
                    GENCOUNT[Track Genomic<br/>Position Count]
                end
                
                SOFTCLIP[Process Soft-clipped<br/>Bases]
            end
            
            PAIRED{Is Read<br/>Paired?}
            
            subgraph OverlapProc["👥 OVERLAP DETECTION"]
                CACHE{Mate Already<br/>Cached?}
                STORE[Store Read<br/>in Cache]
                OVERLAP{Reads<br/>Overlap?}
                PROCOVER[Process Overlap<br/>Region for Both Reads]
                INCONS[Detect & Count<br/>Inconsistencies]
            end
        end
        
        THRESHOLD[Filter Genomic<br/>Mismatches by Threshold]
        MERGE[Merge Local Counts<br/>to Global<br/>locked mutex]
    end
    
    %% Output stage
    subgraph Output["📤 OUTPUT GENERATION"]
        STRUCT[Restructure Data<br/>read,pos → mismatch types]
        VARFILE[Write<br/>potential_variants.tsv<br/>genomic positions ≥ threshold]
        MAINCSV[Write Main CSV<br/>mismatch counts by<br/>read position]
        OVERCSV{Overlaps<br/>Found?}
        OVEROUT[Write Overlap<br/>Mismatch Table]
        INCONSOUT[Write Read Pair<br/>Inconsistency Table]
    end
    
    %% Flow connections
    BAM --> PARSE
    FASTA --> PARSE
    BED --> PARSE
    ARGS --> PARSE
    
    PARSE --> MODELEN
    MODELEN -->|Yes| CALCMODE
    MODELEN -->|No| LOADREF
    CALCMODE --> LOADREF
    
    LOADREF --> BEDPROC
    BEDPROC -->|Yes| BEDMODE
    BEDPROC -->|No| HEADER
    BEDMODE -->|mask| MASKREF
    BEDMODE -->|filter| KEEPBED
    MASKREF --> HEADER
    KEEPBED --> HEADER
    
    HEADER --> REGIONS
    REGIONS --> THREADPOOL
    THREADPOOL --> FETCH
    
    FETCH --> CHECKSTART
    CHECKSTART -->|No| CHECKSTART
    CHECKSTART -->|Yes| BEDFILTER
    BEDFILTER -->|Yes| SKIPREAD
    BEDFILTER -->|No| CIGAR
    SKIPREAD -->|Yes| CHECKSTART
    SKIPREAD -->|No| CIGAR
    
    CIGAR --> ALIGN
    ALIGN --> QUALITY
    QUALITY -->|No| QUALITY
    QUALITY -->|Yes| BEDMASK
    BEDMASK -->|Yes| QUALITY
    BEDMASK -->|No| GETBASES
    GETBASES --> STRAND
    STRAND -->|Yes| COMPLEMENT
    STRAND -->|No| METH
    COMPLEMENT --> METH
    METH -->|Yes| METHADJ
    METH -->|No| COMPARE
    METHADJ --> COMPARE
    COMPARE --> LOCALMM
    LOCALMM --> GENOMIC
    GENOMIC -->|Yes| GENCOUNT
    GENOMIC -->|No| QUALITY
    GENCOUNT --> QUALITY
    
    ALIGN --> SOFTCLIP
    SOFTCLIP --> PAIRED
    
    PAIRED -->|Yes| CACHE
    PAIRED -->|No| CHECKSTART
    CACHE -->|No| STORE
    CACHE -->|Yes| OVERLAP
    STORE --> CHECKSTART
    OVERLAP -->|No| CHECKSTART
    OVERLAP -->|Yes| PROCOVER
    PROCOVER --> INCONS
    INCONS --> CHECKSTART
    
    CHECKSTART -.Region Complete.-> THRESHOLD
    THRESHOLD --> MERGE
    MERGE -.All Regions Complete.-> STRUCT
    
    STRUCT --> VARFILE
    VARFILE --> MAINCSV
    MAINCSV --> OVERCSV
    OVERCSV -->|Yes| OVEROUT
    OVERCSV -->|No| End
    OVEROUT --> INCONSOUT
    INCONSOUT --> End[✅ COMPLETE]
    
    %% Styling
    classDef inputStyle fill:#e1f5ff,stroke:#01579b,stroke-width:2px
    classDef initStyle fill:#fff3e0,stroke:#e65100,stroke-width:2px
    classDef procStyle fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    classDef outputStyle fill:#e8f5e9,stroke:#1b5e20,stroke-width:2px
    
    class BAM,FASTA,BED,ARGS inputStyle
    class VARFILE,MAINCSV,OVEROUT,INCONSOUT outputStyle