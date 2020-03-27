#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
============================================================================
                    DNA DAMAGE NEXTFLOW PIPELINE
============================================================================
# homepage documentation: 

# Authors:
    Brad Langhorst
    Ariel Erijman
============================================================================
*/

params.input = "$baseDir/data/*R{1,2}*.fastq.gz"
params.genome = '' //"$baseDir/data/grch38.fa"   // copy from bwakit the feature to fetch the genome, also allow for w-w/o alt-contigs
params.outdir = "$baseDir/tasmanian_results"
params.date = "Not deided yet"
params.saveReference = false
params.bed = "$baseDir/bedfile"

// understand if the input is fastq or bam. if Bam, flag to skip alignment
params.mode = ( params.input =~/fastq/ ? 'fastq' : 'bam' )

println """
    reads: $params.input
    genome: $params.genome
    data saved in $params.outdir
""".stripIndent()

/*
        ========================================================
            section 1: if necessary, download reference genome
        ========================================================
*/
genome_file = file(params.genome)
if (params.genome == '') {

    process build_genome {

        output:
        file(ref) into genome_file
        file "*.{amb,ann,bwt,pac,sa}" into bwa_index
    
        script:
        """
        bwa.kit/run-gen-ref hs38DH
        bwa.kit/bwa index hs38DH.fa
        """
    }
}

/*
        =======================================
            section 1: input is fastq
        =======================================
*/

if ( params.mode == 'fastq' ) {
    // read_pairs_ch join all R1,R2 that match the prefix
    Channel
        .fromFilePairs( params.input )
        .ifEmpty('Could not find the pair of fastq files')
        .set { read_pairs_ch }

    // two bash codes in the same process stinks. BUT ADD CHECKS!!!
    process AlignTOgenome {
        
        input:  
        file genome from genome_file
        set val(name), file(reads) from read_pairs_ch    
        file(bwa_index) from bwa_index

        output:
        file("${name}.sorted.bam") into aligned_bams

        script:
        """
        //echo ${params.mode_fastq} | grep -q "true"  &&\
        bwa mem -t $task.cpus $genome $reads | samblaster | samtools sort -o ${name}.sorted.bam
        
        #echo ${name} | grep -q "sort" && samtools sort ${name} -o ${name}.sorted.bam
        
        #samtools view ${name}.sorted.bam | wc -l | tr -d " \n" 
        #printf "${name}," && samtools view ${name}.sorted.bam | wc -l 
        """
    }
}

/*
        =======================================
            section 2: input is bam
        =======================================
*/

else { // implicity, params.mode == 'bam'

    process bam_sorted {
        
        input:
        file(bam) from params.input

        output:
        file("${sorted_bam}") into sorted_bams

        shell:
        '''
        sorted_bam=$(echo !{bam} | sed 's/bam$$/sorted\\.bam/')
        
        samtools view -H !{bam} | grep "SO" | sed 's/.*\tSO:\\(.*\\)$/\\1/' |\
        grep -q "coordinate\\|queryname" &&\
        ln -s !{bam} ${sorted_bam} ||\
        samtools sort !{bam} -o ${sorted_bam}
        '''
    }

    process get_bam_readNum {
        input:
        file(bam) from sorted_bams

        output:
        stdout into bam_lengths
        set val(stdout), file(bam) into bam_list
    
        script:
        """
        samtools view $bam | wc -l | tr -d " \n"
        """
    }

    /* Channel holds until all bam_lengths are collected, then stores the min
       number of reads. This is to subssample all files to that number. */
    Channel
        .from bam_lengths
        .min { it }
        .subscribe { println "the value is  $it **" }
        .set { bam_lengths_list }

    process subsample_bams {
        
        input:
        val(lenght) from bam_lengths
        set val(bam_len), file(name) from sorted_bam

        output:
        file("${bam}") into sub_bam
        
        shell:
        '''
        size=$(echo !{length} / !{bam_len} + 4 | bc)  ## added  +4 for the random integer in samtools -s
        samtools view -s ${size} ${name} > ${bam}
        '''
    }
}

/*
        ===========================================
            section 3: Tasmanian Artifact Metrics
        ===========================================
*/

process Tasmanian {
    
    input:
    file(bam) from sub_bam
    file(bed) from file(params.bed)
    file(genome) from genome_file
    //set val(name), file(reads) from read_pairs_ch // This is only to retrieve name
    val(num_reads) from bam_lengths_list

    //output:

    shell:
    '''  
    # INTERESETING: samtools idx into text tells how many reads.
    # instead of samtools head, calculate the percentage and run subsample.
    name=$(echo !{bam} | sed 's/bam//')
    samtools view !{bam} | head -n !{num_reads} | run_intersections -b !{bed} | run_tasmanian -r !{genome} >!{baseDir}/${name}.table.csv
    ''' 
}
