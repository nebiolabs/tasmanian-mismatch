#!/usr/bin/env nextflow
//nextflow.preview.dsl=2

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
        file("hs38DH.fa") into genome_file
        file("*.{amb,ann,bwt,pac,sa}") into bwa_index
    
        script:
        """
        bwa.kit/run-gen-ref hs38DH
        bwa.kit/bwa index hs38DH.fa
        """
    }
}
else {
    process find_genome {
    
        input:
        val(genome) from genome_file

        output:
        file("*.{amb,ann,bwt,pac,sa}") into bwa_index
        
        shell:
        '''
        prefix=$(echo !{genome} | sed 's/fasta$//' | sed 's/fa$//')
        files=$(ls ${prefix}* | grep -v "fasta$\\|fa$")
        ln -s ${files} .
        '''
    }
}

/*
        =======================================
            section 1: input is fastq
        =======================================
*/

if ( params.mode == 'fastq' ) {

    println(" working in fastq mode ")

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
        file(index) from bwa_index

        output:
        file("${name}.sorted.bam") into sorted_bams

        script:
        """
        bwa mem -t $task.cpus $genome $reads | samblaster | samtools sort -o ${name}.sorted.bam
        """
    }
}

/*
        =======================================
            section 2: input is bam
        =======================================
*/

else { // implicity, params.mode == 'bam'

    Channel
        .fromPath(params.input)
        .ifEmpty('Could not find the alignment bam files')
        .set { input_bams }

    process bam_sorted {
        
        input:
        file(bam) from input_bams

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
}

/*
        ==================================================
            section 3: subsample bams to eq num of reads
        ==================================================
*/

process get_bam_readNum {

    input:
    file(bam) from sorted_bams

    output:
    stdout into bam_lengths
    stdout into bam_lengthsB // ToDo: Find a more elegant forking
    tuple stdout, file("${bam}") into bam_list

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
    val(length) from bam_lengths_list
    set val(bam_len), file(name) from bam_list

    output:
    file("${bam_name}") into sub_bam
    
    shell:
    '''
    bam_name=$(echo !{name} | sed 's/bam//')subsampled.bam
    if [ !{length} -ne !{bam_len} ]; then
        size=$(echo "!{length} / !{bam_len} + 4" | bc -l)  ## added  +4 for the random integer in samtools
        samtools view -s ${size} ${name} > ${bam_name}
    else 
        ln -s !{name} ${bam_name}
    fi
    '''
}


sub_bam.println()

/*
        ===========================================
            section 3: Tasmanian Artifact Metrics
        ===========================================
*/

/*process Tasmanian {
    
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
}*/