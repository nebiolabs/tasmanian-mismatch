nextflow.enable.dsl=2

params.input_glob = '*.{1,2}.fastq*' //"$baseDir/data/*R{1,2}*.fastq.gz"
params.reference = '' 
params.outdir = "tasmanian_results"
params.bed = ""
params.boundary = "3"
params.softclips = 0

println """
    INPUTS PROVIDED:
    ================
    reads: $params.input_glob 
    genome: $params.reference 
    Bed file: $params.bed (If not provided, Tasmanian-mismatch tables not splitted)
    boundary: $params.boundary (default 3)
    softclips: $params.softclips (default 0)
""".stripIndent()


process align {
    conda "bwa samtools"
    tag { "$prefix" }
    publishDir "${prefix}"
    cpus 8

    input:
        path(read1)

    output:
        tuple path("*.bam"), path("*.csi"), val(prefix), emit: aligned
    
    shell:

   prefix = read1.baseName //.replaceFirst(/1.fastq/,"")

    '''
    read2=$(readlink !{read1} | sed 's/1.fastq/2.fastq/')
    ln -s $read2
    r2=$(echo $read2 | awk -F"/" '{print $NF}')
    pref=$(echo ${r2} | sed 's/2.fastq//')
    bwa mem -t !{task.cpus} !{params.reference} !{read1} ${r2} 2> bwa.log.err | samtools sort -@!{task.cpus} --write-index -o ${pref}.bam
    '''
}

process tasmanian {
    conda "tasmanian-mismatch samtools"
    tag { "${prefix} tasmanian" }
    publishDir "${prefix}"

    input:
        tuple path(bam), path(bai), val(prefix)

    output:
        tuple path("*csv"), path("*html")

    shell:

    '''
    prefix=$(echo !{bam} | sed 's/.bam//') 
    bed=$([ !{params.bed} == "" && echo "" || echo "-b !{params.bed}"])
    
    samtools view !{bam} | run_intersections ${bed} | run_tasmanian -r !{params.reference} -c !{params.boundary} -s !{params.softclips} > ${prefix}.table.csv
    '''
}




Channel
    .fromPath(params.input_glob)
    .filter{ it.name =~ /(1\.fastq)$/ }
    .ifEmpty {error "${params.input_glob} could not find any required file."}
    //.map { filename -> [input_file: filename]}
    .set{ inputChannel }



workflow {
    main:
        inputChannel.view()
        alignedReads = align( inputChannel )
        tasmanians   = tasmanian ( alignedReads.aligned )
}
