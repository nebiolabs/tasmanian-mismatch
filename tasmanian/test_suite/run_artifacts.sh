input=$1
picard_output=$(echo $1 | sed 's/bam$/picard/')
tasmanian_output=$(echo $1 | sed 's/bam$/tasmanian/')

samtools=/Users/aerijman/Documents/samtools-1.11/samtools


eval "$(conda shell.bash hook)"
conda activate picard
picard CollectSequencingArtifactMetrics I=${input} O=${picard_output} R=genome/small_region.fa

eval "$(conda shell.bash hook)"
conda activate ipython
$samtools view ${input} | python ../tasmanian_script.py -r genome/small_region.fa -o ${tasmanian_output} -z ${tasmanian_output}.summary > ${tasmanian_output}.table
