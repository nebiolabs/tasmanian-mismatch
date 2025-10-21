
python test_runs.py
bwa mem -t4 genome_simul_test simulated_reads_R1.fastq simulated_reads_R2.fastq > aligned.sam
samtools sort --write-index aligned.sam -o aligned_sorted.bam
samtools view aligned_sorted.bam | awk '{arr[$2]++}END{for (i in arr) {print i"\t"arr[i]} }'
samtools view aligned_sorted.bam | run_tasmanian -r test_data/test_genome.fa > aligned_sorted_tasmanianed.csv
