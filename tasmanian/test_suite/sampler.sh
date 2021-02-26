# original reads

# read1
# AAGTGAGCTAGCAAATAAATTTCCCTATGGTCTCAGCTCTGAGTGGAGAGAGAAAATGTTCCCTGTGGAGTTTATA
# read2
# TGAGCTTGCTTTGATGATTTATTTGTCCAGAGAGGATTTTTTTTCCTACCTAGCATTTTGGACTGCTATCAACCTG


bash bwa_and_RG.sh original_read1.fastq original_read2.fastq original.bam

# MUTATIONS in different combinations of READS and ORIENTATIONS
# =============================================================
bash bwa_and_RG.sh <(cat original_read1.fastq | sed 's/AAGTGAGCTAG/AAGTGAGTTAG/') original_read2.fastq C-T_end.read1.fwd.bam
bash bwa_and_RG.sh <(cat original_read1.fastq | sed 's/AAGTGAGCTAG/AAATGAGCTAG/') original_read2.fastq G-A_end.read1.fwd.bam
bash bwa_and_RG.sh <(cat original_read1.fastq | sed 's/TGGTCTCAGCT/TGGTTTCAGCT/') original_read2.fastq C-T.read1.fwd.bam
bash bwa_and_RG.sh <(cat original_read1.fastq | sed 's/TGGTCTCAGCT/TAGTCTCAGCT/') original_read2.fastq G-A.read1.fwd.bam

bash bwa_and_RG.sh original_read1.fastq <(cat original_read2.fastq | sed 's/GAGCTTGCTTTG/GAGTTTGCTTTG/') C-T_end.read2.rev.bam
bash bwa_and_RG.sh original_read1.fastq <(cat original_read2.fastq | sed 's/GAGCTTGCTTTG/AAGCTTGCTTTG/') G-A_end.read2.rev.bam
bash bwa_and_RG.sh original_read1.fastq <(cat original_read2.fastq | sed 's/CCAGAGAGGATT/TCAGAGAGGATT/') C-T.read2.rev.bam
bash bwa_and_RG.sh original_read1.fastq <(cat original_read2.fastq | sed 's/CCAGAGAGGATT/CCAAAGAGGATT/') G-A.read2.rev.bam

# FWD is now READ2 and REV is now READ1
# -------------------------------------
bash bwa_and_RG.sh <(cat original_read2.fastq | sed 's/2:N:0/1:N:0:/' | sed 's/CCAGAGAGGATT/TCAGAGAGGATT/')\
   			       <(cat original_read1.fastq | sed 's/1:N:0:/2:N:0:/') C-T.read1.rev.bam
bash bwa_and_RG.sh <(cat original_read2.fastq | sed 's/2:N:0:/1:N:0:/' | sed 's/CCAGAGAGGATT/CCAAAGAGGATT/')\
   			       <(cat original_read1.fastq | sed 's/1:N:0:/2:N:0:/') G-A.read1.rev.bam
bash bwa_and_RG.sh <(cat original_read2.fastq | sed 's/2:N:0:/1:N:0:/')\
   			       <(cat original_read1.fastq | sed 's/1:N:0:/2:N:0:/' | sed 's/TGGTCTCAGCT/TGGTTTCAGCT/') C-T.read2.fwd.bam
bash bwa_and_RG.sh <(cat original_read2.fastq | sed 's/2:N:0:/1:N:0:/')\
	<(cat original_read1.fastq | sed 's/1:N:0:/2:N:0:/' | sed 's/TGGTCTCAGCT/TAGTCTCAGCT/') G-A.read2.fwd.bam


# SOFTCLIPS
# =========
# Included in the mutation files

# FLAGS
# =====
samtools=/Users/aerijman/Documents/samtools-1.11/samtools
# duplucates
$samtools view -h original.bam | awk 'BEGIN{OFS="\t";}{if ($2==99) $2=1123; else if ($2==147) $2=1171; print $0}' | $samtools view -b > flag-duplicates.bam
# supp-alignment
$samtools view -h original.bam | awk 'BEGIN{OFS="\t";}{if ($2==99) $2=2147; else if ($2==147) $2=2195; print $0}' | $samtools view -b > flag-supp-align.bam
# unmapped
#$samtools view -h original.bam | awk 'BEGIN{OFS="\t";}{if ($2==99) $2=73; else if ($2==147) $2=133; print $0}' | $samtools view -b > flag-unmap.bam
r1=AAGTGAGCTAGCAAATAAATTTCCCTATGGTCTCAGCTCTGAGTGGAGAGAGAAAATGTTCCCTGTGGAGTTTATA
r2=CAGGTTGATAGCAGTCCAAAATGCTAGGTAGGAAAAAAAATCCTCTCTGGACAAATAAATCATCAAAGCAAGCTCA
bash bwa_and_RG.sh <(cat original_read1.fastq | sed "s/$r1/$(echo $r1 | rev)/")\
	<(cat original_read2.fastq original.bam | sed "s/$r2/$(echo $r2 | rev)/") flag-unmap.bam
