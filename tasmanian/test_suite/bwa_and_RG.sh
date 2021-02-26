#!/usr/bin/env bash

r1=$1
r2=$2
bam=$3

eval "$(conda shell.bash hook)"
conda activate aligners
bwa mem genome/small_region.fa $r1 $r2 | head -4 > tmp.bam

conda deactivate && conda activate picard
picard AddOrReplaceReadGroups I=tmp.bam O=$bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

rm tmp.bam
