#!/usr/bin/env bash

# CONSIDERING that the mapping to the genome was done with bowtie2 (-x argument)
# If using bwa, the extraction of the reference from the bam should change.

# INPUTS
input_bam=$1
output_bam=$2
rescaling_matrix=$3


if [ $# -ne 2 ]; then
   echo  "bash rescale.sh <input_bam_path> <output_bam_filename>"
   exit 1
fi


script_full_path=$(realpath ${BASH_SOURCE[0]})
script_path="$(dirname ${script_full_path})/"
tasmanian_path=$(echo ${script_path} | sed 's/scripts\///')

reference="$(samtools view -H ${input_bam} | grep -o "\-x [^ ]*" | sed 's/\-x //').fa"

# generate mismatch-table
samtools view ${input_bam} | python ${tasmanian_path}/tasmanian_script.py -r ${reference} > mismatch_table.csv
cat mismatch_table.csv | awk -F, '{a=""; for (i=35;i<=50;i++) {a=a","$i} print $1","$2""a}' > rescaling_matrix.csv
python ${script_path}/rescale.py # takes rescaling_matrix.csv and outputs rescaling_factors.csv

# write the new bam file with the base qualities updated, based on the rescaling factors.
samtools view -h ${input_bam} | python ${tasmanian_path}/tasmanian_script.py -r ${reference} --rescale-phred-scores rescaling_factors.csv | samtools view -hb -o ${output_bam}

