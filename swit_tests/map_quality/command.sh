#!/bin/bash

samtools view AG.bam | \ 
python tasmanian.py -r genome/small_region.fa -m 21 | \
tail -n +7 | awk -F, '{for(i=3;i<=NF;i++) t+=$i; print t; t=0}'
