#!/usr/bin/bash

bash sampler.sh

for i in *error_summary_metrics; do 
	name=$(echo $i | sed 's/\.picard\.error_summary_metrics//')
	printf "$name \n"
	t=$(echo $i | sed 's/picard\.error_summary_metrics/tasmanian\.table/')
	./summary $t; cat $i | sed -n '8,13p'
done

for i in *bam; do 
	bash run_artifacts.sh $i $(echo $i | sed 's/bam/picard/') $(echo $i | sed 's/bam/tasmanian/')
done 2> log-run_artifacts
