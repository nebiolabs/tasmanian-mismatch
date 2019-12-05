#!/bin/python

import os,sys,re
import numpy as np
import pandas as pd
import pickle

genome_file = '../hg38'
output_file = 'test.output'

for n,i in enumerate(sys.argv):
	if i=='-o':
		output_file = sys.argv[n+1]
	if i=='g':
		genome_file = sys.argv[n+1]



def main():

	genome = read_genome(genome_file)
	print('read genome... DONE!')

	results = {'GC':[],
			   'CpG':[]
			  }

	for line in sys.stdin:

		chrom, start, end, strand = line.strip().split('\t')[:4]
		# ['chr1', '79099745', '79099796', '+', '(ATTTA)n', 'Simple_repeat', 'Simple_repeat']
		start = int(start)
		end = int(end)

		if chrom not in genome:
			continue		

		CpG = -1
		GC = -1

		seq = genome[chrom][start:end]

		if strand == '+':
			CpG = seq.count('CG')
		elif strand=='-':
			CpG = seq.count('GC')
		else:
			print('strand not recognized ',strand)

		GC = seq.count('C') + seq.count('G')

		results['CpG'].append(CpG/len(seq))
		results['GC'].append(GC/len(seq))
		
		
	with open(output_file + '.pkl', 'wb') as f:
		pickle.dump(results, f)






def read_genome(ref_file):
    '''
    reads genome fasta file into dictionary of sequences with keys=chromosomes
    INPUT: fasta filename
    OUTPUT: dictionary. keys=chromsomes, 
            values=sequence of chromosome
    '''
    reference = {}
    with open(ref_file) as f:
        while True:
            try:
                line = next(f).strip()
                if line[0] == ">":
                    if 'current_key' in locals():
                        reference[current_key] = ''.join(tmp)
                    current_key = line.split(' ')[0][1:].replace(' ','')
                    tmp = []

                else: 
                    tmp.append(line.strip())

            except StopIteration:
                reference[current_key] = ''.join(tmp)
                break
    return reference




if __name__ == '__main__': 
	main()
