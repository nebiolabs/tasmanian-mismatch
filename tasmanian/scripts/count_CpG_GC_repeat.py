#!/bin/python

'''
    The idea is to observe patterns within the first 20 nucleotides of reads
    that intersect repeats. Why do we see these artifacts enriched in the 
    5' end of read 2?
    
    The results consist on matrices where the first dimention is 1 to 20, 
    representing the first 20 nucleotides from the splitting point between 
    the region of the read that intersects a repeat and the region that does
    not. The other dimentions are the nucleotide types.

    ALSO, include 20mer map of CG and GC (I don't yet fully understand wether
    methylation on one strand could affect the opposite strand or not). 
'''



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

    print('reading genome...')
	genome = read_genome(genome_file)
    os.system('clear')
	print('read genome... DONE!')

    # save statistics on a category basis on read-intersection category
    # but also on a repeat family category
    repeat_families = [
		'DNA',
		'DNA?',
		'LINE',
		'LTR',
		'LTR?',
		'Low_complexity',
		'RC',
		'RC?',
		'RNA',
		'Retroposon',
		'SINE',
		'SINE?',
		'Satellite',
		'Simple_repeat',
		'Unknown',
		'rRNA',
		'scRNA',
		'snRNA',
		'srpRNA',
		'tRNA'
    ]

    positional = {'A':0,'C':0,'T':0,'G':0}
	basic_dict = {'GC':dict(zip(range(20),[positional]*20)), 
                  'CpG':np.zeros(20)}
    
    results = {i: basic_dict for i in repeat_families}
    
        
'''
THIS HAS BEEN USED TO EVALUATE GC AND CPG CONETNT IN REPEAT REGIONS WITHOUT
ANY CORRELATION TO READS FROM OUR DATA

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
'''





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
