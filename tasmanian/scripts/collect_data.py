#!/bin/python

import sys,os,re
sys.path.append(os.path.abspath('../'))
from utils import read_genome
import gzip
import pandas as pd
import numpy as np

for n,i in enumerate(sys.argv):
    if i=='-i':
        in_file = sys.argv[n+1]


def main():
    
    # read genome into dictionary
    genome = read_genome('../hg38')

    # read repeats into table
    df = load_repeats('../bed_files/RepeatMasker.bedGraph')
    
    # get data from sam.gz file and store important values in a array
    results = get_data(in_file, genome)

    results = pd.DataFrame(results)
    results.columns = ['header','category', 'rep_class', 'read_number', 'gc_fragment', 'gc_read', 'cpg_fragment', 'cpg_read', 'lenght', 'noN_lenght']

    # convert to numeric what's necessary
    results[['read_number','gc_fragment','gc_read','cpg_fragment','cpg_read','lenght','noN_lenght']] = results[[
            'read_number','gc_fragment','gc_read', 'cpg_fragment','cpg_read','lenght','noN_lenght']].astype(int)

    results.to_pickle(in_file[:-6] + 'pkl')



# read repeats and check if there is any enrichment of CpGs or GC-content
def load_repeats(repeats_bed_file):
    repeats = []
    with open(repeats_bed_file,'r') as f:
        while True:
            try:
                line = next(f)
                chrom, start, end, strand, a, b, c = line.strip().split('\t')
                repeats.append([chrom, start, end, strand, a, b, c])
            except StopIteration:
                break

    df = pd.DataFrame(repeats)
    df.columns = ['chrom','start','end','strand','name','class','family']
    df.start = df.start.astype(int)
    df.end = df.end.astype(int)
    return df



def get_data(filename, genome):
    '''
        open sam.gz file and extracts categories and lengths
        into a dictionary with categories as keys
    '''
    
    def revcomp(seq):
        d = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 
             'W':'W', 'M':'M', 'U':'U', 'S':'S', 'K':'K', 
             'R':'R', 'Y':'Y', 'B':'B', 'D':'D', 'H':'H', 'V':'V', 'Z':'Z'}
        return ''.join([d[i] for i in seq][::-1])
        
    n=0
    results = []
    
    with gzip.open(filename, 'rt') as f:
        for line in f.read().strip().split('\n'):
            try:
                header, flag, chrom, coordA, _0, _cigar, _1, coordB, lenght, seq, qual, category, repeat = line.split('\t')
                category = re.sub('^.*:.*:', '', category)
                strand, re_family, rep_class, rep_name = repeat.split(':')
                coordA, coordB, lenght = int(coordA), int(coordB), np.abs(int(lenght))
                start, end = [coordA, coordA+lenght] if flag in ['99','163'] else [coordB, coordB+lenght]
                read_number = 1 if flag in ['99','86'] else 2
                
                fragment_seq = genome[chrom][start:end] if strand == '+' else revcomp(genome[chrom][start:end])
                noN_lenght = len([i for i in seq if i!='N'])
                
                gc_fragment = fragment_seq.count('G') + fragment_seq.count('C')
                gc_read = seq.count('G') + seq.count('C')
                cpg_fragment = fragment_seq.count('CG')
                cpg_read = seq.count('CG')
                
                results.append([header, category, rep_class, read_number, gc_fragment, gc_read, cpg_fragment, cpg_read, lenght, noN_lenght])
            
            except Exception as e:
                print(str(e))
                continue
                
    return results

if __name__=='__main__':
    main()
