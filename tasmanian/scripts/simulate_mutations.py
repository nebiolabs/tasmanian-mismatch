#!/usr/bin/env python

from itertools import product
import numpy as np
import re
seed = 42 
import sys, os
from fast_string_replace import replace_at_positions

'''
The idea is to load a mock reference genome (short. About 1000 lines total) and sample from there
read1 and read2 to generate:
    reads without mismatches from the reference.
    reads with a basal level of noise.
    reads with specific noise. For example, for FFPE, we will get reads where read2 has the 
    characteristic pattern of C->T and read1 G->A.
'''


def load_reference_genome(fasta):
    genome, sequence = {}, []

    with open(fasta, 'r') as f:
        for line in f:
            if line[0] == ">":

                if len(sequence) > 5: # skip very first header
                    genome[contig_name] = ''.join(sequence)

                contig_name = line[1:].strip()
                sequence = []

            else:
                sequence.append(line.strip())

    genome[contig_name] = ''.join(sequence)

    return genome


def make_reference_bisulfite(ref_genome, pct_methylation):

    global seed
    np.random.seed(seed)
    
    bisulfite_reference = {}
    for contig, sequence in ref_genome.items():
        CpG_fwd = np.array([i.span()[0] for i in re.finditer("CG", sequence)])
        CpG_rev = np.array([i.span()[0] for i in re.finditer("GC", sequence)])

        # randomly pick the pct_methylation
        np.random.shuffle(CpG_fwd)
        np.random.shuffle(CpG_rev)
        n = ( CpG_fwd.shape[0] + CpG_rev.shape[0] ) * pct_methylation / 100
        n = int(n//2) # roughly modify both strands equally.
        
        methylation_dict_fwd = { int(i):"T" for i in CpG_fwd[:n] }
        methylation_dict_rev = { int(i):"A" for i in CpG_rev[:n] }

        meth_dict = dict(sorted( (methylation_dict_fwd | methylation_dict_rev).items() ) )

        bisulfite_reference[contig] = replace_at_positions(sequence, meth_dict)

    return bisulfite_reference


def read_proposed_variants(vars_file):
    '''
    vars_file is a tab separated file (TSV) and NO HEADER
    each line in vars_file contains original_base, variant_base, frequency(%), strand
    e.g. C  T   10  +
         G  A   0.8 - (yes. 0.8%)
    '''
    with open(vars_file, 'r') as f:
        pre_vars = [i.strip().split("\t") for i in f]


    vars_list = [
        {
        'original':  i[0],
        'variant':   i[1],
        'frequency': float(i[2]),
        'strand':    i[3]
        }
        for i in pre_vars
    ]

    return vars_list


def complement(sequence):
        d = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
        return ''.join([d[i] for i in sequence])

def reverse_complement(sequence):
    return complement(sequence)[::-1]

def assign_contigs_to_variants(vars_list, genome_dict):
    '''
    goal:
        assign variants to contigs randomly.
    input: 
        list of variants (loaded from the file)
        genome (dict)

    output:
        nested dictionaries
        kyes1   = contigs
        keys2   = positions
        values2 = list of characters containing 100 characters
                  e.g. for a G->T with allele frequency=0.1
                  [T,T,G,T,T,T,T,T,T,T...]
    '''
    global seed
    np.random.seed(seed)

    n_vars = len(vars_list)
    n_contigs = len(genome_dict)

    contig_vars = np.random.randint(0, n_contigs, n_vars)

    variants_dict = { contig: {} for contig in genome_dict.keys() }

    for n, (contig, sequence) in enumerate(genome_dict.items()):
        idx_contig_vars = np.where(contig_vars == n)[0]
        
        for idx_contig_var in idx_contig_vars:
            var_dict = vars_list[idx_contig_var]
            strand = var_dict['strand']
            old_base = [var_dict['original'] if strand == "+" else reverse_complement(var_dict['original'])][0]
            new_base = [var_dict['variant'] if strand == "+" else reverse_complement(var_dict['variant'])][0]
            positions_var = np.array( [i.span()[0] for i in  re.finditer(old_base, sequence)] )
            selected_position = positions_var[ np.random.randint( len(positions_var) ) ]
            floored_frequency = np.max( [var_dict['frequency'], 1] ).astype(int)
            
            # variant and original with their frequencies in a arr[10] 
            variants_dict[contig][ selected_position ] = np.array( [old_base] * (100 - floored_frequency) + [new_base] * floored_frequency )
            np.random.shuffle( variants_dict[contig][ selected_position ] )

            # print(f"Assigned variant {var_dict['original']}->{var_dict['variant']} at contig {contig} position {selected_position} with frequency {var_dict['frequency']}%")
            # print(f"floored frequency: {floored_frequency}")
    
    # e.g. variant_dict['chr1'][1014365] = np.arr( ['G','G','T','G','G','T','T','G','G',....'T'] )
    return variants_dict


def generate_reads(genome_dict, 
                   variants_dict, 
                   library_type='FFPE', 
                   n_reads=100, 
                   read_length=76, 
                   insert_length=100, 
                   insert_length_var=10, 
                   base_noise_level=0.0001):

    global seed
    np.random.seed(seed)
    # insert length = WITHOUT adapter

    simulated_reads = []

    # Distribute the reads across contigs evenly, based on the length of each contig.
    contig_lengths = {name: len(seq) for name, seq in genome_dict.items()}
    total_length = np.sum([v for v in contig_lengths.values()])
    read_unit = n_reads / total_length

    contig_reads = {name: int(length * read_unit) for name, length in contig_lengths.items()}
    previous_n = 0

    # Generate reads for each contig.
    for contig_name, original_sequence in genome_dict.items():
        read_starts = np.random.randint(0, contig_lengths[contig_name] - read_length, contig_reads[contig_name])
        insert_lengths = np.random.randint(0, insert_length_var, contig_reads[contig_name])

        # repeat the sampling N times (Ideally 100 but let's start with 10)
        # split array of starting positions in 10 (to make the sampling form "10 var-genomes" easier.
        read_starts_split = np.array_split(read_starts,10)
        insert_lengths_split = np.array_split(insert_lengths, 10)

        # 10 times, generate a "variant-genome" and sample from it.
        for var_index in range(10):
            contig_vars_input = { int(var_position): str(var_nt[var_index]) for var_position, var_nt in variants_dict[contig_name].items() }

            sequence = replace_at_positions(original_sequence, contig_vars_input)

            # generate the reads. You can only see the ends of the fragments (initially both ends are fwd for simplicity)
            reads  = { 
                n: {
                    'left' : sequence[i                                    : i+read_length                                   ].upper(),
                    'right': sequence[i+insert_lengths_split[var_index][n] : i+insert_lengths_split[var_index][n]+read_length].upper() 
                }
                for n,i in enumerate(read_starts_split[var_index])
            }

            # randomly half of thes e reads will be: read1=rev and read2=fwd (and vice versa)
            n_reads_tmp = len(read_starts_split[var_index])
            idx = np.arange( n_reads_tmp ) 
            np.random.shuffle(idx)
            fr = idx[:n_reads_tmp//2]
            rf = idx[n_reads_tmp//2:]

            # make the necessary reads reverse
            new_reads = {n  :{} for n in range(previous_n, n_reads_tmp + previous_n)} # n*(var_index+1) to have unique read IDs across var-genomes

            sys.stderr.write(f"fr reads: {len(fr)}\trf reads: {len(rf)}\n")

            for n in fr:
                new_reads[ n + previous_n ] = {
                        'R1': reads[n]['left'],
                        'R2': complement( reads[n]['right'][::-1] )
                }
            for n in rf:
                new_reads[ n + previous_n ] = {
                        'R1': complement( reads[n]['left' ] )[::-1],
                        'R2': reads[n]['right']            
                }

            simulated_reads.append(new_reads)

            previous_n += n_reads_tmp

    return {k:v for d in simulated_reads for k,v in d.items()}


def print_reads(reads):
    global seed
    '''
    If needed, this can be adapted to different sequencing platforms, 
    barcodes, tails, and proper coordinates
    '''
    np.random.seed(seed)
    phreds = [ ''.join( ['A'] * 5 + ['E'] * ( len(reads[0]['R1']) - 5 ) ) ] * len(reads)
    tiles = np.sort( np.random.randint( 1,5, len(reads) ) )
    coords = np.random.randint( 100, 14000, (len(reads),3) )
    barcodes = "GTTCTGCA+GCTCCTTC"
    headers = [
        f"@NB552064:NB552064:H3CTYAFXC:{tile}:{coord[0]}:{coord[1]}:{coord[2]} 1:N:0:{barcodes}"
        for tile, coord in zip(tiles, coords)
    ]

    preads1 = '\n'.join( 
            [
                '\n'.join([ header, read['R1'], "+", phred[ :len(read['R1']) ] ]) 
                for header, read, phred in zip(headers, reads.values(), phreds)
            ] 
    ) 
    
    preads2 = '\n'.join( 
            [
                '\n'.join([ header, read['R2'], "+", phred[ :len(read['R2']) ] ]) 
                for header, read, phred in zip(headers, reads.values(), phreds)
            ] 
    ) 

    return preads1, preads2


def generate_probabilities(library_type='FFPE', read_length=76):
    '''
    INPUT:
        type of library / pattern of positional mutations -> str
        length of the read -> int
    OUTPUT:
        2 dictionaries for read1 and read2
        each dict has mismatches as keys (e.g. CT=C->T, GT, etc...)
        and values are numpy arrays with shape = (read_length, 1)
        and values indicating the probabilidad of each position to 
        make it into the library.
        E.g. a CT at position 3 in read 1, with a value of 0.8 will be adopted 
        (by other function) in the read, and the C will become T in 
        ~80% of the positions 3 from read 1.
    '''
    bases = ['A','T','C','G']
    mismatches = [''.join(i) for i in product(bases, repeat=2) if i[1]!=i[0]]
    read1 = {mismatch: np.zeros(read_length) for mismatch in mismatches}
    read2 = read1.copy()

    if library_type == 'FFPE':
        
        # C to T in read 2
        # ---------------

        # the curve seems to comprise 2 distinct regions
        x1 = np.arange(1,12)
        x2 = np.arange(12,77)
        
        # with values within the following limits
        y1 = np.linspace(0.009, 0.004, len(x1)) 
        y2 = np.linspace(0.004,0.0015, len(x2))

        # data has some noise:
        var1 = 0.0004 # defined empirically, guided by data
        var2 = 0.0008
        noise1 = ( np.random.rand(len(y1)) - 0.5 ) * var1
        noise2 = ( np.random.rand(len(y2)) - 0.5 ) * var2

        # So combined all
        # x = np.hstack([x1, x2])
        y = np.hstack([ y1 + noise1, y2 + noise2])

        # Visualize the data if needed
       #plt.figure(figsize=(5,5))
        #plt.scatter(x,y, c='k', s=3)
        #plt.grid()
        #plt.ylabel('mock normalized mismatch counts')
        #plt.xlabel('read position');

        read2['CT'] = y

        
        # G to A in read 1
        # ----------------
        
        # linear with a mild exponential increase
        x = np.arange(1,77)

        base = np.linspace(0.0013, 0.0015,len(x))
        expon = -1 * np.linspace(1.1,1.5, len(x))
        y = np.power(base, expon) / (base[-1]**expon[-1]) * 0.003

        # add a little noise
        var2 = 0.0008
        noise2 = ( np.random.rand(len(y)) - 0.5 ) * var2
        y += noise2  

        # Visualize the data if needed    
        #plt.figure(figsize=(5,5))
        #plt.scatter(x, y, c='k', s=3)
        #plt.ylim(0,0.009)
        #plt.grid()

        read1['GA'] = y
    

    elif library_type == 'sonication':
        pass
    else:
        pass

    return read1, read2


def mutate_single_read(read, read_probs, random_probs):
    # read_list = list(read)

    '''
    probs = dict
        keys   : AT, GA, etc. and
        values : probabilities
    read = numpy (str) array

    copilot benchmarked set operation slightly faster 
    than numpy-indices AND operation for this type of 
    cases
    '''
    for m,p in read_probs.items():
        # looking for positions where nt is X (for mutation XY)
        # and random robability of XY is high enough.
        base_position = set([i.span()[0] for i in re.finditer(m[0], read)])
        read_probs_complement = 1 - p
        prob_position = set( np.where(random_probs >= read_probs_complement)[0] )
        positions = base_position.intersection(prob_position)

        
       # for i in positions:
       #      read_list[i] = m[1] 

        if len(positions) == 0: continue

        read = replace_at_positions( read, {int(i):m[1] for i in positions} )

    return read # ''.join(read_list)


def mutate_reads(reads, probabilities, read_length):
    simulated_reads= {}

    for n, read in enumerate(reads.values()):
        simulated_reads[n] = {
            'R1': mutate_single_read(read['R1'], probabilities[0], np.random.rand(read_length)),
            'R2': mutate_single_read(read['R2'], probabilities[1], np.random.rand(read_length))
        }
    return simulated_reads

