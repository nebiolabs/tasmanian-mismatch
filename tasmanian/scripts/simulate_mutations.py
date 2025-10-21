#!/usr/bin/env python

from itertools import product
import numpy as np
import re
seed = 42 

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


def complement(sequence):
        d = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
        return ''.join([d[i] for i in sequence])


def generate_reads(genome_dict, library_type='FFPE', n_reads=100, read_length=76, insert_length=100, insert_length_var=10, base_noise_level=0.0001):
    global seed
    np.random.seed(seed)
    # insert length = WITHOUT adapter

    # Distribute the reads evenly, based on the length of each contig.
    contig_lengths = {name: len(seq) for name, seq in genome_dict.items()}
    total_length = np.sum([v for v in contig_lengths.values()])
    read_unit = n_reads / total_length

    contig_reads = {name: int(length * read_unit) for name, length in contig_lengths.items()}

    for contig_name, sequence in genome_dict.items():
        read_starts = np.random.randint(0, contig_lengths[contig_name] - read_length, contig_reads[contig_name])
        insert_lengths = np.random.randint(0, insert_length_var, contig_reads[contig_name])

        # generate the reads (initially both ends are fwd)
        reads  = { 
            n: {
                'left' : sequence[i               : i+read_length].upper(),
                'right': sequence[i+insert_lengths[n] : i+ insert_lengths[n] +read_length].upper() 
            }
            for n,i in enumerate(read_starts)
        }
        
        # randomly half of these reads will be: read1=rev and read2=fwd (and vice versa)
        n_reads_tmp = len(read_starts)
        idx = np.arange( n_reads_tmp ) 
        np.random.shuffle(idx)
        fr = idx[:n_reads_tmp//2]
        rf = idx[n_reads_tmp//2:]

        # make the necessary reads reverse
        new_reads = {n:{} for n in range(len(read_starts))}
        
        print(f"fr reads: {len(fr)}\trf reads: {len(rf)}")

        for n in fr:
            new_reads[n] = {
                    'R1': reads[n]['left'],
                    'R2': complement( reads[n]['right'][::-1] )
            }
        for n in rf:
            new_reads[n] = {
                    'R1': complement( reads[n]['left' ] )[::-1],
                    'R2': reads[n]['right']            
            #new_reads[n] = {
            #    'R2': complement( reads[n]['left' ] ),
            #    'R1': reads[n]['right'][::-1]
            }

    return new_reads


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
                '\n'.join([header, read['R1'], "+", phred]) 
                for header, read, phred in zip(headers, reads.values(), phreds)
            ] 
    ) 
    
    preads2 = '\n'.join( 
            [
                '\n'.join([header, read['R2'], "+", phred]) 
                for header, read, phred in zip(headers, reads.values(), phreds)
            ] 
    ) 

    return preads1, preads2


def generate_probabilities(library_type='FFPE', read_length=76):

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
    read_list = list(read)
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

        for i in positions:
            read_list[i] = m[1] 

    return ''.join(read_list)


def mutate_reads(reads, probabilities, read_length):
    simulated_reads = {}
    for n, read in enumerate(reads.values()):
        simulated_reads[n] = {
            'R1': mutate_single_read(read['R1'], probabilities[0], np.random.rand(read_length)),
            'R2': mutate_single_read(read['R2'], probabilities[1], np.random.rand(read_length))
        }
    return simulated_reads

