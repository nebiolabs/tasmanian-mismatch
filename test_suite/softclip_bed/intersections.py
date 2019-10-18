#!/bin/python
'''
    input sam file from stdin
    bedfile as one of the arguments (chromosome, start, end - delimiter=tab)

'''
import os, sys
import numpy as np

# Define outputs: 1- SAM lines including ONLY the intersecting part
#                 2- SAM lines including ONLY the NON intersecting part
#                 3- Statistics about these intersections

intersections = []  # part of reads in common between repeats bed regions and bam read
differences = []    # part of reads NOT in commin between repeats bed regions and bam read
statistics = []     # length of the intersection is kept 
unrelateds = []     # reads without ANY intersection with bed regions are written in another file

# load arguments
for n,i in enumerate(sys.argv):

    if i in ['--bed','-b','--bed-file']:
        bedfile = sys.argv[n+1]

    if i in ['--output', '-o']:
        out_prefix = sys.argv[n+1] + "."


def main():

    # load bed_file
    bed, bed_lens, total_bed_lens = read_bed(bedfile)

    #sys.stderr.write('bed file loaded...')
    #print('bed file loaded... not error')

    # to avoid going over the entire chromosome on each read, 
    # keep an updated index of the bedfile to start from
    bed_index, total_beds = 0,0 
    last_chrom = 'chr1'
    finito=False # when finished with the bed file, FINISH!
    skip_chrom=False # if no more reads in chrom, skip fast reads with that chrom in bam

    # read bam file from stdin and 
    for line in sys.stdin:
        line = line.strip('\n') # Don't repeat same opreation over and over on each iteration.

        # if finished looping over the bed file, don't read whole bam file.
        if finito:
            unrelateds.append(line)
            continue # still collect the unrelateds 
        if len(line) < 50: break # avoid potential empty lines at the end.

        _1, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred = line.split('\t')[:11]
        
        skip_read = False
        seq_len = len(seq)
        mapq = int(mapq)
        flag = int(flag)
        tlen = int(tlen)
        start = int(start) #-1 --> correlated with line #139
        end = start + seq_len

        # if chromosomoe not in bed, read is not intereseting also other poor situations added here
        if chrom not in bed: # or mapq<20 or np.abs(tlen)<100 or np.abs(tlen)>500: 
            unrelateds.append(line)
            continue

        # only consider uniquely mapped and proper pair
        if flag not in [163,99,147,83]: continue

        # ================================================================================================
        # Check-point
        '''
            If read is first on the pair, save it to memory until the paired read is found.
            If read is second on the pair, we already have the first in the memory, in a list of first
            reads (sometimes they are 100s of lines appart) and we compute all what follows for both reads.
        '''
        
        # ================================================================================================

        # new chromosome -> zero bed_index!
        if chrom != last_chrom:
            skip_chrom=False
            bed_index=0
            if total_beds >= total_bed_lens: # Have we finished reading the bed file?
                finito = True
                break

        # if read is ahead of bed fragment, advance index
        while start >= bed[chrom][bed_index,1]:
            if bed_index == bed_lens[chrom]-1: # CAREFULL, stop if finished reading bed file.
                skip_chrom = True
                break
            bed_index +=1
            total_beds+=1

        # don't waste time with bam reads with no partner in bed
        if skip_chrom: 
            unrelateds.append(line)
            continue

        # continue until a read intersects the bed fragment
        if end <= bed[chrom][bed_index,0]: 
            unrelateds.append(line)
            continue

        # capture the intersection
        # 1 read      ------------------- 
        #   bed        a  =========  b
        # 2 read              -----------
        # 3 read      --------
        # 4 read           ------
        ''' 
        cases 2 and 3 with only 1 base in common (...a>=0) is equivalent
        to no intersection since I already added a number to the bed-end
        to include the upper value and the read goes from zero and should
        not include the upper value. a=0 and b=0 is included in case 4.
        '''

        a = bed[chrom][bed_index][0] - start
        b = end - bed[chrom][bed_index][1]

        # cigar is to correlate clips with intersections.
        ''' It might happen that most intersections are not correctly mapped
        and therefore, are softcliped'''
        CIGAR = expand_cigar(cigar) 

        intersect_size=None
        inters_seq, diff_seq, diff_seq2 = None, None, None

        if a>0 and b>0:
            diff_seq = seq[:a] + ''.join(['N']*(76-a))
            diff_seq2 = ''.join(['N']*(76-b)) + seq[-b:]

            inters_seq = ''.join(['N']*a) + seq[a:-b] + ''.join(['N']*b)
            intersect_size = seq_len-b-a

        elif a>0 and b<0:  # elif
            diff_seq = seq[:a] + ''.join(['N']*(76-a))

            inters_seq = ''.join(['N']*a) + seq[a:]
            intersect_size = seq_len-a

        elif b>0 and a<0:
            B = seq_len - b
            diff_seq = ''.join(['N']*B) + seq[B:]

            inters_seq = seq[:B] + ''.join(['N']*b)
            intersect_size = B
        
        else:
            inters_seq = seq
            intersect_size = seq_len

        # write SAM reads containing intersections or non-intersecting sequences
        flag, mapq, tlen, start = str(flag), str(mapq), str(tlen), str(start) #+1)
            
        if diff_seq != None:
            differences.append('\t'.join([_1, flag, chrom, start, mapq, cigar, _2, _3, tlen, diff_seq, phred]))
        if diff_seq2 != None:
            differences.append('\t'.join([_1, flag, chrom, start, mapq, cigar, _2, _3, tlen, diff_seq2, phred]))
        
        intersections.append('\t'.join([_1, flag, chrom, start, mapq, cigar, _2, _3, tlen, inters_seq, phred]))
        
        #if int(intersect_size) > 76: print(seq, bed[chrom][bed_index], flag, start, end, tlen, intersect_size, a, b)
        if len(inters_seq)<1: print(seq, bed[chrom][bed_index], a, b, inters_seq, start, end, chrom)

        statistics.append(str(intersect_size))
        last_chrom = chrom

    # save files:
    for data, filename in zip([intersections, differences, statistics, unrelateds],
                              ['intersections.sam','differences.sam','statistics.txt', 'unrelateds.sam']
                             ):
        with open(out_prefix + filename,'w') as f:
            f.write('\n'.join(data))


def expand_cigar(cigar):
    cigarString = []
    num=''
    for i in cigar:
        if i.isnumeric(): 
            num = num + i
        else:
            cigarString.append(''.join([i]*int(num)))
            num=''

    return ''.join(cigarString)


def read_bed(bedfile):
    '''
        read bed file
        INPUT: bedfile path <name>
        OUTPUT: 1. bed dictionary (chromosome=keys, [#start, #end]=values (np.array. each chrom sorted on #starts) 
                2. dictionary of lengths of arrays for each chromosome
                3. total length of bed elements.
    '''
    bed = {}
    with open(bedfile) as f:
        while True:
            try:
                chrom, start, end = next(f).strip().split('\t')[:3]

                if chrom not in bed:
                    bed[chrom] = []
                values = np.array([start, end]).astype(int)
                values[1] += 1 # include upper value in list ([lower, upper])
                bed[chrom].append(values) 

            except StopIteration:
                break

    # sort the starts on each chromosome
    for chrom in bed.keys():
        X = np.vstack(bed[chrom])
        idx = np.argsort(X[:,0])
        bed[chrom] = X[idx]
        del X

    # have the lengths handy to avoid computing it on every iteration later
    bed_lens = {k:len(v) for k,v in bed.items()}
    total_bed_lens = np.sum([len(v) for k,v in bed.items()])

    return bed, bed_lens, total_bed_lens


if __name__=='__main__': main()
