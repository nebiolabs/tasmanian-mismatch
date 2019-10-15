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

intersections = []
differences = []
statistics = []

# load arguments
for n,i in enumerate(sys.argv):
    if i=='--bed' or i=='-b' or i=='--bed-file':
        bedfile = sys.argv[n+1]


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

    # read bam file from stdin and 
    for line in sys.stdin:
        # if finished looping over the bed file, don't read whole bam file.
        if finito or len(line) < 50: break # also avoid potential empty lines at the end.

        _1, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred = line.split('\t')[:11]
        
        skip_read = False
        seq_len = len(seq)
        mapq = int(mapq)
        flag = int(flag)
        tlen = int(tlen)

        # if chromosomoe not in bed, read is not intereseting also other poor situations added here
        if chrom not in bed or mapq<30 or np.abs(tlen)<100 or np.abs(tlen)>500: continue

        # only consider uniquely mapped and proper pair
        if flag==163 or flag==99: 
            strand = "fwd"
            start = int(start)-1
            end = start + seq_len
        elif flag==147 or flag==83: 
            strand = "rev"
            end = int(_3)
            start = end - seq_len
        else: continue

        # new chromosome -> zero bed_index!
        if chrom != last_chrom:
            if total_beds >= total_bed_lens: # Have we finished reading the bed file?
                finito = True
                break
            bed_index=0

        # if read is ahead of bed fragment, advance index
        while start > bed[chrom][bed_index,1]:
            if bed_index >= bed_lens[chrom]-1: # CAREFULL, stop if finished reading bed file.
                break
            bed_index +=1
            total_beds+=1

        # continue until a read intersects the bed fragment
        if end < bed[chrom][bed_index,0]: 
            continue
        
        # capture the intersection
        # read      ------------------- 
        # bed        a  =========  b
        # read              -----------
        # read      --------
        # read           ------

        a = bed[chrom][bed_index][0] - start
        b = end - bed[chrom][bed_index][1]

        # need to consider including ONLY a true cigar in the output
        CIGAR = expand_cigar(cigar)
        intersect_size=None
    
        diff_seq, diff_cigar, diff_phred, diff_start = None, None, None, None
        diff_seq2, diff_cigar2, diff_phred2, diff_start2 = None, None, None, None
        inters_seq = None

        if a>0 and b>0:
            print(start, end, bed[chrom][bed_index], "class 1\n")       #########
            diff_seq, diff_cigar, diff_phred = seq[:a], cigar[:a], phred[:a]
            diff_start = start

            inters_seq, inters_cigar, inters_phred = seq[a:b], cigar[a:b], phred[a:b]
            inters_start = start + a
            intersect_size = seq_len-b-a
            
            diff_seq2, diff_cigar2, diff_phred2 = seq[-b:], cigar[-b:], phred[-b:]
            diff_start2 = start + intersect_size

        elif a>0 and b<=0:  # elif
            print(start, end, bed[chrom][bed_index], "class 2\n")       #########
            diff_seq, diff_cigar, diff_phred = seq[:a], cigar[:a], phred[:a]
            diff_start = start

            inters_seq, inters_cigar, inters_phred = seq[a:], cigar[a:], phred[a:]
            inters_start = start + a
            intersect_size = seq_len-a
            
        elif b>0 and a<=0:
            print(start, end, bed[chrom][bed_index], "class 3\n")       #########
            B = seq_len - b
            diff_seq, diff_cigar, diff_phred = seq[B:], cigar[B:], phred[B:]
            diff_start = start + B

            inters_seq, inters_cigar, inters_phred = seq[:B], cigar[:B], phred[:B]
            inters_start = start
            intersect_size = B
        
        else:
            print(start, end, bed[chrom][bed_index], "class 4?\n")      #########
            inters_seq, inters_cigar, inters_phred = seq, cigar, phred
            inters_start = start
            intersect_size = seq_len


        # write SAM reads containing intersections or non-intersecting sequences
        flag, mapq, tlen = str(flag), str(mapq), str(tlen)
        if diff_seq != None:
            differences.append('\t'.join([_1, flag, chrom, str(diff_start), mapq, diff_cigar, _2, _3, tlen, diff_seq, diff_phred]))
        if diff_seq2 != None:
            differences.append('\t'.join([_1, flag, chrom, str(diff_start2), mapq, diff_cigar2, _2, _3, tlen, diff_seq2, diff_phred2]))

        intersections.append('\t'.join([_1, flag, chrom, str(inters_start), mapq, inters_cigar, _2, _3, tlen, inters_seq, inters_phred]))
        
        #if int(intersect_size) > 76: print(seq, bed[chrom][bed_index], flag, start, end, tlen, intersect_size, a, b)
    
        if len(inters_seq)<1: print(seq, bed[chrom][bed_index], a, b, inters_seq, start, end)

        statistics.append(str(intersect_size))
        last_chrom = chrom

    # save files:
    for data, filename in zip([intersections, differences, statistics],['intersections.sam','differences.sam','statistics.sam']):
        with open(filename,'w') as f:
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

                bed[chrom].append(np.array([start, end]).astype(int))

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
    
