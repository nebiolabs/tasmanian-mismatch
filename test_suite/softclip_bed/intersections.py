#!/bin/python
'''
    input sam file from stdin
    bedfile as one of the arguments (chromosome, start, end - delimiter=tab)

'''
import os, sys
import numpy as np
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent) + '/utils/')
from sam_reads import reads
from utils import *

# Define outputs: 1- SAM lines including ONLY the intersecting part
#                 2- SAM lines including ONLY the NON intersecting part
#                 3- Statistics about these intersections

''' These are replaced with relateds
intersections = []  # part of reads in common between repeats bed regions and bam read
differences = []    # part of reads NOT in commin between repeats bed regions and bam read
'''
statistics = []     # length of the intersection is kept 
unrelateds = []     # reads without ANY intersection with bed regions are written in another file
relateds = []       # reads with WHICHEVER intersection with bed regions

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

    # define a table of flags of proper paired reads
    flags = {
        99: 'first fwd',
        147:'second fwd',
        83: 'first rev',
        163:'second rev'
    }

    # bufer to keep reads until the pair is found and they can be analyzed together.
    buffer = {}
    n_tester = 0
    # read bam file from stdin and 
    for line in sys.stdin:

        n_tester +=1

        line = line.strip('\n') # Don't repeat same opreation over and over on each iteration.

        # if finished looping over the bed file, don't read whole bam file.
        if finito:
            unrelateds.append(line)
            continue # still collect the unrelateds

        if len(line) < 50: break # avoid potential empty lines at the end.

        # instantiate object "current_read"
        _id, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred = line.split('\t')[:11]
        current_read = reads(_id, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred)
        paired_read = None # assume there is no paired_read yet
        skip_read = False

        # if chromosomoe not in bed, read is not intereseting also other poor situations added here
        if current_read.chrom not in bed: # or mapq<20 or np.abs(tlen)<100 or np.abs(tlen)>500: 
            unrelateds.append(line)
            continue

        # only consider uniquely mapped and proper pair
        if current_read.flag not in [163,99,147,83]: continue

        # Check-point ====================================================================================
        '''
            If read is first on the pair, save it to memory (buffer) until the paired read is found.
            If read is second on the pair, we already have the first in the memory, in a list of first
            reads (sometimes they are 100s of lines appart) and we compute all what follows for both reads.
            99: first -> I prefer not to rely on this. Exceptions could be many more than I think in 
            different samples. 
        '''
        if current_read._id in buffer:
            paired_read = buffer.pop(current_read._id) #, "None")   error is better than None here
        
        # sanity check
        if paired_read != None:
            if current_read.chrom != paired_read.chrom:
                print('paired read and read from different chromosomes?')
                continue

        # Bam Block ====================================================================================== 
        ''' In THIS block bam and bed are updated to same genomic region if there are gaps between the two. '''

        # new chromosome -> zero bed_index!
        if current_read.chrom != last_chrom:
            skip_chrom=False
            bed_index=0
            print('still have {} reads in buffer for chromosome {}. Thrown to trash  --  {}'.format(len(buffer), last_chrom, n_tester))
            buffer = {}
            if total_beds >= total_bed_lens: # Have we finished reading the bed file?
                finito = True
                break

        ### if read is ahead of bed fragment, advance index
        if not skip_chrom: # don't get in here is no more bed fragments for that chrome
            while current_read.start >= bed[chrom][bed_index,1]:
                if bed_index == bed_lens[chrom]-1: # If here, bed file finished for this chrome
                    skip_chrom = True
                    break
                bed_index +=1
                total_beds+=1
        
        ### continue until a read intersects the bed fragment or in case there are no more regions for this chrom in bed
        if skip_chrom or current_read.end <= bed[chrom][bed_index,0]:
            current_read.category = 1
            current_read.bed_id = bed_index

            if paired_read != None:     # paired already in buffer
                if paired_read.category == 1:
                    unrelateds.append(paired_read.print(sequence='original'))
                    unrelateds.append(current_read.print(sequence='original'))
                    last_chrom = current_read.chrom # update last_chrom since we don't reach the bottom
                    continue
            else:
                buffer[current_read._id] = current_read

        # Block for intersecting reads ==================================================================== 
        if current_read.category == None: # I have not assigned 1 to it! -> it intercepts a bed region!
            current_read.bed_id = bed_index

            # capture the intersection into a flag called 
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
            a0 = bed[chrom][bed_index][0] - current_read.start
            b0 = current_read.end - bed[chrom][bed_index][1]
            a = a0 if a0>=0 else 0
            b = -b0 if b0>0 else current_read.seq_len
            
            # add feature to read to classify it later
            current_read.category_positions = np.array([a,b]) # negative if contained else positive.
            
            # mask intersections with "N"    
            intersect_size = len(current_read.seq[a:b])  
            current_read.masked_seq = current_read.seq[:a] + ''.join(['N'] * intersect_size) + current_read.seq[b:]
            current_read.intersect_seq = ''.join(['N'] * a) + current_read.seq[a:b] + ''.join(['N'] * b0)

            # cigar is to correlate clips with intersections. It might happen that most 
            # intersections are not correctly mapped and therefore, are softcliped
            current_read.expand_cigar()     # previously called CIGAR
            intersect_size=None
            '''
            # mask intersections with "N"
            if a>0 and b>0:
                intersect_size = current_read.seq_len-b-a
                current_read.masked_seq = current_read.seq[:a] + ''.join(['N']*intersect_size) + current_read.seq[-b:]
                current_read.intersect_seq = ''.join(['N']*a) + current_read.seq[a:-b] + ''.join(['N']*b)
            elif a>0 and b<0:
                intersect_size = current_read.seq_len-a
                current_read.masked_seq = current_read.seq[:a] + ''.join(['N']*(current_read.seq_len-a))
                current_read.intersect_seq = ''.join(['N']*a) + current_read.seq[a:]
            elif b>0 and a<0:
                intersect_size = current_read.seq_len - b
                current_read.masked_seq = ''.join(['N']*intersect_size) + current_read.seq[intersect_size:]
                current_read.intersect_seq = current_read.seq[:-b] + ''.join(['N']*b)
            else:
                intersect_size = current_read.seq_len
                current_read.masked = ''.join(['N'] * current_read.seq_len)
                current_read.intersect_seq = current_read.seq
            '''
        # read_category will be included as a flag at the end of each read in bam file.
        # we nead both paired reads to calculate this flag.
        if paired_read == None:
            buffer[current_read._id] = current_read
            last_chrom = current_read.chrom   # update last_chrom since we don't reach the bottom
            # AT THIS POINT, current_read.category is either 1 or None.
            continue
        else:
            # is current or paired read1? THe other is read2
            read1, read2 = [current_read, paired_read] if current_read.flag in [99,83] else [paired_read, current_read]
            current_read.category = assign_category(read1, read2)
            paired_read.category = current_read.category # category is shared by the paired reads

            print(current_read.category, current_read.category_positions, 
                  paired_read.category, paired_read.category_positions, 
                  current_read.seq, current_read.intersect_seq, current_read.masked_seq, 
                  paired_read.seq, paired_read.intersect_seq, paired_read.masked_seq)
            
            #if current_read.category != 1:
            #    relateds.append(current_read.print('masked'))
            #else:
            #    pass
        
        # final desicion: print sam file with reads as they are in the original file but
        # add a flag with the category assigned to the read in class reads.
        # write SAM reads containing intersections or non-intersecting sequences
        #### flag, mapq, tlen, start = str(flag), str(mapq), str(tlen), str(start) #+1) Now there is no need and if this works --> DELETE!

        '''
        if diff_seq != None:
            differences.append('\t'.join([_1, flag, chrom, start, mapq, cigar, _2, _3, tlen, diff_seq, phred]))
        if diff_seq2 != None:
            differences.append('\t'.join([_1, flag, chrom, start, mapq, cigar, _2, _3, tlen, diff_seq2, phred]))
        
        intersections.append('\t'.join([_1, flag, chrom, start, mapq, cigar, _2, _3, tlen, inters_seq, phred]))
        
        #if int(intersect_size) > 76: print(seq, bed[chrom][bed_index], flag, start, end, tlen, intersect_size, a, b)
        if len(inters_seq)<1: print(seq, bed[chrom][bed_index], a, b, inters_seq, start, end, chrom)

        statistics.append(str(intersect_size))
        '''

        #print('testing if reached this point with {}'.format(n_tester))

        last_chrom = current_read.chrom

    # save files:
    with open('relateds.sam', 'w') as f:
        f.write('\n'.join(relateds))
    '''
    for data, filename in zip([intersections, differences, statistics, unrelateds],
                              ['intersections.sam','differences.sam','statistics.txt', 'unrelateds.sam']
                             ):
        with open(out_prefix + filename,'w') as f:
            f.write('\n'.join(data))
    '''



if __name__=='__main__': main()
