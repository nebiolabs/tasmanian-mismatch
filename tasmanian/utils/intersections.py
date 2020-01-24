#!/bin/python
'''
    input sam file from stdin
    bedfile as one of the arguments (chromosome, start, end - delimiter=tab)

'''
import os, sys, gzip
import numpy as np
#from pathlib import Path
#sys.path.append(str(Path(__file__).parent.parent.parent) + '/utils/')
sys.path.append(sys.path[0] + '/utils')
from sam_reads import reads
from utils import *

HELP = '''
\t\tsamtools view <bam_file> | python -b <bed_file/bedGraph> -o <output.table> 
'''

# initialize lists to contain reads (sam) and statistics (e.g. length of intersections)
statistics = []      
sam_masked = []
sam_intersections = []

# load global arguments
out_prefix = ''

for n,i in enumerate(sys.argv):
    if i in ['--bed','-b','--bed-file']:
        bedfile = sys.argv[n+1]

    if i in ['--output', '-o']:
        out_prefix = sys.argv[n+1] + "."

    if i in ['-h','--help']:
        print(HELP)
        sys.exit(1)

# define the name of the logging file for debugging
logFileName = out_prefix + 'log'


def main():

    # load bed_file
    try:
        bed, bed_lens, total_bed_lens, bed_other_info = read_bed(bedfile)
        logging.info('bedfile {} was succesfully read'.format(bedfile))
    except Exception as e:
        logging.error('{} happened in excecution of read_bed in main()'.format(str(e)))
        exit("there was a problem reading {}. Make sure is tab delimited and all columns \
             and rows are correct.")

    # to avoid going over the entire chromosome on each read,
    # keep an updated index of the bedfile to start from
    total_beds = 0 
    last_chrom = 'chr0'
    finito=False # know when we are finished with the bed file

    # define a table of flags of proper paired reads
    proper_flags = {
        99: 'first fwd',
        147:'second fwd',
        83: 'first rev',
        163:'second rev'
    }

    # buffer to keep reads until the pair is found, so the paired reads can be analyzed together.
    buffer = {}
    n_tester = 0

    # read bam file from stdin and 
    for line in sys.stdin:
        n_tester +=1

        line = line.strip('\n') 

        if len(line) < 50:  # avoid potential empty lines at the end.
            logger.warning('line {} had less than 50 characters'.format(n_tester))
            continue 

        # instantiate object "current_read"
        try:
            _id, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred = line.split('\t')[:11]
            current_read = reads(_id, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred)

            # assume read is not paired yet
            paired_read = None # assume there is no paired_read yet

        except Exception as e:
            logging.error('read {} could\'t be loaded properly'.format(n_tester))
            continue            

        # only consider uniquely mapped and proper pair
        if current_read.flag not in proper_flags:
            continue

        # if chromosome not in bed, read is non overlapping
        if current_read.chrom not in bed:
            sam_masked.append(current_read.print('original'))
            continue

        # new chromosome -> zero bed_index!
        if current_read.chrom != last_chrom:
            #logging.critical('This should ONLY happen ONCE!!! {} != {}'.format(current_read.chrom, last_chrom))
            logging.info('still have {} reads in buffer for chromosome {}. Thrown to trash  --  {}'.format(len(buffer), last_chrom, n_tester))
            skip_chrom=False                    # in case this was on for the previous chromosome
            bed_index=0                         # new chromosome, new bed_index
            buffer = {}                         # restart buffer for the new chromosome
            if total_beds >= total_bed_lens:    # Have we finished reading the bed file?
                finito = True
            last_chrom = current_read.chrom

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
        if paired_read != None: # TESTED
            if current_read.chrom != paired_read.chrom:
                logging.warning('current_read and paired_read have different chromosomes = chimeras')
                continue

        # Bam Block ====================================================================================== 
        ''' In THIS block bam and bed are updated to same genomic region if there are gaps between the two. '''

        # if finished looping over the bed file or bed is ahead of bam
        if finito or skip_chrom or current_read.end <= bed[chrom][bed_index,0]:
            current_read.category = 1


        # if read is ahead of bed fragment, advance index
        else: # don't get in here if no more bed fragments for that chrom or bed is finished
            while current_read.start >= bed[chrom][bed_index,1]:

                if bed_index == bed_lens[chrom]-1: # If here, bed file finished for this chrome
                    skip_chrom = True
                    break
                else:
                    bed_index +=1
                    total_beds+=1

        # If: 1.bam was lower than bed  2.we skipped bed with no bam coverage 3.we got to a bed that's after 
        # the bam region 4. Even though current_read.category is None, no "ab" (see the following section) 
        # will be assigned. Hence, we can already assign a category=1 here.
        if current_read.end <= bed[chrom][bed_index,0]:
            current_read.category=1
        
        # We can already assign bed_id       
        current_read.bed_id = bed_index
        current_read.bed_extra_info = bed_other_info[chrom][bed_index]


        # Block for intersecting reads ==================================================================== 
        if current_read.category == None: # I have not assigned 1 to it! -> it intercepts a bed region!

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
            current_read.category_positions = [a,b] 
            
            # mask intersections with "N"
            intersect_size = len(current_read.seq[a:b])  
            current_read.masked_seq = current_read.seq[:a] + ''.join(['N'] * intersect_size) + current_read.seq[b:]
            current_read.intersect_seq = ''.join(['N'] * a) + current_read.seq[a:b] + ''.join(['N'] * b0)

            # cigar is to correlate clips with intersections. It might happen that most 
            # intersections are not correctly mapped and therefore, are softcliped
            current_read.expand_cigar()     # previously called CIGAR
            intersect_size=None

        # Current_read.category is already assigned ========================================================
        # read_category will be included as a flag at the end of each read in bam file.
        # we nead both paired reads to calculate this flag.
        if paired_read == None: 
            buffer[current_read._id] = current_read

        else:  
            current_read.category = assign_category(current_read, paired_read)
            paired_read.category = current_read.category # category is shared by the paired reads
            sam_masked.append(paired_read.print('masked'))   # if masked_seq==None goes for 'original' by default.
            sam_masked.append(current_read.print('masked'))  # same here!
    
            for READ in [paired_read, current_read]:
                if READ.intersect_seq != None:
                    if READ.masked_seq == None: logging.critical('This should not be happening, masked is None and intersect is not...?')
                    sam_intersections.append(READ.print('intersect'))


    # save files:
    for mtx, fle in zip([sam_masked, sam_intersections], ['masked','intersections']):
        with gzip.open(out_prefix + fle + '.sam.gz', 'wb') as f:
            f.write('\n'.join(mtx).encode())
    
        print(fle, ' written!')
        del mtx
        del f

    return "Done!"



if __name__=='__main__': main()