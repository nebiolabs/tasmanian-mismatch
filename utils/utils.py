import numpy as np
from sam_reads import reads
import logging 
import sys,os
sys.path.append(os.path.abspath('../'))
from intersections import logFileName


logging.basicConfig(filename = logFileName,
                    format = '%(asctime)s %(message)s',
                    filemode = 'w')

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def read_bed(bedfile):
    '''
        read bed file
        INPUT: bedfile path <name>
        OUTPUT: 1. bed dictionary (chromosome=keys, [#start, #end]=values (np.array. each chrom sorted on #starts) 
                2. dictionary of lengths of arrays for each chromosome
                3. total length of bed elements.
    '''
    bed = {}
    bed_others = {}   # here I keep other fields than the coordinates
    with open(bedfile) as f:
        while True:
            try:
                #chrom, start, end = next(f).strip().split('\t')[:3]
                # Now we include information of repeats in columns>3
                chrom, start, end, strand, rep_name, rep_class, rep_family = next(f).strip().split('\t')[:7]

                if chrom not in bed:
                    bed[chrom] = []
                    bed_others[chrom] = []

                values = np.array([start, end]).astype(int)
                values[1] += 1 # include upper value in list ([lower, upper])

                bed[chrom].append(values)
                bed_others[chrom].append([strand, rep_name, rep_class, rep_family])

            except StopIteration:
                break

    # sort the starts on each chromosome
    for chrom in bed.keys():
        X = np.vstack(bed[chrom])
        idx = np.argsort(X[:,0])
        bed[chrom] = X[idx]
        bed_others[chrom] = [bed_others[chrom][n] for n in idx]
        del X

    # have the lengths handy to avoid computing it on every iteration later
    bed_lens = {k:len(v) for k,v in bed.items()}
    total_bed_lens = np.sum([len(v) for k,v in bed.items()])

    return bed, bed_lens, total_bed_lens, bed_others


def assign_category(read_1, read_2):
    '''
        INPUT: read1 and read2 are reads (sam_reads)
            7 possible categories in sam_reads. The value is initialized as None and 
            both reads 1&2 share the same category 
        OUTPUT: category of both reads
    '''
    # define read1 and read2 as first and second in pair
    if read_1.flag in [99, 83]:
        read1 = read_1
        read2 = read_2
    else:
        read1 = read_2
        read2 = read_1

    # case 1 was unrelated
    if read1.category == 1 and read2.category==1: 
        return 1

    elif read2.category == 1: # it's pairing reads case 2 (both unrelated were discared above)
        try:
            category = '2' + subcategory(read1)
        except Exception as e:
            logging.error('category positions is None and should be a number in read with \
                           id={} and error={}     {}'.format(read1._id, str(e), read1.category_positions))
            return None

    elif read1.category == 1: # it's case 3
        try:
            category = '3' + subcategory(read2)
        except Exception as e:
            logging.error('category positions is None and should be a number in read with \
                           id={} and error={}'.format(read2._id, str(e)))
            return None

    elif read1.category != 1 and read2.category !=1: # it's cases 4 to 7 --> NOW 4 or 5 (e.g. 5 is now 4dd)
        try:
            
            if read1.bed_id == read2.bed_id: # both reads intercept same repeat
                category = '4' + subcategory(read1) + subcategory(read2)

            else: # both reads intercept different repeat
                category = '5' + subcategory(read1) + subcategory(read2)
        except Exception as e:
            logging.error('category positions (1 or 2) is None and should be a number in read with \
                          id={} or {} and error={}'.format(read1._id, read2._id, str(e)))
            return None

    return category


# function to be used inside assign_category
def subcategory(read):
    ''' define this outside of assign_category because is called only once inside the function
        and closures are fast when local function is called multiple times, otherwise slower
    '''
    ab = np.array(read.category_positions) > 0

    if ab[0] and not ab[1]:         return 'a'
    elif not ab[0] and not ab[1]:   return 'b'
    elif ab[0] and ab[1]:           return 'c'
    elif not ab[0] and ab[1]:       return 'd'
