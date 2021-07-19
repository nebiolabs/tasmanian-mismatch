import numpy as np
import logging 
import sys, os, re

try:
    from tasmanian.utils.sam_reads import reads
except Exception:
    # Either tests or base_dir, it's downstream of ../tasmanian/tasmanian/
    p = os.path.abspath(os.path.dirname(__file__))
    #p = re.search("(.*tasmanian/tasmanian/).*",p).group(1)
    p_start = [i for i in re.finditer('/tasmanian',p)][-1].end()
    p = p[:p_start]
    utils_path = p + '/utils'
    sys.path = [utils_path] + sys.path

    #p = re.search("(.*tasmanian/tasmanian/).*",p).group(1)
    #utils_path = p + 'utils'
    #sys.path = [utils_path] + sys.path
    from sam_reads import reads 

#sys.path.append(os.path.abspath('../'))


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
    wrong_syntax_lines = 0

    with open(bedfile) as f:
        while True:
            try:
                line =  next(f).strip().split('\t')
                num_fields = len(line)

                # Now we include information of repeats in columns>3
                if num_fields >=7:
                    chrom, start, end, strand, rep_name, rep_class, rep_family = line[:7]

                elif num_fields == 3: 
                    chrom, start, end = line

                    # worth the overhead if we use this data in intersections.py when exists
                    strand, rep_name, rep_class, rep_family = '','','',''                  
                else:
                    wrong_syntax_lines+=1
                    if wrong_syntax_lines>50:
                        sys.stderr.write('bed_file is not formatted as expected... see help menu')
                        exit(1)

                if not start.isnumeric() or not end.isnumeric():  # this migh be the header or empty lines
                    continue

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


def load_reference(ref_file):
    '''
        open reference file into dictionary of chromosomes and their sequences
        INPUT: name of reference file
        OUTPUT: reference dictionary
    '''
    reference = {}
    with open(ref_file) as f:
        while True:
            try:
                line = next(f).strip()

                if len(line)<2: continue

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
    

# Alert in case the inputs are incorrect
class cols:
    bold = '\x1b[0;30;47m'
    normal = '\x1b[0m'

def revcomp(base):
    d={'A':'T','T':'A','C':'G','G':'C',
       'a':'t','t':'a','c':'g','g':'c'}
    return d[base]

def simple_deltas_is_this_garbage(s1,s2,fraction_correct=0.66): # since strings should have same length, this is simple
    l=len(s1)
    if len(s2) != l:
        logging.error('{} and {} have different lengths'.format(s1,s2))
        return

    if np.sum([1 for i in range(l) if s1[i]==s2[i]]) >= fraction_correct*l:
        return False
    else:
        return True


def init_artifacts_table(READ_LENGTH):
    '''
        Initialize a dictionary to store artifact results
        INPUT: READ_LENGTH (the length of the reads)
        OUTPUT: errors table of all possible missmatches
    '''
    errors = {}
    for read in [1,2]:
        errors[read] = {}
        for pos in np.arange(READ_LENGTH):
            errors[read][pos] = {}
            for ref in ['A','C','G','T']:
                errors[read][pos][ref] = {}
                for alt in ['A','C','G','T']:
                    errors[read][pos][ref][alt] = 0
    return errors


def trim_table(table, mode_length):
    '''
        table is initialized with a read_length that is most ofter
        way longer than the reads actually are. Here we trim the 
        excess length.
        INPUT: artifact_table
        OUTPU: artifact_table
    '''
    trimmed_table = {}
    
    for read in [1,2]:
        trimmed_table[read] = {}
        for pos in range(mode_length):
            trimmed_table[read][pos] = {}
            for ref in ['A','C','G','T']:
                trimmed_table[read][pos][ref] = {}
                for alt in ['A','C','G','T']:
                    trimmed_table[read][pos][ref][alt] = table[read][pos][ref][alt]

    return trimmed_table


# define a table of flags of proper paired reads
proper_flags = {
    99: 'first fwd',
    147:'second fwd',
    83: 'first rev',
    163:'second rev'
}

def initialize_PFM(flanking_n=5):
    '''
    Initialize a position frequency matrix
    with flanking_n * 2 columns and 4 rows. (not including middle-base)
    Matrix looks loke this:
    A       -n ... 0 ... n   -> matrix[0]
    C        .     .     .
    T        .     .     .
    G        .     .     .
    N        .     .     .
    '''
    return np.zeros(shape=(flanking_n*2+1, 5))

def fill_PFM(ref_seq, matrix):
    '''
    The matrix should represent the mismatch. This should be checked
    before calling this function.
    We could check the flanking_n here, but we can speed up a tincy bit 
    having the value instead of calculating it.
    ALSO the seqeuence from reference IS PROVIDED, not the ref_position.

    E.g. ref_seq = ATTGCTTAG -> flanking_n=4 and base=C
    '''
    #matrix2 = matrix.copy()
    #assert len(ref_seq) <= matrix.shape[0]

    loc_m = {'A':0, 'C':1, 'T':2, 'G':3, 'N':4} # position (location) in the matrix
    for n,i in enumerate(ref_seq):
        matrix[n,loc_m[i]] +=1
    
    #return matrix2
    return

def pfm2ppm(matrix):
    return matrix / matrix.sum(axis=1).reshape(-1,1)

def ppm2pwm(ppm1, ppm2):
    '''
    ppm1 = measured distribution
    ppm2 = base distribution
    e.g. ppm1='c_t' then ppm2='c_c' 
    '''
    return np.log2(ppm1 / ppm2)

