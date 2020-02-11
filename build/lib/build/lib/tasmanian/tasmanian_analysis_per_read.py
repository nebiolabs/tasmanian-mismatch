#!/bin/python

# TODO: add context nucleotides? 

HELP = """

    required:
    --------
    -r|--reference-fasta

    optional:
    --------
    -u|--unmask-genome (convert masked bases to upper case and include them in the calculations - default=False)
    -q|--base-quality (default=20)
    -f|--filter-indel (exclude reads with indels default=False)
    -l|--filter-length (include only reads with x,y range of lengths, default=0, 76)
    -s|--soft-clip-bypass (Decide when softclipped base is correct(0). Don\'t use these bases(1). Force use them(2).  default=0)
    -m|--mapping-quality (minimum allowed mapping quality -defailt=0)
    -h|--help
    -g|--fragment-length (use fragments withi these lengths ONLY)

"""

import sys, os, re, time
import numpy as np
from itertools import product
import uuid
randLogName = str(uuid.uuid4())
print('log filename: errors_'+randLogName+'.log')

logging = []
SKIP_READS = {
    'indel':0,
    'mapquality':0,
    'fragLength':0,
    'softclips':0,
    'readlength':0
}

# set default parameters
_UNMASK_GENOME=False  #don't unmask the genome, aka don't use lower letter from genome but rather keep them out in the error.log file'
READ_LENGTH=76
MinMapQuality = 0
PHRED = 20 + 33
SOFTCLIP_BYPASS=0
SKIP_INDEL=False
MIN_LENGTH=0
MAX_LENGTH=76
TLEN=np.array([0,10000])

# if there are arguments get them
for n,i in enumerate(sys.argv):
    if i in ['-r','--reference-fasta']:
        ref_file = sys.argv[n+1]
        print('using {} as reference genome'.format(ref_file))
    if i in ['-u', '--unmask-genome']:
        _UNMASK_GENOME = True 
        print('unmasking genome. You are potentially including repetitive regions...')
    if i in ['-q', '--base-quality']:
        PHRED = int(sys.argv[n+1])+ 33
        print('base queality set to {}'.format(PHRED))
    if i in ['-f', '--filter-indel']: 
        SKIP_INDEL=True
        print('Skipping all reads with INDELS')
    if i in ['-l', '--filter-length']: # min and max separated by a comma
        MIN_LENGTH, MAX_LENGTH = np.array(sys.argv[n+1].strip('\n').split(','), dtype=np.uint16)
        print('range of fragment lengths allowed is {}-{}'.format(MIN_LENGTH, MAX_LENGTH))
    if i in ['-s', '--soft-clip-bypass']:
        SOFTCLIP_BYPASS=int(sys.argv[n+1])
        print('softclip bypass set to {}'.format(SOFTCLIP_BYPASS))
    if i in ['-m', '--mapping-quality']:
        MinMapQuality = int(sys.argv[n+1])
        print('Filterung out reads with Mapping quality < {}'.format(MinMapQuality))
    if i in ['-g','--fragment-length']:
        TLEN = np.array(sys.argv[n+1].split(',')).astype(int)
        print('Only reads comming from fragments of lengths {}-{} are considered here'.format(TLEN[0], TLEN[1]))
    if i in ['-h', '--help']:
        exit(HELP)
    if i in ['--counts_per_read']:
        counts_filename = sys.argv[n+1] + '.npy'
        print('printing counts per read to {}'.format(counts_filename))

if len(sys.argv)==1: exit('\n-h|--help\n') # there should be at least one argument = '--reference-genome'

#if os.path.isfile(outputFile): os.remove(outputFile) # avoid clashes and remove file.
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
        logging.append('{} and {} have different lengths'.format(s1,s2))
        return

    if np.sum([1 for i in range(l) if s1[i]==s2[i]]) >= fraction_correct*l:
        return False
    else:
        return True

                                                #########################
                                                ## LOAD REFERENCE DATA ##
                                                #########################
t0 = time.time()
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

t1 = time.time()
print('read reference in {} minutes'.format((t1-t0)/60.))


                                                #########################
                                                ## LOAD ALIGNMENT DATA ##
                                                #########################

# define hash table to collect artifacts
errors = {}
for read in [1,2]:
    errors[read] = {}
    for pos in np.arange(READ_LENGTH):
        errors[read][pos] = {}
        for ref in ['A','C','G','T']:
            errors[read][pos][ref] = {}
            for alt in ['A','C','G','T']:
                errors[read][pos][ref][alt] = 0
    

# array of artifacts per read
artifactsANDposition_per_read = [] 


# read bam from "samtools view" and pipe 
for line in sys.stdin:

    artifacts_in_this_read = [] # each read will report a number of artifacts to "artifacts_per_read"

    _, flag, chrom, start, mapq, cigar, _, _, tlen, seq, phred = line.split('\t')[:11]
    
    skip_read = False
    seq_len = len(seq)
    mapq = int(mapq)
    flag = int(flag)
    start = int(start)-1
    end = start + seq_len
    tlen = int(tlen)

    # reasons to skip reads
    # =====================
    # problems with CIGAR (* is data unavailable)
    # unwanted read length
    # Bad mapping quality
    # unrecognized chromosome

    if cigar=="*" or mapq==255 or chrom not in reference: 
        continue  
    elif SKIP_INDEL and ("I" in cigar or "D" in cigar):
        SKIP_READS['indel']+=1      
    elif seq_len < MIN_LENGTH or seq_len > MAX_LENGTH:
        SKIP_READS['readlength']+=1
    elif mapq<MinMapQuality:
        SKIP_READS['mapquality']+=1
    elif tlen<TLEN[0] or tlen>TLEN[1]:
        SKIP_READS['fragLength']+=1
        continue

    # INSTEAD OF EXPAND THE CIGAR, I PARSED MORE EFFICIENTLY THE CONTENT INTO A NUMPY INDEX
    # if indels were not excluded, these are still different from the genome and have to be 
    # removed from the sequence and CIGAR.
    # I first include the "S" and then check if these are garbage, shifting the start or end in seq and if indeed are garbage
    # Remember that bwa mem (and possibly bowtie2) indicate start as the first "M" in the cigar, after the "S"!!
    # deal on the fly with softclips and use binary indexes for reference and sequence to deal with indels at the end. 
    seq_idx, ref_idx = np.ones(len(seq)), np.ones(len(seq))
    last, position = 0, 0

    for current, i in enumerate(cigar):
        if not i.isnumeric():
            number = int(cigar[last:current])

            if i=='S':
                if position==0: # beginning
                    start -= number # This already covers SOFTCLIP_BYPASS==2
                    end -= number
                    if SOFTCLIP_BYPASS==0:
                        s1, s2 = seq[:number], reference[chrom][start:start+number]
                        if simple_deltas_is_this_garbage(s1,s2):
                            seq = 'N' * number + seq[number:]
                    elif SOFTCLIP_BYPASS==1:
                        seq = 'N' * number + seq[number:]
                        SKIP_READS['softclips']+=1 
                else: # end
                    if SOFTCLIP_BYPASS==0:
                        s1, s2 = seq[position:position+number], reference[chrom][start+position:start+position+number]
                        if simple_deltas_is_this_garbage(s1,s2):
                            seq = seq[:position] + 'N' * number
                    elif SOFTCLIP_BYPASS==1:
                        seq = seq[:position] + 'N' * number
                        SKIP_READS['softclips']+=1 

            elif i=='I':
                seq_idx[position:position+number]=0
            
            elif i=='D':
                ref_idx[position:position+number]=0

            position += number
            last = current+1
    
    ref = reference[chrom][start:end]
    l = int(np.sort([np.sum(seq_idx), np.sum(ref_idx)])[0])
    try:
        seq = [seq[i] for i in range(len(seq_idx)) if seq_idx[i]==1][:l]
        ref = [ref[i] for i in range(len(ref_idx)) if ref_idx[i]==1][:l]
    except Exception as e:
        logging.append('error: {} occurred while fixing for different lengths of seq and ref'.format(str(e)))

    # if bin(int(flag))[2:][-5]=='1' or make it easy for now...
    if flag==99: 
        strand='fwd'; read=1
    elif flag==163: 
        strand='fwd'; read=2
    elif flag==83:
        strand='rev'; read=1
    elif flag==147: 
        strand='rev'; read=2
    else: continue


    if _UNMASK_GENOME: 
        ref = ''.join(ref).upper() # re-think before doing this.

    if len(seq) != len(ref):
        logging.append('seq ={} and ref={} have different lengths'.format(seq,ref))
        continue

    for pos,base in enumerate(seq):
        if pos > len(ref) or pos > len(phred): 
            logging.append('ERROR processing read {}, with ref={}, phred={} and seq={}'.format(_,ref,phred,seq))
            continue

        if base=='N' or ref[pos]=='N' or ord(phred[pos]) < PHRED: # or (CIGAR[pos] not in ['M','X','=']):
            continue
        else:
            read_pos = [pos if strand=='fwd' else seq_len-pos-1][0]
            try:
                # evaluate this..
                if strand=='fwd': 
                    ref_pos = ref[pos]
                    Base = base
                elif strand=='rev': 
                    ref_pos = revcomp(ref[pos])
                    Base = revcomp(base)
            
                errors[read][read_pos][ref_pos][Base] += 1 

                if ref_pos != Base: 
                    artifacts_in_this_read.append(read_pos)

            except Exception as e:
                logging.append('error:{} in chr:{}, position:{}, read:{}, base:{}, seq:{}, start:{} and ref_pos:{}'.format(\
                                                             str(e), chrom, pos, read, base, ''.join(seq), start, ref[pos]))

    # fill results table.
    if artifacts_in_this_read == []:
        artifacts_in_this_read = [-1]
    artifactsANDposition_per_read.append(artifacts_in_this_read)

# serialize counts array into file
np.save(counts_filename, np.array(artifactsANDposition_per_read))


t2=time.time()
print(t2-t0, "seconds")
print('reads discarded')
print('===============')
for k,v in SKIP_READS.items():
    print('{:<15} {}'.format(k,v))  

# print header
sys.stdout.write('read,position,a_a,a_t,a_c,a_g,t_a,t_t,t_c,t_g,c_a,c_t,c_c,c_g,g_a,g_t,g_c,g_g\n')
prnt = []
for read in [1,2]:
    for pos in np.arange(76): 

        prnt.append(str(read) + ',' + str(pos+1) + ',' + ','.join([str(errors[read][pos][i][j]) for i,j in \
               zip(['A','A','A','A','T','T','T','T','C','C','C','C','G','G','G','G'], ['A','T','C','G']*4)]))

sys.stdout.write('\n'.join(prnt))

f = open('errors_'+randLogName+'.log','w',)
f.write('\n'.join(logging))
f.close()
