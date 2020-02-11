#!/bin/python

# TODO: add context nucleotides? 

import sys, os, re, time
import numpy as np
import pandas as pd
from itertools import product
import uuid
#sys.path.append(os.path.abspath(os.path.dirname(__file__)) + '/utils/')
from tasmanian.utils.utils import revcomp, simple_deltas_is_this_garbage, init_artifacts_table, load_reference
from tasmanian.utils.plot import plot_html

###############################################################################
# In order to make the binary scripts work, make all these scripts modular.   # 
# Then, I can import the function from thee module into the scripts in bin/   #
###############################################################################

def analyze_artifacts(Input, Args):

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
        -o|--output-prefix (use this prefix for the output and logging files)
    """

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
    READ_LENGTH=200
    MinMapQuality = 0
    PHRED = 20 + 33
    SOFTCLIP_BYPASS=0
    SKIP_INDEL=False
    MIN_LENGTH=0
    MAX_LENGTH=200
    TLEN=np.array([0,10000])
    randLogName = str(uuid.uuid4())

    # if there are arguments get them
    for n,i in enumerate(Args):
        if i in ['-r','--reference-fasta']:
            ref_file = sys.argv[n+1]
            sys.stderr.write('using {} as reference genome\n'.format(ref_file))
        if i in ['-u', '--unmask-genome']:
            _UNMASK_GENOME = True 
            sys.stderr.write('unmasking genome. You are potentially including repetitive regions...\n')
        if i in ['-q', '--base-quality']:
            PHRED = int(sys.argv[n+1])+ 33
            sys.stderr.write('base queality set to {}\n'.format(PHRED))
        if i in ['-f', '--filter-indel']: 
            SKIP_INDEL=True
            sys.stderr.write('Skipping all reads with INDELS\n')
        if i in ['-l', '--filter-length']: # min and max separated by a comma
            MIN_LENGTH, MAX_LENGTH = np.array(sys.argv[n+1].strip('\n').split(','), dtype=np.uint16)
            sys.stderr.write('range of fragment lengths allowed is {}-{}\n'.format(MIN_LENGTH, MAX_LENGTH))
        if i in ['-s', '--soft-clip-bypass']:
            SOFTCLIP_BYPASS=int(sys.argv[n+1])
            sys.stderr.write('softclip bypass set to {}\n'.format(SOFTCLIP_BYPASS))
        if i in ['-m', '--mapping-quality']:
            MinMapQuality = int(sys.argv[n+1])
            sys.stderr.write('Filterung out reads with Mapping quality < {}\n'.format(MinMapQuality))
        if i in ['-g','--fragment-length']:
            TLEN = np.array(sys.argv[n+1].split(',')).astype(int)
            sys.stderr.write('Only reads comming from fragments of lengths {}-{} are considered here\n'.format(TLEN[0], TLEN[1]))
        if i in ['-h', '--help']:
            exit(HELP)
        if i in ['-o','--output-prefix']:
            randLogName = sys.argv[n+1]

    if len(sys.argv)==1: exit('\n-h|--help\n') # there should be at least one argument = '--reference-genome'

    #if os.path.isfile(outputFile): os.remove(outputFile) # avoid clashes and remove file.
    sys.stderr.write('log filename: errors_'+randLogName+'.log\n')

                                                    #########################
                                                    ## LOAD REFERENCE DATA ##
                                                    #########################
    t0 = time.time()
    reference = load_reference(ref_file)
    t1 = time.time()
    sys.stderr.write('read reference in {} minutes\n'.format((t1-t0)/60.))


                                                    #########################
                                                    ## LOAD ALIGNMENT DATA ##
                                                    #########################

    # define hash table to collect artifacts
    # table 1: all intersecting parts on reads that intersect bed file
    # table 2: all non-intersecting parts on reads that intersect bed file
    # table 3: all reads without interceptions.
    errors_intersection = init_artifacts_table(READ_LENGTH)
    errors_complement = init_artifacts_table(READ_LENGTH)
    errors_unrelated = init_artifacts_table(READ_LENGTH)
    
    # to check that there are tasmanian flags in sam file and if not, use unrelateds table only.
    bed_tag = False
    bed_tag_counter = 0
     
    # read bam from "samtools view" and pipe 
    for line in sys.stdin:
        
        # bypass header
        if line[0]=="@":
            continue

        _, flag, chrom, start, mapq, cigar, _, _, tlen, seq, phred = line.split('\t')[:11]

        # If sam data does not have the tasmanian tag, there is no intersection with bed and proceed to a single output table
        if bed_tag_counter<10 and bed_tag==False:
            try:
                tm_tag = np.array(re.search('tm:Z:([-.0-9]*;)', line).group(1).replace(';','').split('.')).astype(int)
                bed_tag = True
            except AttributeError:
                tm_tag = -1
            bed_tag_counter+=1

        elif bed_tag==True:
             tm_tag = np.array(re.search('tm:Z:([-.0-9]*;)', line).group(1).replace(';','').split('.')).astype(int)

        else:
            tm_tag = -1

        skip_read = False
        seq_len = len(seq)
        mapq = int(mapq)
        flag = int(flag)
        start = int(start)-1
        end = start + seq_len
        tlen = abs(int(tlen))

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
    
        #if strand=='rev':
        #    rev_seq = revcomp(seq)
        #    sys.stderr.write(seq, rev_seq)

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
                
                    #if pos == 0:
                    #    sys.stderr.write(seq[0], strand, ref_pos, Base)

                    if tm_tag[0] == -1:
                        assert base in ['A','C','G','T'], "{} should be upper case".format(base)
                        errors_unrelated[read][read_pos][ref_pos][Base] += 1
 
                    elif pos >= tm_tag[0] and pos < tm_tag[1]:
                        assert base in ['a','c','t','g'], "{} should be lower case".format(base) 
                        errors_intersection[read][read_pos][ref_pos][Base.upper()] += 1

                    else:
                        assert base in ['A','C','G','T'], "{} should be upper case".format(base)
                        errors_complement[read][read_pos][ref_pos][Base] += 1

                except Exception as e:
                    logging.append('error:{} in chr:{}, position:{}, read:{}, base:{}, seq:{}, start:{} and ref_pos:{}'.format(\
                                                                 str(e), chrom, pos, read, base, ''.join(seq), start, ref[pos]))

    #######################
    ## REPORTING SECTION ##
    #######################

    # Report performance and less relevant statistics
    t2=time.time()
    sys.stderr.write('tasmanian finished the analysis in {} seconds \n'.format(str(t2-t0)))
    sys.stderr.write('reads discarded\n')
    sys.stderr.write('===============\n')
    for k,v in SKIP_READS.items():
        sys.stderr.write('{:<15} {}\n'.format(k,v))  

    # print table header
    sys.stdout.write('read,position,' + ','.join \
        ([','.join([
            Sample+mut for mut in 'a_a,a_t,a_c,a_g,t_a,t_t,t_c,t_g,c_a,c_t,c_c,c_g,g_a,g_t,g_c,g_g'.split(',')
            ]) for Sample in ['I','C','N']
        ]) + '\n')

    # print rows
    prnt = []
    for read in [1,2]:
        for pos in np.arange(READ_LENGTH): 
            prnt.append(str(read) + ',' + str(pos+1) + ',' + ','.join([ \
                ','.join([str(Table[read][pos][i][j]) for i,j in zip(['A','A','A','A','T','T','T','T','C','C','C','C','G','G','G','G'], ['A','T','C','G']*4)])
            for Table in [errors_intersection, errors_complement, errors_unrelated]]))

    sys.stdout.write('\n'.join(prnt))

    # also save table as numpy arrays in dicts for ploting into the html report
    # we already printed the table as an array of strings. Just parse it into pandas and 
    # re-assign data-types
    col_names = ['read','position'] + ['_'.join([i[0].lower(), i[1].lower()]) for i in  
                    zip(['A','A','A','A','T','T','T','T','C','C','C','C','G','G','G','G'], ['A','T','C','G']*4)]

    df = pd.DataFrame([ i.split(',') for i in prnt]) 

    ids_intersection    = np.hstack([np.array([0,1]), np.arange(2,2+16)])
    ids_complement      = np.hstack([np.array([0,1]), np.arange(18,18+16)])
    ids_nonintersection =  np.hstack([np.array([0,1]), np.arange(34,34+16)]) 

    dfi = df.iloc[:, ids_intersection]
    dfc = df.iloc[:, ids_complement]
    dfn = df.iloc[:, ids_nonintersection]    

    table = {}
    for DF, dfName in zip([dfi, dfc, dfn],['intersection','complement','non_intersection']):
        DF.columns = col_names
        table[dfName] = DF.astype(int)
    
    return table



if __name__=='__main__':

    # logging in debug mode
    if '--debugging-mode' in sys.argv:
        logging.basicConfig(filename = logFileName,
                            format = '%(asctime)s %(message)s',
                            filemode = 'w')
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

    # run tasmanian
    table = analyze_artifacts(sys.stdin, sys.argv)

    # avoid overrighting report file before saving.
    report_filename, file_num = 'Tasmanian_artifact_report.html', 0
    while os.path.isfile(report_filename):
        file_num += 1
        report_filename = 'Tasmanian_artifact_report-' + str(file_num) + '.html'

    # create html template with plot      
    report_html = plot_html(table)
    
    with open(report_filename, 'w') as f:
        f.write(report_html)
    
    # Report errors to a logging file
    f = open('errors_'+randLogName+'.log','w',)
    f.write('\n'.join(logging))
    f.close()
