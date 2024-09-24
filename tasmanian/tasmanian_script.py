#!/usr/bin/env python

# TODO: add context nucleotides? 

import sys, os
project_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
sys.path.insert(0, project_root)
from tasmanian.utils.constants import constants

# sanity check before loading all the modules for faster execution
if len(sys.argv)==1 or '-h' in sys.argv or '--help' in sys.argv:
    exit(constants.HELP)

import os, re, time, pickle
import numpy as np
import pandas as pd
from itertools import product, permutations
import logging, uuid
from scipy.stats import mode

from tasmanian.utils.utils import revcomp, simple_deltas_is_this_garbage, init_artifacts_table, load_reference, trim_table
from tasmanian.utils.utils import fill_PFM, initialize_PFM, pfm2ppm, ppm2pwm
from tasmanian.utils.plot import plot_html

###############################################################################
# In order to make the binary scripts work, make all these scripts modular.   # 
# Then, I can import the function from thee module into the scripts in bin/   #
###############################################################################

def analyze_artifacts(Input, Args):
    globals().update(vars(constants))

    MinMapQuality = 20
    randLogName = str(uuid.uuid4())
    check_lengths = [constants.READ_LENGTH] # this is to rezise the results table from READ_LENGTHS to the mode of the lengths
    check_lengths_counter = 0
    confidence = 20
    debug = False
    picard = False
    flanking_n = False, 5

    # if there are arguments get them
    for n,i in enumerate(Args):
        if i in ['-r','--reference-fasta']:
            ref_file = sys.argv[n+1]
            sys.stderr.write('using {} as reference genome\n'.format(ref_file))
        if i in ['-u', '--unmask-genome']:
            constants._UNMASK_GENOME = True 
            sys.stderr.write('unmasking genome. You are potentially including repetitive regions...\n')
        if i in ['-q', '--base-quality']:
            constants.PHRED = int(sys.argv[n+1])+ 33
            sys.stderr.write('minimum base quality set to {}\n'.format(constants.PHRED))
        if i in ['-f', '--filter-indel']: 
            constants.SKIP_INDEL=True
            sys.stderr.write('Skipping all reads with INDELS\n')
        if i in ['-l', '--filter-length']: # min and max separated by a comma
            constants.MIN_LENGTH, constants.MAX_LENGTH = np.array(sys.argv[n+1].strip('\n').split(','), dtype=np.uint16)
            sys.stderr.write('range of fragment lengths allowed is {}-{}\n'.format(constants.MIN_LENGTH, constants.MAX_LENGTH))
        if i in ['-s', '--soft-clip-bypass']:
            constants.SOFTCLIP_BYPASS=int(sys.argv[n+1])
            sys.stderr.write('softclip bypass set to {}\n'.format(constants.SOFTCLIP_BYPASS))
        if i in ['-m', '--mapping-quality']:
            MinMapQuality = int(sys.argv[n+1])
            sys.stderr.write('Filtering out reads with Mapping quality < {}\n'.format(MinMapQuality))
        if i in ['-g','--fragment-length']:
            constants.TLEN = np.array(sys.argv[n+1].split(',')).astype(int)
            sys.stderr.write('Only reads comming from fragments of lengths {}-{} are considered here\n'.format(constants.TLEN[0], constants.TLEN[1]))
        if i in ['-o','--output-prefix']:
            randLogName = sys.argv[n+1]
        if i in ['-c', '--confidence']:
            confidence = int(sys.argv[n+1])
        if i in ['-d','--debug']:
            debug = True
        if i in ['-O','--ont']:
            constants.ONT = True
            constants.READ_LENGTH=100000
            constants.MAX_LENGTH=100000
            constants.TLEN=np.array([0,100000])
            check_lengths = constants.READ_LENGTH
        if i in ['-p','--picard-logic']:
            picard = True
        if i in ['-P','--include-pwm']:
            PWM = True
            if len(sys.argv)>n+1: 
                flanking_n = int(sys.argv[n+1]) if sys.argv[n+1].isnumeric() else flanking_n


    if debug:
        # if debugging create this logfile    
        logging.basicConfig(filename = randLogName + '.log',
                            format = '%(asctime)s %(message)s',
                            filemode = 'w')
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)


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
    errors_intersection = init_artifacts_table(constants.READ_LENGTH)
    errors_complement = init_artifacts_table(constants.READ_LENGTH)
    errors_unrelated = init_artifacts_table(constants.READ_LENGTH)

    # in thee tables the CONFIDENCE reads ONLY. These are reads with >= CONFIDENCE bases in the complement
    errors_intersectionB = init_artifacts_table(constants.READ_LENGTH)
    errors_complementB = init_artifacts_table(constants.READ_LENGTH)

    # initialize PFM (later on converted to PWM)    
    if constants.PWM:
        All_combinations = [
            ''.join(i) for i in permutations(['A','C','T','G'], 2)] +\
                 ['AA','CC','GG','TT']

        PFM = {i: initialize_PFM(flanking_n=flanking_n) for i in All_combinations}

    # to check that there are tasmanian flags in sam file and if not, use unrelateds table only.
    bed_tag = False
    bed_tag_counter = 0
     
    # read bam from "samtools view" and pipe 
    for line in sys.stdin:
        
        # bypass header
        if line[0]=="@":
            continue

        read_id, flag, chrom, start, mapq, cigar, _, _, tlen, seq, phred = line.split('\t')[:11]

        # If sam data does not have the tasmanian tag, there is no intersection with bed and proceed to a single output table
        if bed_tag_counter<10 and bed_tag==False:
            try:
                tm_tag = np.array(re.search('tm:Z:([-.0-9]*;)', line).group(1).replace(';','').split('.')).astype(int)
                bed_tag = True
            except AttributeError:
                tm_tag = [-1]

            bed_tag_counter+=1

        elif bed_tag==True:
             tm_tag = np.array(re.search('tm:Z:([-.0-9]*;)', line).group(1).replace(';','').split('.')).astype(int)

        else:
            tm_tag = [-1]

        seq_len = len(seq)
        mapq = int(mapq)
        flag = int(flag)
        start = int(start)-1
        end = start + seq_len
        tlen = abs(int(tlen))

        # If there is no tasmanian tag, confidence intersections table will include ALL intersections and complement nothing.
        # Here tag reads based on their level of confidence. If there was an intersection, use this level of confidence
        # include the bases in the 4th table.
        if bed_tag:
            confidence_value = int(re.search('tc:i:([0-9]*)', line).group(1))
        else:
            confidence_value = 0

        # reasons to skip reads
        # =====================
        # problems with CIGAR (* is data unavailable)
        # unwanted read length
        # Bad mapping quality
        # unrecognized chromosome

        if cigar=="*" or mapq==255 or chrom not in reference: 
            constants.SKIP_READS['others']+=1
            continue  
        elif constants.SKIP_INDEL and ("I" in cigar or "D" in cigar):
            constants.SKIP_READS['indel']+=1
            continue
        elif seq_len <constants.MIN_LENGTH or seq_len > constants.MAX_LENGTH:
            constants.SKIP_READS['readlength']+=1
            continue
        elif mapq<MinMapQuality:
            constants.SKIP_READS['mapquality']+=1
            continue
        elif tlen<constants.TLEN[0] or tlen>constants.TLEN[1]:
            constants.SKIP_READS['fragLength']+=1
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
                        if constants.SOFTCLIP_BYPASS==0:
                            s1, s2 = seq[:number], reference[chrom][start:start+number]
                            if simple_deltas_is_this_garbage(s1,s2):
                                seq = 'N' * number + seq[number:]
                                constants.SKIP_BASES['softclips']+=number
                        elif constants.SOFTCLIP_BYPASS==1:
                            seq = 'N' * number + seq[number:]
                            constants.SKIP_BASES['softclips']+=number 
                    else: # end
                        if constants.SOFTCLIP_BYPASS==0:
                            s1, s2 = seq[position:position+number], reference[chrom][start+position:start+position+number]
                            if simple_deltas_is_this_garbage(s1,s2):
                                seq = seq[:position] + 'N' * number
                                constants.SKIP_BASES['softclips']+=number
                        elif constants.SOFTCLIP_BYPASS==1:
                            seq = seq[:position] + 'N' * number
                            constants.SKIP_BASES['softclips']+=number

                elif i=='I':
                    seq_idx[position:position+number]=0
                
                elif i=='D':
                    ref_idx[position:position+number]=0
                
                elif i=='H':
                    ref_idx[position:position+number]=0
                    seq_idx[position:position+number]=0

                position += number
                last = current+1

        ref = reference[chrom][start:end]
        l = int(np.sort([np.sum(seq_idx), np.sum(ref_idx)])[0])
        try:
            seq = [seq[i] for i in range(len(seq_idx)) if seq_idx[i]==1][:l]
            ref = [ref[i] for i in range(len(ref_idx)) if ref_idx[i]==1][:l]
        except Exception as e:
            if debug:
                logger.warning('error: {} occurred while fixing for different lengths of seq and ref'.format(str(e)))
            pass

        # if bin(int(flag))[2:][-5]=='1' or make it easy for now...
        if flag==99 or (constants.ONT and flag in [0, 2048]): # for ont, not considering "secondary alignments", only "supp"
            strand='fwd'; read=1
        elif flag==163: 
            strand='fwd'; read=2
        elif flag==83 or (constants.ONT and flag in [16, 2064]):
            strand='rev'; read=1
        elif flag==147: 
            strand='rev'; read=2
        else: continue

        # At this point only accepted reads are analyzed. incorporate the fist X read length and later get mode
        if check_lengths_counter<100 and not constants.ONT:
            check_lengths.append(seq_len)
            check_lengths_counter +=1
        elif constants.ONT:
            if seq_len > check_lengths: # If ONT, check_lengths is a number not a list
                check_lengths = seq_len # THIS IS SAFE: if read length is > TLEN, it will be discarded

        if constants._UNMASK_GENOME: 
            ref = ''.join(ref).upper() # re-think before doing this.

        if len(seq) != len(ref):
            if debug:
                logger.warning('seq ={} and ref={} have different lengths'.format(seq,ref))
            continue    

        #if strand=='rev':
        #    rev_seq = revcomp(seq)
        #    sys.stderr.write(seq, rev_seq)

        # If there are more than "N" bases in the complement area and (of course) the read intersects a bed region
        # include complement and interections in differents tables
        


        for pos,base in enumerate(seq):

            if pos > len(ref) or pos > len(phred): 
                if debug:
                    logger.warning('ERROR processing read {}, with ref={}, phred={} and seq={}'.format(_,ref,phred,seq))
                continue

            if base=='N' or ref[pos]=='N' or ord(phred[pos]) < constants.PHRED: # or (CIGAR[pos] not in ['M','X','=']):
                #SKIP_BASES['quality']+=1  # report this in stderr. 
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
                        if confidence_value >= confidence:
                             errors_intersectionB[read][read_pos][ref_pos][Base.upper()] += 1
                        else:
                            errors_intersection[read][read_pos][ref_pos][Base.upper()] += 1

                    else:
                        assert base in ['A','C','G','T'], "{} should be upper case".format(base)
                        if confidence_value >= confidence:
                             errors_complementB[read][read_pos][ref_pos][Base.upper()] += 1
                        else:
                             errors_complement[read][read_pos][ref_pos][Base] += 1
                    
                    # write somewhere read_id, flag(read-number + strand), 
                    # read_position, chrom, genomic-position, mismatch-type.
                    if debug:
                        if ref_pos != Base and Base in ['A','C','T','G']:
                            sys.stderr.write('{},{},{},{},{},{},{}\n'.format(read_id, flag, read_pos, chrom, pos+start, ref_pos, Base))
                    
                    if constants.PWM:
                        # We have to fix this. For now, avoid positions too close to the ends of the reads
                        #if pos <=flanking_n:
                        #    this_seq = ref[0:pos*2+1] # keep it symetrical (I am not sure though)
                        #elif pos+flanking_n > len(ref):
                        #    continue
                        #else:
                        this_seq = ref[pos-flanking_n:pos+flanking_n+1] # Assuming 0-based index
                        
                        if strand == 'rev':
                            this_seq = [revcomp(b) for b in this_seq][::-1]
                        
                        # This will be replaced with a smarter option.
                        if len(this_seq) < flanking_n*2+1:
                            continue

                        this_seq = ''.join(this_seq)
                        fill_PFM(this_seq, PFM[''.join([ref_pos,Base])])

                except Exception as e:
                    if debug:
                        logger.warning('error:{} in chr:{}, position:{}, read:{}, base:{}, seq:{}, start:{} and ref_pos:{}'.format(\
                                                                    str(e), chrom, pos, read, base, ''.join(seq), start, ref[pos]))
    PWM = {}
    if PWM:
        for k,v in constants.PFM.items():
            PWM[k] = pfm2ppm(v)

        for k,v in PWM.items():
            if k[0]==k[1]: # careful dont make CC, GG, TT, AA to matrices of ones!!
                continue
            else:
                base_dist = PWM[k[0]+k[0]]
                PWM[k] = ppm2pwm(v, base_dist)

    # fix tables on length
    constants.READ_LENGTH = mode(check_lengths).mode if not constants.ONT else check_lengths
    #READ_LENGTH = np.max(check_lengths)
    
    if debug:
        logger.info('MODE READ LENGTH: {}'.format(str(constants.READ_LENGTH)))


    errors_intersection = trim_table(errors_intersection, constants.READ_LENGTH)
    errors_complement = trim_table(errors_complement, constants.READ_LENGTH)
    errors_unrelated = trim_table(errors_unrelated, constants.READ_LENGTH)
 
    errors_intersectionB = trim_table(errors_intersectionB, constants.READ_LENGTH)
    errors_complementB = trim_table(errors_complementB, constants.READ_LENGTH)

    #######################
    # REPORTING SECTION ##
    #######################
    
    # Report performance and less relevant statistics
    t2=time.time()
    sys.stderr.write('tasmanian finished the analysis in {} seconds \n'.format(str(t2-t0)))
    sys.stderr.write('reads discarded\n')
    sys.stderr.write('===============\n')
    for k,v in constants.SKIP_READS.items():
        sys.stderr.write('{:<15} {}\n'.format(k,v))
    sys.stderr.write('\nBases discarded\n')
    sys.stderr.write('===============\n')
    for k,v in constants.SKIP_BASES.items():
        sys.stderr.write('{:<15} {}\n'.format(k,v))

    # print table header
    #sys.stdout.write('read,position,' + ','.join \
    #    ([','.join([
    #        Sample+mut for mut in 'a_a,a_t,a_c,a_g,t_a,t_t,t_c,t_g,c_a,c_t,c_c,c_g,g_a,g_t,g_c,g_g'.split(',')
    #        ]) for Sample in ['I','C','N','cI','cC']
    #    ]) + '\n')

    # print rows
    prnt = []
    for read in [1,2]:
        for pos in np.arange(constants.READ_LENGTH): 
            prnt.append(str(read) + ',' + str(pos+1) + ',' + ','.join([ \
                ','.join([str(Table[read][pos][i][j]) for i,j in zip(['A','A','A','A','T','T','T','T','C','C','C','C','G','G','G','G'], ['A','T','C','G']*4)])
            for Table in [errors_intersection, errors_complement, errors_unrelated, errors_intersectionB, errors_complementB]]))

    #sys.stdout.write('\n'.join(prnt))
    
    # also save table as numpy arrays in dicts for ploting into the html report
    # we already printed the table as an array of strings. Just parse it into pandas and 
    # re-assign data-types
    col_names = ['read','position'] + ['_'.join([i[0].lower(), i[1].lower()]) for i in  
                    zip(['A','A','A','A','T','T','T','T','C','C','C','C','G','G','G','G'], ['A','T','C','G']*4)]

    df = pd.DataFrame([ i.split(',') for i in prnt]) 

    ids_intersection    = np.hstack([np.array([0,1]), np.arange(2,2+16)])
    ids_complement      = np.hstack([np.array([0,1]), np.arange(18,18+16)])
    ids_nonintersection =  np.hstack([np.array([0,1]), np.arange(34,34+16)]) 
    ids_intersectionB    = np.hstack([np.array([0,1]), np.arange(50,50+16)])
    ids_complementB      = np.hstack([np.array([0,1]), np.arange(66,66+16)])

    dfi = df.iloc[:, ids_intersection]
    dfc = df.iloc[:, ids_complement]
    dfn = df.iloc[:, ids_nonintersection] 
    dfiB = df.iloc[:, ids_intersectionB]
    dfcB = df.iloc[:, ids_complementB]

    table = {}
    for DF, dfName in zip([dfi, dfc, dfn, dfiB, dfcB],['intersection','complement','non_intersection', 'intersection_C','complement_C']):
        DF.columns = col_names
        table[dfName] = DF.astype(int)

    # Report errors to a logging file
    #f = open('errors_'+randLogName+'.log','w',)
    #f.write('\n'.join(logging))
    #f.close()
    
    return table, constants.PWM



if __name__=='__main__':

    # logging in debug mode
    if '--debugging-mode' in sys.argv:
        Args = sys.argv + ['--debugging-mode']
    else:
        Args = sys.argv

    # run tasmanian
    table, constants.PWM = analyze_artifacts(sys.stdin, Args) 

    # print results
    table_print = pd.concat(
        [i if n==0 else i.iloc[:,2:] for n,i in enumerate(table.values())]
        , axis=1)
        
    table_print.columns = ['read','position'] + ','.join([
        ','.join([t+i for i in table['intersection'].columns[2:]]) for t in ['I','C','N','cI','cC']
        ]).split(',')

    table_print.to_csv(sys.stdout, index=False)

    # avoid overrighting report file before saving.
    report_filename = 'Tasmanian_artifact_report.html'
    for n,i in enumerate(Args):
        if i in ['-o','--output-prefix']:
            report_filename = Args[n+1] + '.html'

    file_num = 0
    while os.path.isfile(report_filename) and report_filename == 'Tasmanian_artifact_report.html':
        file_num += 1
        report_filename = 'Tasmanian_artifact_report-' + str(file_num) + '.html'

    # create html template with plot
    report_html = plot_html(table)

    with open(report_filename, 'w') as f:
        f.write(report_html)

    # save PWM into pickle
    pwm_filename = 'Tasmanian_pwm_file' + str(file_num)
    with open(pwm_filename, 'wb') as f:
        pickle.dump(constants.PWM, f, protocol=pickle.HIGHEST_PROTOCOL)

    sys.stderr.write('\n' + report_filename + " and " + pwm_filename + " related files created\n")
