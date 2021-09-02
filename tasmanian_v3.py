#!/usr/bin/env python

from utils_v3 import *
import  multiprocessing as mp

GENOME = 'grch38_full.fa'
OVERLAPP_CONFIDENCE = 30
MAP_QUALITY = 20
SOFTCLIPS_OPTION = 'fix'

def _sam2table():

    def check1(i):
        # limit to paired, mapped in proper pair.
        # if CIGAR is unavailable or there are indels or suspicious chromosomes
        if i[1] in ['83','99','147','163'] and not \
            re.search(r'[ID]',i[5]) and i[2][:3]=="chr":
            return True
        else: 
            return False

    # read std input into list split by tabs
    tmp = [i.strip().split("\t") for i in open('test.sam','r')] # Changed this line stdin to filehandler.

    # Unfortunately, some flags are not present in all reads and hence,
    # the position of a columns may vary but if intersections was executed,
    # it's tag should be the -2 since tc:i:... is the -1.
    intersections_executed = False
    for n in range(10):
        if tmp[n][-2][:5] == 'tm:Z:':
            intersections_executed = True
    if intersections_executed:
        tmp  = [[i[n] for n in [1,2,3,4,5,9,10,11,-2]] for i in tmp if check1(i)] # if i[21][:5] == 'tm:Z:']
    else:
        tmp  = [[i[n] for n in [1,2,3,4,5,9,10,11]]+['tm:Z:-1'] for i in tmp if check1(i)]

    # Make it DataFrame for easier manipulation
    sam = pd.DataFrame(tmp)
    sam.columns = ['flag','chrom','position','qual','cigar','seq','bqual','md','tm']

    # save space in ram and make manupulations faster
    sam['flag'] = sam.loc[:,'flag'].astype(np.int16)
    sam['qual'] = sam.loc[:,'qual'].astype(np.int8)
    sam['position'] = sam.loc[:,'position'].astype(np.uint32)
    sam['chrom'] = sam['chrom'].astype('category')
    sam['md'] = sam['md'].astype('category')
    sam['tm'] = sam['tm'].astype('category') # lots of "tm:Z:-1"

    # I am not sure this helps but it shouldn't hurt
    del tmp

    # index for filtering by map-quality above the cutoff
    qual_idx = (sam['qual']>=MAP_QUALITY) & (sam['qual']!=255)

    # restrict reads to correctly mapped and in proper pair.
    # flag_id = (sam['flag']==83) | (sam['flag']==99) | (sam['flag']==147) | (sam['flag']==163)

    ## This is not strictly necessary if we can find the genome
    ## that was used to generate the bam file.
    # restrict to known chromosomes (we might change this later?)
    sam = sam.loc[sam['chrom'].str.contains("chr")]
    # EVEN HARSHER PURGING WOULD BE THIS
    sam = sam.loc[~sam['chrom'].str.contains("Un|random|decoy|alt", regex=True)]

    # filter using the above indices.
    # have to reset index to merge md and sam later on
    sam = sam.loc[qual_idx & flag_id].reset_index()

    # update sequence length to report reverse mismatches later
    sam.loc[:, 'seq_len'] = [len(i) for i in sam.loc[:, 'seq']]

    # mask bad qualities
    sam.loc[:, 'seq'] = [mask_bad_quality(q,s,l) for q,s,l in sam.loc[:, ['bqual','seq','seq_len']].values]

    return sam



def run_parallel(fxin):
    if fxin==1:
        return _sam2table()
    else:
        return load_genome(GENOME)



def main():

    inputs = [1,2]
    with mp.Pool(2) as p:
        [sam, genome] = p.map(run_parallel, inputs)


    # include sequence from the reference
    sam['ref_seq'] = [genome[i][j-1:j+k-1] for i,j,k in sam.loc[:,['chrom','position', 'seq_len']].values]

    # also restrict the reads and genome to ATCGN as with degenerate bases will be difficult to predict what the original base is.
    contain_weird_characters_in_seq = sam['ref_seq'].str.contains(r'[^ACTGN]')
    sam = sam.loc[~contain_weird_characters_in_seq,:]


    check_cpu_and_mem_usage()



    # r'^(\d+)S.*[A-Z](\d+)S$' will only extract when both groups are found ONLY.
    soft_idx = sam['cigar'].str.contains("S").values

    
    if SOFTCLIPS_OPTION=='eliminate':
        sam = sam.loc[~soft_idx,:]

    elif SOFTCLIPS_OPTION=='fix':
        S5 = sam.loc[soft_idx, 'cigar'].str.extract(r'^(\d+)S')
        S3 = sam.loc[soft_idx, 'cigar'].str.extract(r'(\d+)S$')

        temp = pd.concat([S5, S3], axis=1).fillna(0)
        temp.columns = ['S5', 'S3']

        # sys.getsizeof(new_sam) = sys.getsizeof(new_sam.fillna(0))
        sam = pd.concat([sam, temp], axis=1) #.fillna(0) --> can't because some categorical column already had NaN but not zeros.

        for k in ['S5','S3']:
            sam[k] = sam[k].fillna(0).astype(np.int8)

        # User should be able to choose on using or not softclips. If using them, we should decide if they are useful.
        # here assuming correction.

        sam.loc[soft_idx, 'seq'] = [
            correct_softclips(i,j,k,l) for i,j,k,l in sam.loc[soft_idx, ['seq','ref_seq','S5','seq_len']].values
        ]

        # also correct for S3 softclips
        sam.loc[soft_idx, 'seq'] = [
            correct_softclips(i,j,-k,l) for i,j,k,l in sam.loc[soft_idx, ['seq','ref_seq','S3','seq_len']].values
        ]


    check_cpu_and_mem_usage()


    #if OVERLAP_CONFIDENCE !=0:

    tms = sam['tm'].str.extractall(r'tm:Z:(\d+)\.(\d+);.*').astype(np.int16)
    tms.index = tms.index.droplevel(1)
    tms.columns = ['tm5','tm3']
    tms = pd.concat([sam.loc[tms.index, ['seq','seq_len']],tms], axis=1)


    overlapp_confidence = (tms['seq_len'] - (tms['tm3'] - tms['tm5']) >= OVERLAPP_CONFIDENCE).values
    overlapp_confidence_idx = tms.loc[overlapp_confidence].index

    tms['seq'] = tms.loc[overlapp_confidence, 'seq']

    sam.loc[overlapp_confidence_idx, 'seq'] = [
        k[0:i+1] + ''.join( (j-i) * ['N'] ) + k[j+1:] for i,j,k in tms.loc[overlapp_confidence_idx, ['tm5', 'tm3', 'seq']].values
    ]


    # all sequences aligned and padded so that align from the beginning
    read_length = np.unique([len(i) for i in sam['seq']], return_counts=True)
    read_length = read_length[0][ np.argmax(read_length[1]) ]

    sequences = sam.loc[:,'seq'].str.pad(read_length, side='right')
    sam_fwd = (sam.loc[:,'flag'].values == 163) | (sam.loc[:,'flag'].values == 99 )
    sam_rev = (sam.loc[:,'flag'].values == 83 ) | (sam.loc[:,'flag'].values == 147)
    sam_r1 =  (sam.loc[:,'flag'].values == 99 ) | (sam.loc[:,'flag'].values == 83 )
    sam_r2 =  (sam.loc[:,'flag'].values == 163) | (sam.loc[:,'flag'].values == 147)

    # also reverse complement
    sam.loc[sam_rev, 'ref_seq'] = [revcomp_seq(i) for i in sam.loc[sam_rev, 'ref_seq'].values]
    sam.loc[sam_rev, 'seq'] = [revcomp_seq(i) for i in sam.loc[sam_rev, 'seq'].values] # we also need to revcomp the read sequences

    # We also need to pad these sequences
    sam.loc[:, 'ref_seq'] = sam.loc[:, 'ref_seq'].str.pad(read_length, side='right')
    sam.loc[:, 'seq'] = sam.loc[:, 'seq'].str.pad(read_length, side='right')

    # revcomp removes ' ' so we have to padd them at the right now
    sam.loc[sam_rev, 'ref_seq'].str.pad(read_length, side='right')

    # update sequences to reflect reverse complements
    # sequences is a Series with same index as sam
    sequences.loc[sam_rev] = [revcomp_seq(i) for i in sequences.loc[sam_rev]] # IF INCLUDED REF_SEQ, WE NEED TO REVCOMP THAT TOO!!!!!
    sequences.loc[sam_rev] = sequences.loc[sam_rev].str.pad(read_length, side='right')

    # generate a dictionary with the total counts of each nucleotide at each position
    totals = { i:{'A':[], 'C':[], 'G':[], 'T':[]} for i in ['r1','r2']}

    table_inputs = [
        [read, values, base] for read,values,base in \
            product(['r1','r2'], [sam_r1,sam_r2], ['A','C','G','T'])
        ]

    def generate_tables(inputs_index):
        global table_inputs
        readNum, values, base = table_inputs[inputs_index]
        print(readNum, len(values), base)
        global sequences
        global read_length
        #for base in totals[readNum].keys(): # I could have chosen r2. Arbitrary decision.
        return [
            np.sum( np.array([seq[position] for seq in \
            sequences.loc[values].values]) == base ) for \
            position in range(read_length)
        ]

    with mp.Pool(len(table_inputs)) as p:
        results = p.map(generate_tables, [i for i in range(len(table_inputs))]) # [totals['r1'], totals['r2']] = 
    #results = map(generate_tables, inputs)

    for n,result in enumerate(results):
        read, val, base = inputs[n]
        totals[read][base] = result

    #for (read, values, base),result in zip(inputs):
    #    totals[read][base] = results[n]

    '''
    for base in totals['r1'].keys(): # I could have chosen r2. Arbitrary decision.
        totals['r1'][base] = [np.sum( np.array([seq[position] for seq in sequences.loc[sam_r1].values]) == base ) \
                            for position in range(read_length)]


        totals['r2'][base] = [np.sum( np.array([seq[position] for seq in sequences.loc[sam_r2].values]) == base ) \
                            for position in range(read_length)]
    '''
    '''
    # now excluding softclipped regions
    # ---------------------------------
    cigar_idx = ~sam['cigar'].str.contains("S").values # TRY including softclipped bases

    # generate a dictionary with the total counts of each nucleotide at each position
    totals_S = { i:{'A':[], 'C':[], 'G':[], 'T':[]} for i in ['r1','r2']}

    for base in totals_S['r1'].keys(): # I could have chosen r2. Arbitrary decision.
        totals_S['r1'][base] = [np.sum( np.array([seq[position] for seq in sequences.loc[sam_r1 & cigar_idx].values]) == base ) \
                            for position in range(read_length)]


        totals_S['r2'][base] = [np.sum( np.array([seq[position] for seq in sequences.loc[sam_r2 & cigar_idx].values]) == base ) \
                            for position in range(read_length)]
    '''

    sam.loc[:,'diff'] = [diff(i,j) for i,j in sam.loc[:,['ref_seq','seq']].values]

    f,ax = plt.subplots(1,2, figsize=(10,5));

    # including soft-clips
    positions = sam.loc[sam_r2,'diff'].str.extractall("C(\d+)T").values
    positions = np.hstack(positions).astype(int)
    x,y = np.unique(positions, return_counts=True)
    base_distribution = np.array(totals['r2']['C'])[x]

    ax[0].scatter(x,y/base_distribution)
    ax[0].grid()


    plt.show()


if __name__ == '__main__': main()
