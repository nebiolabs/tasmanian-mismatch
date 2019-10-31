import numpy as np
from sam_reads import reads

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



def assign_category(read1, read2):
    '''
        INPUT: read1 and read2 are reads (sam_reads)
            7 possible categories in sam_reads. The value is initialized as None and 
            both reads 1&2 share the same category 
        OUTPUT: category of both reads
    '''
    # case 1 was unrelated

    if read2.category == 1: # it's pairing reads case 2 (both unrelated were discared above)
        try:
            ab1 = read1.category_positions > 0
        except Exception as e:
            print('category positions is None and should be a number in read with id={} and error={}'.format(read1._id, str(e)))
            return None
        
        if ab1[0] and ab1[1]:           category = '2d'
        elif not ab1[0] and ab1[1]:     category = '2b'
        elif ab1[0] and not ab1[1]:     category = '2c'
        elif not (ab1[0] and ab1[1]):   category = '2a'

    elif read1.category==1: # it's case 3
        try:
            ab2 = read2.category_positions > 0
        except Exception as e:
            print('category positions is None and should be a number in read with id={} and error={}'.format(read2._id, str(e)))
            return None

        if ab2[0] and ab2[1]:           category = '3d'
        elif not ab2[0] and ab2[1]:     category = '3b'
        elif ab2[0] and not ab2[1]:     category = '3c'
        elif not (ab2[0] and ab2[1]):   category = '3a'

    elif read1.category != 1 and read2.category !=1: # it's cases 4 to 7
        try:
            ab1, ab2 = read1.category_positions > 0, read2.category_positions > 0
        except Exception as e:
            print('category positions (1 or 2) is None and should be a number in read with id={} or {} and error={}'.format(read1._id, read2._id, str(e)))
            return None

        if ab1[0] and ab2[1] and not (ab1[1] and ab2[0]): 
            if read1.bed_id == read2.bed_id:    category = '4'
            else :                              category = '7'

        elif not (ab1[0] and ab2[1] and ab1[1] and ab2[0]): 
           if read1.bed_id == read2.bed_id:     category = '5'
           else :                               category = '6'

    return category