import numpy as np

class reads:
    '''
        read attrubutes are it's sam columns
    '''
    def __init__(self, _id, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred, tags):
        self._id = _id
        self.flag = int(flag)
        self.chrom = chrom
        self.start = int(start) #-1 --> correlated with line #139
        self.mapq = int(mapq)
        self.cigar = cigar
        self._2 = _2
        self._3 = _3
        self.tlen = int(tlen)
        self.seq = seq
        self.phred = phred
        self.tags = tags

        self.seq_len = len(seq)
        self.end = self.start + self.seq_len

        # How many bases in the complement side of an intersection?
        # if there is no intersection, complement will be length of the read.
        self.complement = self.seq_len 

        self.extended_cigar = None

        # according to the category of the read, get 2 informative positions: a & b
        # a. where the read intersects the bed fragment at the left (5')
        # b. where the read intersects the bed fragment at the right(3')
        self.category_positions = [None, None]

        # categories to classify the read as follows
        self.category = None
        ''' categories could be:
            1 - Unrelated
            2 - intersects (only read 1)
                a. completely
                b. partially left
                c. partially right
                d. both sides
            3 - intersects (only read 2)
                a. completely
                b. partially left
                c. partially right
                d. both sides
            4 - intersects BOTH READS
            5 - completely contained (mapq should be low?)
            6 - both reads contained in different repeats
            7 - both reads intersect different repeats
        '''
        # In masked_seq, all intersecting bases are replaced with "N"
        # intersect_seq is the complement (oposite) of masked
        self.masked_seq = None # update when reading bed and bam or leave it in tasmanian.
        self.intersect_seq = None # 
        self.junction = '-1;'  # updated in intersections

        # to know if the read and it's paired read intersect the same bed region
        # or not, we need to keep the bed id in the read
        self.bed_id = None

        # extra information to incorporate in the sam output file to correlate later 
        self.bed_extra_info = None


    def expand_cigar(self):
        ''' converts cigar into string to slice as needed'''
        cigarString = []
        num=''
        for i in self.cigar:
            if i.isnumeric(): 
                num = num + i
            else:
                cigarString.append(''.join([i]*int(num)))
                num=''

        self.expanded_cigar = ''.join(cigarString)
        return


    def print(self, sequence='original'):
        ''' prints sam read '''

        # This will allow to use 'masked' as argument as 'original'. One sentence instead of 
        # repeating this for read1 and read2 in the main.
        if self.masked_seq == None: 
            self.masked_seq = self.seq

        if sequence == 'original':      
            seq = self.seq
            tag = '0 '
        elif sequence == 'masked':      
            seq = self.masked_seq
            tag = '1 '
        elif sequence == 'intersect':   
            seq = self.intersect_seq
            tag = '2 '

        if self.category_positions == [None, None]: self.category_positions = ['nan','nan']
        categories = "categories " + ':'.join([str(i) for i in self.category_positions] + [str(self.category)])
        self.tags = self.tags + '\ttm:Z:' + self.junction  + '\ttc:i:' + str(self.complement)#tag + str(self.category) +' '+ self.junction

        try: 
            return '\t'.join([self._id, str(self.flag), self.chrom, str(self.start), 
                              str(self.mapq), self.cigar, self._2, self._3, 
                              str(self.tlen), seq, self.phred, self.tags
                             ])
        except Exception as e:
            tmp = '\t'.join([self._id, str(self.flag), self.chrom, str(self.start), 
                              str(self.mapq), self.cigar, self._2, self._3, 
                              str(self.tlen), seq, self.phred, self.tags
                            ])

            sys.stderr.write(self.bed_extra_info, 'problem: ', str(e), tmp)
