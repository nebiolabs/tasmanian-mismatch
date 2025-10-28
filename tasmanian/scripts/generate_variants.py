import sys, os
from simulate_reads import load_reference_genome

reference = sys.argv[1]
genome = load_reference_genome(reference)


def read_proposed_variants(vars_file):
    '''
    vars_file is a tab separated file (TSV) and NO HEADER
    each line in vars_file contains original_base, variant_base, frequency(%), strand
    e.g. C  T   10  +
         G  A   0.8 - (yes. 0.8%)
    '''
    with open(vars_file, 'r') as f:
        pre_vars = [i.strip().split("\t") for i in f]

    
    variants = [
        {
        'original':  i[0],
        'variant':   i[1],
        'frequency': i[2],
        'strand':    i[3]
        }
        for i in pre_vars
    ]

    return variants


def initialize_genome_probabilities(genome_dict):
    ''' All positions have only one option initially'''
    return {contig: np.ones(len(sequence)) for contig, sequence in genome_dict.items()}


def generate_single_variants(genome_dict, genome_probs, single_variant, contig, locus):
    '''
    The frequencies are managed in this script,  through Copy Number Alterations (CNAs).
    Initially, we will assume that the amplified loci (arbitrarily set here) are placed
    at the end of the contig (not true in reality).
    Then, the genome_probabilities (to be sequenced) will be reduced by half every time 
    we duplicate these gneomic positions.
    INPUT:
        genome_dict,
        genome_probs (probabilities for each position to be sequenced),
        single_variant (dictionary with orig, var, freq and strand),
        contig, 
        locus (random start position). 
    OUTPUT:
        genome_dict  -> modified! (omitted)
        genome_probs -> modified! (omitted)
        {locus: [base, var, freq, length]}
    '''
    # randombly pick a "locus" (and assign a random length <100bp for now), find any "original" base and mutate it to "variant"
    # then replicate this region and assign probabilities of being selected 
    # appropriately.
    
    # define locus and randomly select the base to mutate
    np.random.seed(seed)
    locus_length = np.random.randint(50,100,1)
    sequence = genome_dict[contig][locus:locus+locus_length+1]
    idx = np.array( [i for i in re.finditer( single_variant['original'], sequence )] )
    np.random.shuffle(idx)
    mut_pos = idx[0]

    if single_variant['strand'] == "+":
        var = single_variant['variant']
    else:
        var = complement( single_variant['variant'] )
    
    # update probabilities and generate new genome including the mutated version of locus.
    genome_probs[contig][mut_pos:mut_pos+locus_length+1] /=2
    l = len(genome_probs[contig])
    genome_probs[contig][l:l+locus_length+1] = genome_probs[contig][mut_pos:mut_pos+locus_length+1]

    new_sequence = sequence[:mut_pos] + single_variant[variant] + sequence[mut_post+1:]
    genome_dict[contig] = genome_dict[contig] + new_sequence

    locus_dict = {
        '_'join([contig, 
                 str(locus), 
                 str(l), 
                 str(mut_pos),
                ':'.join(single_variant.values())
                ])
    }

    return {contig, str(locus), str



def generate_variants():
    global seed
    np.random.seed(seed)
    # insert length = WITHOUT adapter

    # Distribute the reads evenly, based on the length of each contig.
    contig_lengths = {name: len(seq) for name, seq in genome_dict.items()}
    total_length = np.sum([v for v in contig_lengths.values()])
    read_unit = n_reads / total_length

    contig_reads = {name: int(length * read_unit) for name, length in contig_lengths.items()}
