import sys, os
sys.path.append(os.path.realpath("./"))
from simulate_mutations import *

read_length = 76
genome_dict = load_reference_genome("../../test_data/test_genome.fa")
variants_list = read_proposed_variants('variants_proposed.tsv')
variants_dict = assign_contigs_to_variants(variants_list, genome_dict)

reads = generate_reads(
    genome_dict=genome_dict,
    variants_dict=variants_dict,
    read_length=read_length,
    n_reads=1000000,
    insert_length=100,
    insert_length_var=35
)

probs = generate_probabilities(read_length=read_length)

mutated_reads = mutate_reads(reads=reads, probabilities=probs, read_length=read_length)

assert reads != mutated_reads, "HEY!!! reads and mutated reads are the same!!!"


reads_dict = print_reads(mutated_reads)

with open('simulated_reads_R1.fastq', 'w') as f:
    f.write(reads_dict[0])

with open('simulated_reads_R2.fastq', 'w') as f:
    f.write(reads_dict[1])
