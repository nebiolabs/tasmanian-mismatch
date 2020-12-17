#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include "uthash.h"

// Function to release arrays of pointers or strings from heap
// INPUT: Name of pointer to a array of strings
// OUTPUT: Nothing. It deallocates the array from memory
void free_buffer(char **buffer);

// Function to split an input line.
// INPUT: string (line), char (delimiter, optional) 
// OUTPUT: array of strings. Also, by reference, gives the 
//         number of strings in the array (n_items).
char **split_line(char *line, char *delim, size_t *n_items);

// Structure of a bedfile.
// It's a struct with all components of each fragment.
// Fragments should be sorted by chromosome and start coordinate
// in a linked list of bed_fragments.
// (added a type alias bed_fragment in the global namespace)
typedef struct bed_fragment{
    char chrom[6], rep_name[15], rep_class[15], rep_family[15], *strand;
    long int start, end;
} bed_fragment_t; 

// Function generates a bed_fragment_t from input file.
// INPUT: string (line), char (delimiter, optional)
// OUTPUT: bed_fragment_t struct.
bed_fragment_t new_fragment(char *line, char *delim);

// Structure, where the KEY is the CHROMOSOME
// each CHROMOSOME is a 2d-sorted-array with START and END values
// each chrX-start-end is the key for a second struct BED_FEATURES
// with the data for rep-family -class -name and strand
typedef struct chromosome_coords_bed{
    char *chrom; 
    long int start, end;
}chromosome_coords_bed_t;

typedef struct bed_features{
    char *key_chrom_start_end;
    char *rep_family, *rep_class, *rep_name, strand;
    UT_hash_handle hh; /* makes this structure hashable */
}bed_features_t;

#endif