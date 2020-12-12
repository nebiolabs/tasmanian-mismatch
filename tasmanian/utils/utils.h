#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

// Function to release arrays of pointers or strings from heap
// INPUT: Name of pointer to a array of strings
// OUTPUT: Nothing. It deallocates the array from memory
void free_buffer(char **buffer);

// Function to split an input line.
// INPUT: string (line), char (delimiter, optional) 
// OUTPUT: array of strings. Also, by reference, gives the 
//         number of strings in the array (n_items).
char **split_line(char *line, char *delim, size_t *n_items);

// Function to read the bedfile.
// It returns a struct with all components of each fragment.
// Fragments shouldbe sorted by chromosome and start coordinate
// in a linked list of bed_fragments.
// (added a type alias bed_fragment in the global namespace)
typedef struct bed_fragment{
    char chrom[6], rep_name[15], rep_class[15], rep_family[15], *strand;
    long int start, end;
} bed_fragment_t; 

// Function generated a bed_fragment_t from input file.
// INPUT: string (line), char (delimiter, optional)
// OUTPUT: bed_fragment_t struct.
bed_fragment_t new_fragment(char *line, char *delim);

#endif