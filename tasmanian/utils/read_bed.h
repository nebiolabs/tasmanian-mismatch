#ifndef READ_BED_H_INCLUDED
#define READ_BED_H_INCLUDED

#include "uthash.h"
#include "utils.h"

/*  Bed file is expected to have the following 7 fields:
 * chr1	19972	20727	(248935695)	+	L3	LINE/CR
 * corresponding to char, start, end, ?, rep_name, strand, rep_class/family
 * Bed file will be stored in two structures:
 * bed_coords: chromosomes are keys and start and end are values
 * bed_data: containes all other data with keys=chrom_start_end.
 */
typedef struct bed_coords_s{
	char *chrom; //key
	size_t size, counter; // n_size or count start end
	long unsigned int *start, *end;
	UT_hash_handle hh; /* makes this structure hashable key=chrom*/	
} bed_coords_t;

typedef struct bed_data_s{
	char *chrom_start_end; // key 
	char *rep_name, *rep_class_family, strand;
	UT_hash_handle hh; 
} bed_data_t;

// Function populates a bed_coords_t and a bed_data_t
// structs from input file.
// INPUT: BED filename (str), pointers to coords and data
// OUTPUT: bed_fragment_t struct.
void bed_to_struct(FILE *bedfile, bed_coords_t **coords, bed_data_t **data);

#endif