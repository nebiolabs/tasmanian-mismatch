#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "uthash.h"
#include "read_bed.h"

void bed_to_struct(FILE *bedfile, bed_coords_t **coords, bed_data_t **data){

	if (!bedfile) {fprintf(stderr, "Error: Bedfile not found\n"); return;}

	char *line=NULL;
	size_t n_characters=0, line_size=0;
	char *non_numeric_part; // to be used in strtoul
	bed_coords_t *s, *bed_coords = NULL;
	bed_data_t *s2, *bed_data = NULL;
	char **tokens;
	size_t n_tokens=0; 
	// reusables
	char *chrom, *rep_name, strand, *rep_class, chrom_start_end[30];
	unsigned long int start, end;

	while ((n_characters = getline(&line, &line_size, bedfile)) != -1) {
		if (line[n_characters - 1] == '\n') {
			line[n_characters - 1] = 0;
		}
		// some guess for empty lines. skip them	
		if (strlen(line) < 10) continue;

		// split the line into tokens separated by "\t"
		tokens = split(line, &n_tokens, "\t");

		// title should be skiped
		if (!isNumeric(tokens[1]) || !isNumeric(tokens[2])) continue;

		// collect tokens
		chrom     = tokens[0];
		start     = strtoul(tokens[1], &non_numeric_part, 10);
		end       = strtoul(tokens[2], &non_numeric_part, 10);
		rep_name  = tokens[5];
		strand    = *tokens[4];	 		
		rep_class = tokens[6];
		sprintf(chrom_start_end, "%s%s%ld%s%ld", chrom, "_", start,"_", end); 

		// if key exists in bed_coords, allocate memory for the arrays, start, stop in bed_data_.
		s = (bed_coords_t *) malloc(sizeof *s);
		HASH_FIND_STR(bed_coords, chrom, s);

		// if key is not yet in the hash table
		if (s == NULL){
			s = (bed_coords_t *) malloc(sizeof *s);

			s->chrom   = (char *) malloc(strlen(chrom) * sizeof(char*)); strcpy(s->chrom, chrom);
			s->size    = 100000;
        	s->start   = (unsigned long int *) malloc(s->size * sizeof s->start[0]);
        	s->end     = (unsigned long int *) malloc(s->size * sizeof s->end[0]);
			s->counter = 0;
			HASH_ADD_STR(bed_coords, chrom, s);
		}
		else {  // if we need to allocate more memory for start and end arrays
			if (s->size == s->counter){
				s->size  = s->size + 10000;
				s->start = realloc(s->start, s->size * sizeof (s->start[0]));
				s->end   = realloc(s->end, s->size * sizeof (s->end[0]));
			}
			// Only add to start and end arrays			
			s->counter++;
			s->start[ s->counter ] = start;
			s->end[ s->counter ] = end;
		}

		// FILL BED_DATA SECTION
		s2 = (bed_data_t *) malloc(sizeof *s2);
		s2->chrom_start_end = (char *) malloc( (strlen(chrom_start_end) + 1) * sizeof(char*));
		s2->rep_class_family = (char *) malloc( (strlen(rep_class) + 1) * sizeof(char*));
		s2->rep_name = (char *) malloc( (strlen(rep_name) +1) * sizeof(char*));

		strcpy(s2->chrom_start_end, chrom_start_end);
		strcpy(s2->rep_class_family, rep_class);
		strcpy(s2->rep_name, rep_name);
		s2->strand = strand;
		HASH_ADD_STR(bed_data, chrom_start_end, s2);

		for (int i=0;i<n_tokens; i++) memset(tokens[i],0,strlen(tokens[i]));
	}
	// FREE MEMORY SECTION
	free(tokens);
	free(s);
	free(s2);

	*coords = bed_coords;
	*data   = bed_data;
	//printf("\nadd***%p***ADD\n", *coords);
}

/*
int main(int argc, char **argv) {
	FILE *bedfile;
	bed_coords_t *coords;
	bed_data_t   *data;

	if ((bedfile = fopen(argv[1], "r")) != NULL) {
		bed_to_struct(bedfile, &coords, &data);
	}
	else {
		fprintf(stderr, "file %s not found. Aborting!\n", argv[1]);
		return 0;
	}
	//////
	printf("\nadd***%p***ADD..%zu\n",coords, coords->counter);
	bed_coords_t *ss;
	HASH_FIND_STR(coords, "chr1", ss);
	printf("%ld - %zd -- ",ss->start[10], ss->size);
	//////
	free(coords->start);
	free(coords->end);
	
	fclose(bedfile);
	return 0;
}
*/