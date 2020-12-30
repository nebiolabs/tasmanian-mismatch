#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "uthash.h"


/*  Bed file is expected to have the following 7 fields:
 * chr1	19972	20727	(248935695)	+	L3	LINE/CR
 * corresponding to char, start, end, ?, rep_name, strand, rep_class/family
 * Bed file will be stored in two structures:
 * bed_coords: chromosomes are keys and start and end are values
 * bed_data: containes all other data with keys=chrom_start_end.
 */
typedef struct bed_coords_s{
	char chrom[25]; //key
	size_t size, counter; // n_size or count start end
	long unsigned int *start, *end;
	UT_hash_handle hh; /* makes this structure hashable key=chrom*/	
} bed_coords_t;

typedef struct bed_data_s{
	char chrom_start_end[25]; // key 
	char rep_name[15], rep_class_family[15], strand[3];
	UT_hash_handle hh; 
} bed_data_t;


char **split (char *in, size_t *n, char *delimiter){
	if (delimiter == NULL) delimiter="\t"; // default
	
	char **tokens = malloc(20 * sizeof(char*));

	char *token = strtok(in, delimiter);
	while (token != NULL) {
		tokens[*n] = malloc(30 * sizeof(char*)); // <- why not *tokens[0]
		strcpy(tokens[*n], token);
    	(*n)++;
		token = strtok(NULL, "\t");
	}
	return tokens;
}

bool isNumeric(const char *str){
	while (*str != '\0') {
		if (*str<'0' || *str>'9'){
			return false;
		}
		str++;
	}
	return true;
}


int main(int argc, char** argv){

	char *line=NULL;
	size_t n_characters=0, line_size=0;
	char *non_numeric_part; // to be used in strtoul
	bed_coords_t *s, *bed_coords = NULL;
	bed_data_t *s2, *bed_data = NULL;

	while ((n_characters = getline(&line, &line_size, stdin)) != -1) {
		if (line[n_characters - 1] == '\n') {
			line[n_characters - 1] = 0;
		}
		
		// some guess for empty lines. skip them	
		if (strlen(line) < 10) continue;

		// split the line into tokens separated by "\t"
		char **tokens;
		size_t n=0;
		tokens = split(line, &n, "\t");
	
		// title should be skiped
		if (!isNumeric(tokens[1]) || !isNumeric(tokens[2])) continue;
	
		char chrom[10]; strcpy(chrom,tokens[0]);
		unsigned long int start = strtoul(tokens[1], &non_numeric_part, 10);
		unsigned long int end   = strtoul(tokens[2], &non_numeric_part, 10);
		char *rep_name = tokens[4];
		char *strand   = tokens[5];	 		
		char *rep_class= tokens[6];
		char chrom_start_end[25]; 

		//printf("%s\t%ld\t%ld\t%s\t%s\t%s\n",chrom, start, end, rep_name, strand, rep_class);

		// if key exists in bed_coords, allocate memory for the arrays, start, stop in bed_data_.
		s = (bed_coords_t *) malloc(sizeof *s);
		HASH_FIND_STR(bed_coords, chrom, s);

		// if key is not yet in the hash table
		if (s == NULL){
			//free(s);
			s = (bed_coords_t *) malloc(sizeof *s);
			strcpy(s->chrom, chrom);
			s->size = 10000;
        	s->start = (unsigned long int *) malloc(s->size * sizeof s->start[0]);
        	s->end   = (unsigned long int *) malloc(s->size * sizeof s->end[0]);
			s->counter = 0;
			HASH_ADD_STR(bed_coords, chrom, s);
			printf("era NULL");	
		}
		// if we need to allocate more memory for start and end arrays
		else {
			if (s->size == s->counter){
				s->size = s->size + 1000;
				s->start = realloc(s->start, s->size * sizeof (s->start[0]));
				s->end = realloc(s->end, s->size * sizeof (s->end[0]));
			}
			// Only add to start and end arrays			
			s->counter++;
			s->start[ s->counter ] = start;
			s->end[ s->counter ] = end;
		}
/* */
		//free(s);
		s = (bed_coords_t *) malloc(sizeof *s);
		HASH_FIND_STR(bed_coords, chrom, s);
		printf("%ld - %zd -- ",s->start[s->counter-1], s->counter-1);

		//free(s);
		s = (bed_coords_t *) malloc(sizeof *s);
		HASH_FIND_STR(bed_coords, chrom, s);
		printf("%ld",s->start[5]);
/* */
		// FILL BED_DATA SECTION
		s2 = (bed_data_t *) malloc(sizeof *s2);
		sprintf(chrom_start_end, "%s%s%ld%s%ld", chrom, "_", start,"_", end);
		strcpy(s2->chrom_start_end, chrom_start_end);
		strcpy(s2->rep_class_family, rep_class);
		strcpy(s2->rep_name, rep_name);
		strcpy(s2->strand, strand);
		HASH_ADD_STR(bed_data, chrom_start_end, s2);

		//free(s2);
		s2 = (bed_data_t *) malloc(sizeof *s2);
		HASH_FIND_STR(bed_data, chrom_start_end, s2);
		printf("\n\n%s ********************************************\n",s2->rep_class_family);

		// FREE MEMORY SECTION
		for  (int i=0;i<n;i++){
			printf("\n%s, %d\n", tokens[i],i);
			free(tokens[i]);
		}
		free(tokens);
		//free(s);
		//free(s2);
		//printf("%s is %zu characters long", line, n_characters);

	}

	return(0);
}
