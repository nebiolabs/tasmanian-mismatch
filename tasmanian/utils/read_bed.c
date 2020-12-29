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
	//bed_coords = (bed_coords_t *) malloc(sizeof *bed_coords);
	
	bed_data_t *bed_data_tmp, *bed_data = NULL;


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

		printf("%s\t%ld\t%ld\t%s\t%s\t%s\n",chrom, start, end, rep_name, strand, rep_class);

		// if key exists in bed_coords, allocate memory for the arrays, start, stop in bed_data_.
		s = (bed_coords_t *) malloc(sizeof *s);
		HASH_FIND_STR(bed_coords, chrom, s);

		// if key is not yet in the hash table
		if (s == NULL){
			s = (bed_coords_t *) malloc(sizeof *s);
			strcpy(s->chrom, chrom);
			s->size = 10;
        	s->start = (unsigned long int *) malloc(s->size * sizeof s->start[0]);
        	s->end   = (unsigned long int *) malloc(s->size * sizeof s->end[0]);
			s->counter = 0;
			HASH_ADD_STR(bed_coords, chrom, s);
			printf("era NULL");	
		}
		// if we need to allocate more memory for start and end arrays
		else {
			if (s->size == s->counter){
				printf("size=%ld - counter=%ld",s->size, s->counter);
				s->size = s->size + 10;
				printf("size=%ld - counter=%ld",s->size, s->counter);
				s->start = realloc(s->start, s->size * sizeof (s->start[0]));
				s->end = realloc(s->end, s->size * sizeof (s->end[0]));
			}
			// Only add to start and end arrays			
			s->counter++;
			s->start[ s->counter ] = start;
			s->end[ s->counter ] = end;

			printf("yeag!!!\n");
		}
/*
		s = (bed_coords_t *) malloc(sizeof *s);
		HASH_FIND_STR(bed_coords, chrom, s);
		printf("%ld - %zd -- ",s->start[s->counter-1], s->counter-1);

		a = (bed_coords_t *) malloc(sizeof *s);
		HASH_FIND_STR(bed_coords, chrom, s);
		printf("%ld",s->start[5]);
*/

/*
		bed_coords_tmp = (bed_coords_t *)malloc(sizeof *bed_coords_tmp);
		bed_data_tmp   = (bed_data_t *)  malloc(sizeof *bed_data_tmp);		


		sprintf(bed_data_tmp->chrom_start_end, "%s%s%s%s%s", bed_coords_tmp->chrom, "_", start,"_", end);
		bed_coords_tmp->start = (unsigned long int *) strtoul(start, &non_numeric_part, 10);
		//bed_coords_tmp->end = (unsigned long int *) strtoul(end, &non_numeric_part, 10);

		printf("RRRRRRRRRRRRRRRRRRR %s RRRRRRRRRRRRRRRR\n", start); //bed_coords_tmp->start);

		//printf("%s \n", bed_data_tmp->chrom_start_end);

		//HASH_ADD_STR(bed_coords, bed_coords_tmp->chrom, bed_coords_tmp);
		//HASH_ADD_STR(bed_data, bed_data_tmp->chrom_start_end, bed_data_tmp);

			for  (int i=0;i<n;i++){
				printf("\n%s, %d\n", tokens[i],i);
				free(tokens[i]);
			}
			free(tokens);
			//printf("%s is %zu characters long", line, n_characters);
*/
		}

	return(0);
}
