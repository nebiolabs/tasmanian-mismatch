#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

void free_buffer(char **buffer){
    register int i=0;
 
    while(buffer[i]) {
        free(buffer[i]);
        i++;
    }
    free(buffer);
}

char **split_line(char *line, char *delim, size_t *n_items){
    if (delim == NULL) delim = "\t"; // default is tab delimited
    char **results = malloc(15 * sizeof *results);
    *n_items=0;

    char *ptr = strtok(line, delim);
    while (ptr != NULL) {
        results[*n_items] = malloc(150 * sizeof *results[0]);  // double check why did I use results[0]...?
        strcpy(results[*n_items], ptr);
        ptr = strtok(NULL, delim);
        *n_items = *n_items + 1;
    }
    return results;
}

bed_fragment_t new_fragment(char *line, char *delim) {
    bed_fragment_t bf;
    if (delim == NULL) delim = "\t"; // default is tab separated values
    int n=0;

    // The order should be: chrom, start, end, strand, rep_name, rep_class, rep_family
    char *ptr = strtok(line, delim);
    while (ptr != NULL) {
        if (n==0) strcpy(bf.chrom, ptr);
        else if (n==1) bf.start = strtol(ptr, (char** )NULL, 10);
        else if (n==2) bf.end = strtol(ptr, (char** )NULL, 10);
        else if (n==3) strcpy(bf.strand, ptr);
        else if (n==4) strcpy(bf.rep_name, ptr);
        else if (n==5) strcpy(bf.rep_class, ptr);
        else if (n==6) strcpy(bf.rep_family, ptr);
        n++;
        ptr = strtok(NULL, delim);
    }
    return bf; 
}

// chromosomes has keys=chromosome name
// and values = 2d array [start, end]
chromosome_coords_bed_t* create_chromosome(int ROWS){
    chromosome_coords_bed_t *chromosome = malloc(sizeof(chromosome_coords_bed_t) + sizeof(int *));
    chromosomes->n_fragment = 0;
    chromosome->start_end = malloc(sizeof(int[ROWS][2]);
    // No need for initialization, as we keep the n_fragment
    // index or counter.
    // REMEMBER TO FREE THE HEAP!!!
    return chromosomes;
}

void increase_chromosome(chromosome_coords_bed_t **chromosome){
    *(chromosome)->start_end = (int *) realloc(*(chromosome)->start_end, sizeof(*(chromosome)->start_end) * 2);
}


chromosome_coords_bed_t chromosomes;
chromosomes = 

