#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "uthash.h"

/* We are reading a sam file. We store X number of reads in a array of strings
 * in the heap with calloc and we have to release the space in ram heap.
 * We will read MAX_LINES lines and will add more reads as needed with realloc.
 */
#define MAX_LINES 255 
#define LENGTH_LINE 500  // Chances that a line has more than 500 characters are very low.
char* HELP_MESSAGE = "\nScript should be called: ...\n";


int main(int argc, char* argv){

    /*
    if (argc<2){
        fprintf(stderr, "%s", HELP_MESSAGE);
        return 1;
    }
    */

    char *line = NULL; 
    ssize_t n_chars_line = 0;                   // characters read by getline
    size_t n = 0;                               // set to zero = no length restrictions.
    char **matrix = NULL;       
    char **matrix_tmp = NULL;
    int n_lines = 0;                            // # lines in matrix so far
    size_t limit_lines = MAX_LINES;
    size_t new_limit_lines = 0;
    // variables to use in splitting the input strings into arrays of strings 
    char** splitted_line;
    size_t n_items;

    // each line will be processed as a string array and dynamically appended to matrix.
    matrix = calloc(MAX_LINES, sizeof(*matrix)); 

    while ((n_chars_line = getline(&line, &n, stdin)) != -1) {
        if (line[n_chars_line-1] == 0xa) {      // if feed line
            line[n_chars_line-1] = 0;           // strip newline
            n_chars_line--;                 
        }
        if (n_lines >= limit_lines -1) {
            new_limit_lines = limit_lines * 2;  // ToDo: decide which increase is better

            if (line == NULL) continue;  

            if (matrix_tmp = calloc(new_limit_lines, sizeof(*matrix))) {
                // copy matrix into first half of tmp_matrix
                if (memcpy(matrix_tmp, matrix, n_lines * sizeof(*matrix)) == matrix_tmp) {
                    free(matrix);
                    matrix = matrix_tmp;
                }
                else {
                    fprintf(stderr, "Error: memcpy failed in line 63\n");
                    return 1;
                }
            }
            else {
                fprintf(stderr, "Error: matrix_tmp allocation failed in line 61\n");
                return 1;
            }
            limit_lines = new_limit_lines;
        }
        matrix[n_lines] = strdup(line);
        n_lines++;
    }
    if (line) free(line);

    n_lines=0;
    while(matrix[n_lines]){
        splitted_line = split_line(matrix[n_lines], "\t", &n_items);
        for (int i=0; i<n_items; i++) printf("%s - ", splitted_line[i]);
        printf("\n");

        printf("line[%d] says %s\n", n_lines, matrix[n_lines]);
        free(matrix[n_lines]);
        n_lines++;
    }
    free(matrix);

    FILE *fp = fopen("./test.bed", "r");
    char *line_buf=NULL;
    size_t line_buf_size=0;
    ssize_t line_size;
    bed_fragment_t bed_line;

    if (!fp) {fprintf(stderr, "Error: file not found"); return 1;}

    line_size = 1;
    while (line_size>=0) {
        line_size = getline(&line_buf, &line_buf_size, fp);
        line_buf[line_size-1]=0;
        printf("bedfile %s - LALALALLA\n", line_buf);
        bed_line = new_fragment(line_buf, "\t");
        printf("%s - %ld - %ld - %s - %s - %s - %s\n", bed_line.chrom, bed_line.start, 
                                                     bed_line.end, bed_line.strand, 
                                                     bed_line.rep_name, bed_line.rep_class, 
                                                     bed_line.rep_family);
    }
    free(line_buf);



    // CREATE STRUCT
    struct my_struct {
        char *id;            /* we'll use this field as the key */
        char name[10];
        UT_hash_handle hh; /* makes this structure hashable */
    };

    struct my_struct *users = NULL;

    // ADD USER
    void add_user(char *id, char *name) {
        struct my_struct *s;
        s = malloc(sizeof(struct my_struct));
        strcpy(s->id, id);
        strcpy(s->name, name);
        HASH_ADD_STR( users, id, s );  /* id: name of key field */
    }

    add_user("testing-cualquiera-1234", "whatever");


    // FIND USER
    struct my_struct *find_user(char *user_id) {
        struct my_struct *s;

        HASH_FIND_STR( users, user_id, s );
        return s;
    }

    struct my_struct *s = find_user("testing-cualquiera-1234");

    printf("*********************************************** %s",s->id);
    


    return 0;
}