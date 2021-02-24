#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "read_bed.h"
#include "uthash.h"

/* We are reading a sam file. We store X number of reads in a array of strings
 * in the heap with calloc and we have to release it.
 * We will read MAX_LINES lines and will add more reads as needed with realloc.
 */
#define MAX_LINES 255 
#define LENGTH_LINE 500  // Chances that a line has more than 500 characters are very low.
char* HELP_MESSAGE = "\nScript should be called: ...\n";


int main(int argc, char **argv){

    /*
    if (argc<2){
        fprintf(stderr, "%s", HELP_MESSAGE);
        return 1;
    }
    */

    char *line = NULL; 
    ssize_t n_chars_line = 0;                   // characters read by getline
    size_t line_size = 0, n_lines = 0;          // set to zero = no length restrictions.
    char **buffer = NULL;       
    char **buffer_tmp = NULL;
    size_t limit_lines = MAX_LINES, new_limit_lines = 0;
    int finished_header = 0; 
    // variables to use in splitting the input strings into arrays of strings 
    char** splitted_line;
    size_t n_items;

    // each line will be processed as a string array and dynamically appended to matrix.
    buffer = calloc(MAX_LINES, sizeof(*buffer)); 

    while ((n_chars_line = getline(&line, &line_size, stdin)) != -1) {
        if (line[n_chars_line-1] == 0xa) {      // if feed line 0xa = '\n'
            line[n_chars_line-1] = 0;           // strip newline
            n_chars_line--;
        }
        if (line == NULL) continue;
        else if (line[0]=='@') { 
            finished_header++;
            printf("%s\n",line);
            continue;
        } else if (finished_header>0) printf("@PG\tID:tasmanian-mismatch\tVN:0.0.1\n");

        if (n_lines == limit_lines -1) {
            new_limit_lines = limit_lines + MAX_LINES;  // ToDo: decide which increase is better

            if ((buffer_tmp = calloc(new_limit_lines, sizeof(*buffer)))) {
                // copy matrix into first part of tmp_matrix
                if (memcpy(buffer_tmp, buffer, n_lines * sizeof(*buffer)) == buffer_tmp) {
                    free(buffer);
                    buffer = buffer_tmp;
                }
                else {
                    fprintf(stderr, "Error: memcpy failed in line 53\n");
                    return 1;
                }
            }
            else {
                fprintf(stderr, "Error: buffer_tmp allocation failed in line 51\n");
                return 1;
            }
            limit_lines = new_limit_lines;
        }
        buffer[n_lines] = strdup(line); // if no posix, use malloc and strcpy
        n_lines++;
    }
    if (line) free(line);

    n_lines=0;
    while(buffer[n_lines]){
        splitted_line = split(buffer[n_lines], &n_items, "\t");
        for (int i=0; i<n_items; i++) printf("%s - ", splitted_line[i]);
        printf("\n");

        printf("line[%zu] says %s\n", n_lines, buffer[n_lines]);
        free(buffer[n_lines]);
        n_lines++;
    }
    free(buffer);

    return 0;
}