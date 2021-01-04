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

// look https://stackoverflow.com/questions/65526234/
// function-in-c-referencing-memory-externally-allocated
char **split (char *in, size_t *num_tokens, char *delimiter){
    if (delimiter == NULL) delimiter="\t"; // default

    size_t token_size=1;
    *num_tokens=0;
    for (char *token=strtok(in,delimiter); token; token=strtok(NULL,delimiter)){
        (*num_tokens)++;
    }
	// double loop to allocate the EXACT needed memory 
    char **tokens = malloc(*num_tokens * sizeof(*tokens));
    char *token = in;
    for (size_t i=0; i < *num_tokens;i++){
        tokens[i] = token;
        token += strlen(token) + 1;
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

