#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <stdbool.h>
#include "uthash.h"

// Function to release arrays of pointers or strings from heap
// INPUT: Name of pointer to a array of strings
// OUTPUT: Nothing. It deallocates the array from memory
void free_buffer(char **buffer);

// Function to split an input line.
// INPUT: string (in), pointer to size_t, char (delimiter, optional) 
// OUTPUT: array of strings. Also, by reference, gives the 
//         number of strings in the array (num_tokens).
char **split (char *in, size_t *num_tokens, char *delimiter);

// Function checks if string is a number
// INPUT = string
// OUTPUT = Boolean (True or False)
bool isNumeric(const char *str);

#endif