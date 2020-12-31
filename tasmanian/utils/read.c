#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void split(char *in, size_t *n, char *delimiter, char ***ptokens){
	if (delimiter == NULL) delimiter="\t"; // default
	
	char tokens[20][20]; // = malloc(20 * sizeof(char*));

	char *token = strtok(in, delimiter);
	while (token != NULL) {
		//tokens[*n] = malloc(30 * sizeof(char*)); // <- why not *tokens[0]
		strcpy(tokens[*n], token);
    	(*n)++;
		token = strtok(NULL, "\t");
	}
	
	*ptokens = tokens;
}

int main(int argc, char** argv){
	FILE *fp;
	
	fp = fopen(argv[1], "r");

	printf("abriendo %s",argv[0]);	

	char **tokens;
	*tokens = malloc(20 * sizeof *tokens);
	char *line=NULL;
	size_t n_characters, line_size=0;
	size_t n=0;
	while ((n_characters = getline(&line, &line_size, fp)) != -1) {
		if (line[n_characters - 1] == '\n') {
			line[n_characters - 1] = 0;
		}

		split(line, &n, "\t", &tokens);

		n++;
		printf("%ld -> %s\n", n, line);
	}
	fclose(fp);
	
	return 0;
}
