#include <stdlib.h>
#include <stdio.h>
#include <string.h>


char **split (char *in, size_t *num_tokens, char *delimiter){
	if (delimiter == NULL) delimiter="\t"; // default

	size_t token_size=1;	
	*num_tokens=0;
	for (char *token=strtok(in,delimiter); token; token=strtok(NULL,delimiter)){
		(*num_tokens)++;
	}

	char **tokens = malloc(*num_tokens * sizeof(*tokens));
	char *token = in;
	for (size_t i=0; i < *num_tokens;i++){
		tokens[i] = token;
		token += strlen(token) + 1;
	}
	return tokens;
}


int main(int argc, char** argv){
	FILE *fp;
	
	fp = fopen(argv[1], "r");

	printf("abriendo %s",argv[0]);	

	char **tokens;
	char *line=NULL;
	size_t n_characters, line_size=0;
	size_t n=0;
	while ((n_characters = getline(&line, &line_size, fp)) != -1) {
		if (line[n_characters - 1] == '\n') {
			line[n_characters - 1] = 0;
		}

		tokens = split(line, &n, "\t");

		printf("%ld -> %s\n", n, line);
		for (int i=0; i<n;i++){
			printf("\t%s",tokens[i]);
		}
		printf("\n");
		//free(tokens);
		memset(&tokens,0,sizeof tokens);
	}
	fclose(fp);
	free(tokens);
	return 0;
}
