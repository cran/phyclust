/* This file contains functions to call seq_gen_main(). in "seq-gen.c". */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

FILE *R_temp_file_pointer;
int seq_gen_main(int argc, char **argv);

SEXP R_seq_gen_main(SEXP R_argv, SEXP R_temp_file_name){
	int argc, i;
	const char **argv;
	const char *temp_file_name;

	argc = length(R_argv);
	argv = (const char**) malloc(argc * sizeof(char *));
	if(argv == NULL){
		error("Memory allocation fails!\n");
	}

	for(i = 0; i < argc; i++){
		argv[i] = CHAR(STRING_ELT(R_argv, i)); 
	}

	temp_file_name = CHAR(STRING_ELT(R_temp_file_name, 0));
	R_temp_file_pointer = fopen(temp_file_name, "w");

	GetRNGstate();		/* Get the seed from R. */
	seq_gen_main(argc, (char**) argv);
	PutRNGstate();		/* Update the seed of R. */

	fclose(R_temp_file_pointer);
	free(argv);
	return(R_NilValue);
} /* End of R_seq_gen_main(). */
