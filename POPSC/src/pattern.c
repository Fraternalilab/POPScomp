/*==============================================================================
pattern.c : pattern matching routines
Copyright (C) 2009 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "pattern.h"

/*____________________________________________________________________________*/
/** compile single pattern */
void compile_pattern(regex_t *regex, char *pattern)
{
	/* regcomp outside assert(): it must also run when compiled with NDEBUG */
    if (regcomp(regex, pattern, REG_EXTENDED) != 0) {
		fprintf(stderr, "Error: Failed to compile pattern '%s'\n", pattern);
		exit(1);
	}
}

/*____________________________________________________________________________*/
/** compile multiple patterns */
void compile_patterns(regex_t *regex, char (*pattern)[32], int nPattern)
{
	int i;

	for (i = 0; i < nPattern; ++ i)
		compile_pattern(&(regex[i]), &(pattern[i][0]));
}

/*____________________________________________________________________________*/
/** match single pattern */
int match_pattern(regex_t *regex, char *searchString)
{
    return regexec(regex, searchString, 0, NULL, 0); 
}

/*____________________________________________________________________________*/
/** match multiple patterns */
int match_patterns(regex_t *regex, int nPattern, char *searchString)
{
	int i;

	for (i = 0; i < nPattern; ++ i) 
		if (regexec(&(regex[i]), searchString, 0, NULL, 0) == 0)
			return i; 

	return -1;
}

/*____________________________________________________________________________*/
/** free single pattern */
void free_pattern(regex_t *regex)
{
	regfree(regex);
}

/*____________________________________________________________________________*/
/** free multiple patterns */
void free_patterns(regex_t *regex, int nPattern)
{
	int i;

	for (i = 0; i < nPattern; ++ i)
		regfree(&(regex[i])); 
}

