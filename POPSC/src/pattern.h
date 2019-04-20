/*==============================================================================
pattern.h : pattern matching routines
Copyright (C) 2009 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef PATTERN_H
#define PATTERN_H

#include <assert.h>
#include <regex.h>
#include <safe.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*___________________________________________________________________________*/
/* compile single or multiple patterns */
 void compile_pattern(regex_t *regex, char *pattern);
 void compile_patterns(regex_t *regex, char (*pattern)[32], int nPattern);

/*___________________________________________________________________________*/
/* match single or multiple patterns */
int match_pattern(regex_t *regex, char *searchString);
int match_patterns(regex_t *regex, int nPattern, char *searchString);

/*___________________________________________________________________________*/
/* free single or multiple patterns */
void free_pattern(regex_t *regex);
void free_patterns(regex_t *regex, int nPattern);

#endif

