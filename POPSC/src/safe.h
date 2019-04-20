/*==============================================================================
safe.h : safe standard routines
Copyright (C) 2004 John Romein
Read the COPYING file for license information.
==============================================================================*/

#ifndef SAFE_H
#define SAFE_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/*___________________________________________________________________________*/
/* file */
FILE *safe_open(const char *name, const char *mode);

/*___________________________________________________________________________*/
/* allocation */
void *check_non_null(void *ptr);
void *safe_malloc(size_t size);
void *safe_realloc(void *ptr, size_t size);

#endif

