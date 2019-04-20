/*==============================================================================
safe.c : safe standard routines
Copyright (C) 2004 John Romein
Read the COPYING file for license information.
==============================================================================*/

#include "safe.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*___________________________________________________________________________*/
/** safe file opening */
FILE *safe_open(const char *name, const char *mode)
{
    FILE *file = fopen(name, mode);
    if (file) {
		return file;
	} else {
		fprintf(stderr, "Error: Failed accessing file '%s'\n", name);
		exit(1);
	}
}

/*___________________________________________________________________________*/
/** safe memory allocation */
void *check_non_null(void *ptr)
{
    assert (ptr != 0);
    return ptr;
}

void *safe_malloc(size_t size)
{
    return check_non_null(malloc(size));
}

void *safe_realloc(void *ptr, size_t size)
{
    return check_non_null(realloc(ptr, size));
}

