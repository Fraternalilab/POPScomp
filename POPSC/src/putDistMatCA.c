/*==============================================================================
putDistMatCA.c : routines for printing Calpha distance matrix 
Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "arg.h"
#include "putDistMatCA.h"
#include "topol.h"

/*___________________________________________________________________________*/
/** print Calpha distance matrix */
void print_distMatCA(Arg *arg, Topol *topol)
{
	int i, j;

	arg->distMatCAOutFile = safe_open(arg->distMatCAOutFileName, "w");

	/* print residue numbers as column names */
	for (j = 0; j < topol->nCA2; ++ j) {
		if (j > 0) {
			fprintf(arg->distMatCAOutFile, " %d", topol->resCA2[j]);
		} else {
			fprintf(arg->distMatCAOutFile, "%d", topol->resCA2[j]);
		}
	}
    fprintf(arg->distMatCAOutFile, "\n");
		
	 for (i = 0; i < topol->nCA1; ++ i) {
        if (i > 0) fprintf(arg->distMatCAOutFile, "\n");

		/* print residue numbers as row names */
		fprintf(arg->distMatCAOutFile, "%d ", topol->resCA1[i]);
    	
		for (j = 0; j < topol->nCA2; ++ j) {
			if (j < topol->nCA2 - 1) {
       			fprintf(arg->distMatCAOutFile, "%3.2f ", topol->distMatCA[i][j]);
			} else {
       			fprintf(arg->distMatCAOutFile, "%3.2f", topol->distMatCA[i][j]);
			}
       	}
    }
    fprintf(arg->distMatCAOutFile, "\n");
	
	fclose(arg->distMatCAOutFile);
}

