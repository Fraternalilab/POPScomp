/*==============================================================================
putCAdistMat.c : routines for printing Calpha distance matrix 
Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "arg.h"
#include "putCAdistMat.h"
#include "topol.h"

/*___________________________________________________________________________*/
/** print Calpha distance matrix */
void print_CAdistMat(Arg *arg, Topol *topol)
{
	unsigned int i, j;

	safe_open(arg->CAdistMatOutFileName, "w");

	for (i = 0; i < topol->nCA; ++ i) { 
		for (j = 0; j < topol->nCA; ++ j) { 
			print_mat2D_float(arg->CAdistMatOutFileName,
								topol->distMatCA, i, j);
		}
	}

	fclose(arg->CAdistMatOutFile);
}

