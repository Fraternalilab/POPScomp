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
	arg->distMatCAOutFile = safe_open(arg->distMatCAOutFileName, "w");
	print_mat2D_float(arg->distMatCAOutFileName, topol->distMatCA, topol->nCA1, topol->nCA2);
	fclose(arg->distMatCAOutFile);
}

