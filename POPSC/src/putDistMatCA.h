/*===============================================================================
putdistMatCA.h : routines for printing Calpha distance matrix
Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef PUTCADISTMAT_H
#define PUTCADISTMAT_H

#include <stdio.h>

#include "arg.h"
#include "getpdb.h"
#include "topol.h"

/*____________________________________________________________________________*/
/* print Calpha distance matrix */
void print_distMatCA(Arg *arg, Topol *topol);

#endif
