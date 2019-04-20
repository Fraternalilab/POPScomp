/*===============================================================================
putsigma.h : routines for printing sigma and SFE values
Copyright (C) 2011 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef PUTSIGMA_H
#define PUTSIGMA_H

#include <stdio.h>

#include "arg.h"
#include "argpdb.h"
#include "getpdb.h"
#include "sigma.h"
#include "sigma_const.h"
#include "type.h"

/*____________________________________________________________________________*/
/* print sigma and Solvation Free Energy values */
void print_sfe(Arg *arg, Argpdb *argpdb, Str *pdb, Type *type, Topol *topol, \
				MolSFE *molSFE, ConstantSigma *constant_Sigma, int frame);

#endif
