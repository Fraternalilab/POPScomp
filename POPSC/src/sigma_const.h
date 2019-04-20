/*==============================================================================
sigma_const.h : sigma values for Solvation Free Energy (SFE) computation
Copyright (C) 2011 Franca Fraternali
Copyright (C) 2011 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef SIGMACONST_H
#define SIGMACONST_H

#include <stdlib.h>

typedef struct
{
	char residueName[4];
	char atomName[4];
	float sigma_type; /* sigma atomtype parameter */
	float sigma_group; /* sigma atomgroup parameter */
} AtomdataSigma;

/*___________________________________________________________________________*/
/* structures */
/* constant SIGMA data */
typedef struct
{
    int nResidueType; /* number of residue types */
    int nAtomResidue[64]; /* number of atoms per residue */
    AtomdataSigma atomDataSigma[64][64]; /* atom-specific data */
} ConstantSigma;

#endif
