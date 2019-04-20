/*==============================================================================
sasa_const.h : SASA computation
Copyright (C) 2002-2018 Franca Fraternali
Copyright (C) 2008-2018 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin and Valerie Hindie
Read the COPYING file for license information.
==============================================================================*/

#if !defined SASACONST_H
#define SASACONST_H

#include <stdlib.h>

typedef struct
{
	char residueName[4];
	char atomName[4];
	float radius; /* atom radius */
	float parameter; /* atom-specific  SASA parameter */
	int polarity; /* hydrophobic = 0, hydrophilic = 1 */
	int ring; /* not in ring = 0, in ring = 1 */
	double surface; /* SASA of fully exposed atom */
} AtomdataSasa;

/*___________________________________________________________________________*/
/* structures */
/* constant SASA data */
typedef struct {
    int nResidueType; /* number of residue types */
    int nAtomResidue[64]; /* number of atoms per residue */
    AtomdataSasa atomDataSasa[64][64]; /* atom-specific data */
	float connect_12_parameter; /* connectivity parameter for 1-2 interactions */
	float connect_13_parameter; /* connectivity parameter for 1-3 interactions */
	float connect_14_parameter; /* connectivity parameter for 1-4 interactions */
	float connect_15_parameter; /* connectivity parameter for >(1-4) (Lennard-Jones) interactions */
	float rSolvent; /* solvent radius */
} ConstantSasa;

#endif
