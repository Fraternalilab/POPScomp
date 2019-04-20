/*=============================================================================
topol : molecular topology
Copyright (C) 2002-2018 Franca Fraternali
Copyright (C) 2008-2018 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin and Valerie Hindie
Read the COPYING file for license information.
=============================================================================*/

#ifndef TOPOL_H
#define TOPOL_H

#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include "getpdb.h"
#include "matrix.h"
#include "sasa_const.h"
#include "type.h"

/*____________________________________________________________________________*/
/* structures */
/* topology */
typedef struct  
{
	int nMol; /* number of molecules */
	int nResidue; /* number of residues */
	int nBond; /* number of bonds */
	int nAngle; /* number of angles */
	int nTorsion; /* number of torsions */
	int nNonBonded; /* overlapping L-J (non-bonded) atoms */
	int *ib; /* atoms forming the bonds: array of size bond-number */
	int *jb;
	int *it; /* atoms forming the angles: array of size 3*bond-number */
	int *jt;
	int *kt;
	int *ip; /* atoms forming the torsions: array of size 5*bond-number */
	int *jp;
	int *kp;
	int *lp;
	int *in; /* atoms forming the L-J overlaps: array of size ~30*bond-number */
	int *jn; 
	Type (*angle)[3]; /* atom and residue types for each angle */
	int **bondState; /* records bonded atom pairs */
	int **neighbourState; /* records the neighboured (non-bonded) atom pairs */
	float **neighbourPar; /* records the POPS parameters 'p_ij * b_ij' of neighbours */
	int *interfaceNn; /* nearest neighbour on separate chain */
	float *interfaceNnDist; /* distance to nearest neighbour on separate chain */
} Topol;

/*____________________________________________________________________________*/
/* prototypes */
float atom_distance(Str *pdb, int i, int j);
float cutoff_radius(Type *type, ConstantSasa *constant_sasa, int i, int j, float rSolvent);
int get_bonds(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa, Argpdb *argpdb); /* calculate bonds (from atoms) */
int get_angles(Str *pdb, Topol *topol); /* calculate angles (from bonds) */
int get_torsions(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa); /* calculate torsions (from angles) */
int nonbonded_overlaps(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa, Arg *arg); /* calculate overlapping atoms */
void init_topology(Str *pdb, Topol *topol);
void free_topology(Str *pdb, Topol *topol);
int get_topology(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa, Argpdb *argpdb, Arg *arg); /* call topology routines */

#endif

