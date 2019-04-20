/*==============================================================================
sasa.h : SASA computation
Copyright (C) 2002-2018 Franca Fraternali
Copyright (C) 2008-2018 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin and Valerie Hindie
Read the COPYING file for license information.
==============================================================================*/

#if !defined SASA_H
#define SASA_H

#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "topol.h"
#include "type.h"
#include "sasa_const.h"
#include "vector.h"

/*___________________________________________________________________________*/
/* structures */
/* Atomic Solvent Accessible Surface Area */
typedef struct  
{
    double sasa; /* SASA */
	double surface; /* surface area of isolated atom */
	int nOverlap; /* number of overlaps */
	int polar; /* polarity */
    double phobicbSasa; /* hydrophobic bSASA */
    double philicbSasa; /* hydrophilic bSASA */
    double bSasa; /* bSASA */
} AtomSasa;

/* Residuic Solvent Accessible Surface Area */
typedef struct  
{
	double surface; /* surface area of isolated atoms in residue */
	double phobicSasa; /* hydrophobic SASA */
	double philicSasa; /* hydrophilic SASA */
    double sasa; /* SASA */
	int nOverlap; /* number of overlaps */
	int atomRef; /* first atom of residue to refer to residue number and name */
    double phobicbSasa; /* hydrophobic bSASA */
    double philicbSasa; /* hydrophilic bSASA */
    double bSasa; /* bSASA */
} ResSasa;

/* Chainic Solvent Accessible Surface Area */
typedef struct  
{
	double surface; /* surface area of isolated atoms in chain */
	double phobicSasa; /* hydrophobic SASA */
	double philicSasa; /* hydrophilic SASA */
    double sasa; /* SASA */
	int first, last; /* first/last atom of chain */
    double phobicbSasa; /* hydrophobic bSASA */
    double philicbSasa; /* hydrophilic bSASA */
    double bSasa; /* bSASA */
} ChainSasa;

/* Molecular Solvent Accessible Surface Area */
typedef struct  
{
	AtomSasa *atomSasa; /* atomic values */
	ResSasa *resSasa; /* molecular values */
	ChainSasa *chainSasa; /* chain values */
    double philicSasa; /* hydrophilic SASA */
    double phobicSasa; /* hydrophobic SASA */
    double sasa; /* SASA */
    double phobicbSasa; /* hydrophobic bSASA */
    double philicbSasa; /* hydrophilic bSASA */
    double bSasa; /* bSASA */
} MolSasa;

/*___________________________________________________________________________*/
/* prototypes */
int init_sasa(Str *pdb, Type *type, MolSasa *molSasa, ConstantSasa *constant_sasa, Arg *arg);
void free_sasa(MolSasa *molSasa);
void compute_sasa(Str *pdb, Topol *topol, Type *type, MolSasa *molSasa, \
		ConstantSasa *constant_sasa, ConstantSasa *res_sasa, Arg *arg);

#endif
