/*==============================================================================
sigma.h : Solvent Free Energy (SFE) computation
Copyright (C) 2011 Franca Fraternali
Copyright (C) 2011 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef SIGMA_H
#define SIGMA_H

#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "topol.h"
#include "type.h"
#include "sasa.h"
#include "sigma_const.h"
#include "vector.h"

/*___________________________________________________________________________*/
/* structures */
/* Atomic solvation free energy */
typedef struct  
{
    double sfe_type;
	double sfe_group;
} AtomSFE;

/* Residuic solvation free energy */
typedef struct  
{
    double sfe_type;
	double sfe_group;
	int atomRef; /* first atom of residue to refer to residue number and name */
} ResSFE;

/* Chainic solvation free energy */
typedef struct  
{
    double sfe_type;
	double sfe_group;
	int first, last; /* first/last atom of chain */
} ChainSFE;

/* Molecular solvation free energy */
typedef struct  
{
	AtomSFE *atomSFE; /* atomic values */
	ResSFE *resSFE; /* molecular values */
	ChainSFE *chainSFE; /* chain values */
    double sfe_type;
	double sfe_group;
} MolSFE;

/*___________________________________________________________________________*/
/* prototypes */
int init_sfe(Str *pdb, Type *type, MolSFE *molSigma, ConstantSigma *constant_sigma, Arg *arg);
void free_sfe(MolSFE *molSFE);
void compute_sfe(Str *pdb, Type *type, MolSasa *molSasa, MolSFE *molSFE, ConstantSigma *constant_sigma, Arg *arg);

#endif
