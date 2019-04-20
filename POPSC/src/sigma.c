/*=============================================================================
sigma.c : Solvation Free Energy (SFE) computation
Copyright (C) 2011-2017 Franca Fraternali
Copyright (C) 2011-2017 Jens Kleinjung
Read the COPYING file for license information.
=============================================================================*/

#include "config.h"
#include "sigma.h"

/*___________________________________________________________________________*/
/** initialise all SFEs */
int init_sfe(Str *pdb, Type *type, MolSFE *molSFE, ConstantSigma *constant_sigma, Arg *arg)
{
	unsigned int i;

	molSFE->atomSFE = safe_malloc(pdb->nAtom * sizeof(AtomSFE));
	molSFE->resSFE = safe_malloc(pdb->nAllResidue * sizeof(ResSFE));
	molSFE->chainSFE = safe_malloc(pdb->nChain * sizeof(ChainSFE));

	for (i = 0; i < pdb->nAtom; ++ i) {
		/*___________________________________________________________________________*/
		/* initialise atom SFEs */
		molSFE->atomSFE[i].sfe_type = 0.;
		molSFE->atomSFE[i].sfe_group = 0.;
/*
#if DEBUG>1
		fprintf(stderr, "%s:%d: atomSFE %d\t%f\t%f\n",
			__FILE__, __LINE__, i,
			molSFE->atomSFE[i].sfe_type,
			molSFE->atomSFE[i].sfe_group);
#endif
*/
	}

	/*___________________________________________________________________________*/
	/* initialise residue SFEs */
	for (i = 0; i < pdb->nAllResidue; ++ i) {
			molSFE->resSFE[i].sfe_type = 0.;
			molSFE->resSFE[i].sfe_group = 0.;
			molSFE->resSFE[i].atomRef = -1;
/*
#if DEBUG>1
		fprintf(stderr, "%s:%d: resSFE %d\t%f\t%f\n",
			__FILE__, __LINE__, i,
			molSFE->resSFE[i].sfe_type,
			molSFE->resSFE[i].sfe_group);
#endif
*/
	}

	/*___________________________________________________________________________*/
	/* initialise chain SFEs */
	for (i = 0; i < pdb->nChain; ++ i) {
			molSFE->chainSFE[i].sfe_type = 0.;
			molSFE->chainSFE[i].sfe_group = 0.;
/*
#if DEBUG>1
		fprintf(stderr, "%s:%d: chainSFE %d\t%f\t%f\n",
			__FILE__, __LINE__, i,
			molSFE->chainSFE[i].sfe_type,
			molSFE->chainSFE[i].sfe_group);
#endif
*/
	}

	/*___________________________________________________________________________*/
	/* initialise molecular SFE */
	molSFE->sfe_type = 0.;
	molSFE->sfe_group = 0.;

	return 0;
}

/*___________________________________________________________________________*/
/** free SFEs */
void free_sfe(MolSFE *molSFE)
{
	free(molSFE->atomSFE); /* atoms */
	free(molSFE->resSFE); /* residues */
	free(molSFE->chainSFE); /* chains */
}

/*___________________________________________________________________________*/
/** atomic SFE calculation */
/* scale SASA from Angstrom to nm (factor 1/100) */
static int compute_atom_sfe(Str *pdb, Type *type, MolSasa *molSasa,\
	MolSFE *molSFE, ConstantSigma *constant_sigma, Arg *arg)
{
	unsigned int i;
   
	for (i = 0; i < pdb->nAtom; ++ i) {
		molSFE->atomSFE[i].sfe_type = molSasa->atomSasa[i].sasa / 100 * \
				constant_sigma->atomDataSigma[type->residueType[i]][type->atomType[i]].sigma_type;
		molSFE->atomSFE[i].sfe_group = molSasa->atomSasa[i].sasa / 100 * \
				constant_sigma->atomDataSigma[type->residueType[i]][type->atomType[i]].sigma_group;
	}

	return(0);
}

/*___________________________________________________________________________*/
/** residuic and molecular SFE calculation */
static int compute_res_chain_mol_sfe(Str *pdb, Type *type, MolSFE *molSFE)
{
    unsigned int i, j, k;

	/*___________________________________________________________________________*/
	/* record first residue number in PDB file */
	/* atom '0' is reference for first residue */
	molSFE->resSFE[0].atomRef = 0;
	molSFE->chainSFE[0].first = 0;

    for (i = 0, j = 0, k = 0; i < pdb->nAtom; ++ i) { 
		/*___________________________________________________________________________*/
		/* first (==0) residue index */
		/* first atom of each residue is reference for residue type */
		if (i == 0) {
			molSFE->resSFE[j].atomRef = i;
		}
		/*___________________________________________________________________________*/
		/* increment residue index */
		/* first atom of each residue is reference for residue type */
		if (i > 0 && (pdb->atom[i].residueNumber != pdb->atom[i - 1].residueNumber || 
					 strcmp(pdb->atom[i].icode, pdb->atom[i - 1].icode) != 0)) {
			++ j; /* increment residue index */
			molSFE->resSFE[j].atomRef = i; /* assign atom reference */
		}

		/*___________________________________________________________________________*/
		/* increment chain index */
        if (i > 0 && pdb->atom[i].chainIdentifier[0] != pdb->atom[i - 1].chainIdentifier[0]) {
			++ k;
			molSFE->chainSFE[k-1].last = i-1;
			molSFE->chainSFE[k].first = i;
		}

		/*___________________________________________________________________________*/
		/* sum atomic SFE to residue, chain and molecule SFE */
		/* sum over atom type sigmas */
		molSFE->resSFE[j].sfe_type += molSFE->atomSFE[i].sfe_type;
		molSFE->chainSFE[k].sfe_type += molSFE->atomSFE[i].sfe_type;
		molSFE->sfe_type += molSFE->atomSFE[i].sfe_type;
		/* sum over atom group sigmas */
		molSFE->resSFE[j].sfe_group += molSFE->atomSFE[i].sfe_group;
		molSFE->chainSFE[k].sfe_group += molSFE->atomSFE[i].sfe_group;
		molSFE->sfe_group += molSFE->atomSFE[i].sfe_group;

	}
	molSFE->chainSFE[k].last = i - 1;

	return(0);
}

/*___________________________________________________________________________*/
/** compute SFEs */
void compute_sfe(Str *pdb, Type *type, MolSasa *molSasa, \
	MolSFE *molSFE, ConstantSigma *constant_sigma, Arg *arg)
{
	compute_atom_sfe(pdb, type, molSasa, molSFE, constant_sigma, arg); /* compute SFE per atom */
	compute_res_chain_mol_sfe(pdb, type, molSFE); /* sum up atom SFEs */
}

