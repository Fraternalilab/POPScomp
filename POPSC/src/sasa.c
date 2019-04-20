/*=============================================================================
sasa.h : SASA computation
Copyright (C) 2002-2018 Franca Fraternali
Copyright (C) 2008-2018 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin and Valerie Hindie
Read the COPYING file for license information.
=============================================================================*/

#include "config.h"
#include "sasa.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*___________________________________________________________________________*/
/** compute sphere surface */
__inline__ static double sphere_surface(double atomRadius, float rSolvent)
{
	return (4. * PI * pow((atomRadius + rSolvent), 2));
}
/*___________________________________________________________________________*/
/** compute c[ij]1 */
__inline__ static double compute_c1(double atomRadius, float rSolvent)
{
    return (PI * (atomRadius + rSolvent));
}

/*___________________________________________________________________________*/
/** compute c[ij]3 */
__inline__ static double compute_c3(double atomRadius_k, double atomRadius_l, \
	double atomDistance)
{
    return (1. + (atomRadius_k - atomRadius_l) / atomDistance);
}

/*___________________________________________________________________________*/
/** compute atomic SASA */
__inline__ static double atom_sasa(MolSasa *molSasa, int k, double connectivityParameter, \
	double bkl, double atomParameter_k)
{
	return (molSasa->atomSasa[k].sasa * (1.0 - \
		(connectivityParameter * bkl * atomParameter_k / molSasa->atomSasa[k].surface)));
}

/*___________________________________________________________________________*/
/** compute atomic bSASA */
__inline__ static double atom_bsasa(MolSasa *molSasa, int k, double connectivityParameter, \
	double bkl, double atomParameter_k)
{
	return (molSasa->atomSasa[k].sasa * \
		(connectivityParameter * bkl * atomParameter_k / molSasa->atomSasa[k].surface));
}

/*___________________________________________________________________________*/
/** initialise all SASAs */
int init_sasa(Str *pdb, Type *type, MolSasa *molSasa, ConstantSasa *constant_sasa, Arg *arg)
{
	unsigned int i;

	molSasa->atomSasa = safe_malloc(pdb->nAtom * sizeof(AtomSasa));
	molSasa->resSasa = safe_malloc(pdb->nAllResidue * sizeof(ResSasa));
	molSasa->chainSasa = safe_malloc(pdb->nChain * sizeof(ChainSasa));

	for (i = 0; i < pdb->nAtom; ++ i) {
		/*___________________________________________________________________________*/
		/* initialise atom SASAs */
		/* start value of atom SASA is surface of isolated atom */
		sphere_surface(constant_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].radius, arg->rProbe);
		molSasa->atomSasa[i].surface = \
			sphere_surface(constant_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].radius, arg->rProbe);
		molSasa->atomSasa[i].sasa = molSasa->atomSasa[i].surface;
		molSasa->atomSasa[i].nOverlap = 0; /* no overlaps yet (isolated atom) */
		molSasa->atomSasa[i].phobicbSasa = 0.; /* hydrophobic buried SASA */
		molSasa->atomSasa[i].philicbSasa = 0.; /* hydrophilic buried SASA */
		molSasa->atomSasa[i].bSasa = 0.; /* buried SASA */
	}

	/*___________________________________________________________________________*/
	/* initialise residue SASA */
	for (i = 0; i < pdb->nAllResidue; ++ i) {
		molSasa->resSasa[i].surface = 0.;
		molSasa->resSasa[i].phobicSasa = 0.;
		molSasa->resSasa[i].philicSasa = 0.;
		molSasa->resSasa[i].sasa = 0.;
		molSasa->resSasa[i].nOverlap = 0;
		molSasa->resSasa[i].atomRef = -1;
		molSasa->resSasa[i].phobicbSasa = 0.;
		molSasa->resSasa[i].philicbSasa = 0.;
		molSasa->resSasa[i].bSasa = 0.;
	}

	/*___________________________________________________________________________*/
	/* initialise chain SASA */
	for (i = 0; i < pdb->nChain; ++ i) {
		molSasa->chainSasa[i].surface = 0.;
		molSasa->chainSasa[i].phobicSasa = 0.;
		molSasa->chainSasa[i].philicSasa = 0.;
		molSasa->chainSasa[i].sasa = 0;
		molSasa->chainSasa[i].first = -1;
		molSasa->chainSasa[i].last = -1;
		molSasa->chainSasa[i].phobicbSasa = 0.,
		molSasa->chainSasa[i].philicbSasa = 0.;
		molSasa->chainSasa[i].bSasa = 0.;
	}

	/*___________________________________________________________________________*/
	/* initialise molecular SASA */
	molSasa->phobicSasa = 0.;
	molSasa->philicSasa = 0.;
	molSasa->sasa = 0.;
	molSasa->phobicbSasa = 0.,
	molSasa->philicbSasa = 0.;
	molSasa->bSasa = 0.;

	return 0;
}

/*___________________________________________________________________________*/
/** free SASAs */
void free_sasa(MolSasa *molSasa)
{
	free(molSasa->atomSasa); /* atoms */
	free(molSasa->resSasa); /* residues */
	free(molSasa->chainSasa); /* chains */
}

/*___________________________________________________________________________*/
/** atom SASA modification from one contact */
__inline__ static int mod_atom_sasa(Str *pdb, Topol *topol, Type *type, \
	MolSasa *molSasa, ConstantSasa *constant_sasa, double connectivityParameter, \
	int i, int j, float rSolvent)
{
    double atomRadius_i, atomRadius_j;
	double atomParameter_i, atomParameter_j;
    double ci1, cj1, cc2, ci3, cj3, bij, bji;
	float atomDistance = 0.;
	float cutoffRadius = 0.;

	/*___________________________________________________________________________*/
	/* safety check */   
	if (i == j) {
		fprintf(stderr, "Problematic conformation at atoms %d %d\n",
			pdb->atom[i].atomNumber, pdb->atom[i+1].atomNumber);
		/*Warning("Atom distances too short? Probably incorrect POPS results!");*/
		Error("Serious conformational problems in input structure.\nCheck the atom numbers listed above for steric clashes.");
	}

	/*___________________________________________________________________________*/
	/* initialise */   
	/* atom specific radius */
    atomRadius_i = constant_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].radius;
    atomRadius_j = constant_sasa->atomDataSasa[type->residueType[j]][type->atomType[j]].radius;
	/* atom specific parameter */
    atomParameter_i = constant_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].parameter;
    atomParameter_j = constant_sasa->atomDataSasa[type->residueType[j]][type->atomType[j]].parameter;

	/*___________________________________________________________________________*/
	/* skip overlap area computation if the two atoms i,j do not overlap */
    if ((cutoffRadius = cutoff_radius(type, constant_sasa, i, j, rSolvent)) < \
		 (atomDistance = atom_distance(pdb, i, j)))
		return(1);

	/* shortest atomic bond length is .74 A in hydrogen molecule H_2 */
	if (atomDistance < .74) {
		fprintf(stderr, "Atom distance %d %d = %f Angstrom\n",
			pdb->atom[i].atomNumber, pdb->atom[j].atomNumber, atomDistance);
		Warning("Atom distance too short! Probably incorrect POPS results!");
	}

	/*___________________________________________________________________________*/
	/* compute c[ij]1 */
	ci1 = compute_c1(atomRadius_i, rSolvent);
	cj1 = compute_c1(atomRadius_j, rSolvent);

	/* compute cc2 */
    cc2 = cutoffRadius - atomDistance;

	/* compute c[ij]3 */
    ci3 = compute_c3(atomRadius_j, atomRadius_i, atomDistance);
    cj3 = compute_c3(atomRadius_i, atomRadius_j, atomDistance);
    
	/* compute b[ij,ji] */
    bij =  ci1 * cc2 * ci3;
    bji =  cj1 * cc2 * cj3;

	/* compute c[ij]; jk: moved to atom_sasa formula */
    /*ci  = pCon * bij;
    cj  = pCon * bji;*/

	/* count overlaps */
    ++ molSasa->atomSasa[i].nOverlap;

	/* count overlaps */
    ++ molSasa->atomSasa[i].nOverlap;
    ++ molSasa->atomSasa[j].nOverlap;

	/* compute atom SASA for atoms i and j */
    molSasa->atomSasa[i].sasa = atom_sasa(molSasa, i, connectivityParameter, bij, atomParameter_i);
    molSasa->atomSasa[j].sasa = atom_sasa(molSasa, j, connectivityParameter, bji, atomParameter_j);

	/* compute atom bSASA for atoms i and j */
	/* select side-chain (including CA) atoms and
		determine polarity of neghbour (overlap) atom */
	if (pdb->atom[i].residueNumber != pdb->atom[j].residueNumber) {
		if (((type->atomType[i] == 1) || (type->atomType[i] > 3)) && (constant_sasa->atomDataSasa[type->residueType[j]][type->atomType[j]].polarity == 0)) {
			molSasa->atomSasa[i].phobicbSasa += atom_bsasa(molSasa, i, connectivityParameter, bij, atomParameter_i);
		} else if (((type->atomType[i] == 1) || (type->atomType[i] > 3)) && (constant_sasa->atomDataSasa[type->residueType[j]][type->atomType[j]].polarity == 1)) {
			molSasa->atomSasa[i].philicbSasa += atom_bsasa(molSasa, i, connectivityParameter, bij, atomParameter_i);
		} else {
			molSasa->atomSasa[i].philicbSasa += 0.;
		}

		if ((type->atomType[j] > 3) && (constant_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].polarity == 0)) {
			molSasa->atomSasa[j].phobicbSasa += atom_bsasa(molSasa, j, connectivityParameter, bij, atomParameter_j);
		} else if ((type->atomType[j] > 3) && (constant_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].polarity == 1)) {
			molSasa->atomSasa[j].philicbSasa += atom_bsasa(molSasa, j, connectivityParameter, bij, atomParameter_j);
		} else {
			molSasa->atomSasa[i].philicbSasa += 0.;
		}
	}

	molSasa->atomSasa[i].bSasa = molSasa->atomSasa[i].phobicbSasa + molSasa->atomSasa[i].philicbSasa;
	molSasa->atomSasa[j].bSasa = molSasa->atomSasa[j].phobicbSasa + molSasa->atomSasa[j].philicbSasa;

	/* record parameters: increment neighbour index */
	++ topol->neighbourPar[i][0];
	++ topol->neighbourPar[j][0];
	/* record parameter product 'p_ij * b_ij' */
	topol->neighbourPar[i][(int)topol->neighbourPar[i][0]] = connectivityParameter * bij;
	topol->neighbourPar[j][(int)topol->neighbourPar[j][0]] = connectivityParameter * bji;

	return(0);
}

/*___________________________________________________________________________*/
/** atomic SASA calculation:
	modify initial atom SASA (=total surface) for each atom interaction; */
/** atomic bSASA calculation (last subroutine):
	compute buried SASA due to neighbour atoms */
static int compute_atom_sasa(Str *pdb, Topol *topol, Type *type, MolSasa *molSasa, \
	ConstantSasa *constant_sasa, Arg *arg)
{
	unsigned int i;
   
	/*___________________________________________________________________________*/
	/** 1-2 interactions along bonds 
		using 1-2 connectivity parameters and identifiers of bond-forming atoms */
	for (i = 0; i < topol->nBond; ++ i) {
        mod_atom_sasa(pdb, topol, type, molSasa, constant_sasa,
			constant_sasa->connect_12_parameter,
			topol->ib[i], topol->jb[i], arg->rProbe);
	}

	/*___________________________________________________________________________*/
    /** 1-3 interactions along angles */
    for (i = 0; i < topol->nAngle; ++ i) {
        mod_atom_sasa(pdb, topol, type, molSasa, constant_sasa,
			constant_sasa->connect_13_parameter,
			topol->it[i], topol->kt[i], arg->rProbe);
	}

	/*___________________________________________________________________________*/
    /* 1-4 interactions along torsion angles */
    for (i = 0; i < topol->nTorsion; ++ i) {
		mod_atom_sasa(pdb, topol, type, molSasa, constant_sasa,
			constant_sasa->connect_14_parameter,
			topol->ip[i], topol->lp[i], arg->rProbe);
	}

	/*___________________________________________________________________________*/
    /*  >(1-4) interactions */
    for (i = 0; i < topol->nNonBonded; ++ i) {
        mod_atom_sasa(pdb, topol, type, molSasa, constant_sasa,
			constant_sasa->connect_15_parameter,
			topol->in[i], topol->jn[i], arg->rProbe);
	}
	return(0);
}

/*___________________________________________________________________________*/
/** residuic and molecular SASA calculation */
static int compute_res_chain_mol_sasa(Str *pdb, Type *type, MolSasa *molSasa, \
	ConstantSasa *constant_sasa, ConstantSasa *res_sasa)
{
    unsigned int i, j, k;

	/*___________________________________________________________________________*/
	/* record first residue number in PDB file */
	/* atom '0' is reference for first residue */
	molSasa->resSasa[0].atomRef = 0;
	molSasa->chainSasa[0].first = 0;

    for (i = 0, j = 0, k = 0; i < pdb->nAtom; ++ i) { 
		/*___________________________________________________________________________*/
		/* first (==0) residue index */
		/* first atom of each residue is reference for residue type */
		if (i == 0) {
			molSasa->resSasa[j].atomRef = i;
			molSasa->resSasa[j].surface = res_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].surface;
		}
		/*___________________________________________________________________________*/
		/* increment residue index */
		//if (i > 0 && (pdb->atom[i].residueNumber != pdb->atom[i - 1].residueNumber)) { 
		if (i > 0 && (pdb->atom[i].residueNumber != pdb->atom[i - 1].residueNumber ||
				      strcmp(pdb->atom[i].icode, pdb->atom[i - 1].icode) != 0)) {
			++ j;
			molSasa->resSasa[j].atomRef = i; /* assign atom reference */
			/* set reference residue surface*/
			molSasa->resSasa[j].surface = res_sasa->atomDataSasa[type->residueType[i]][type->atomType[i]].surface;
		}
		/*___________________________________________________________________________*/
		/* increment chain index */
        if (i > 0 && pdb->atom[i].chainIdentifier[0] != pdb->atom[i - 1].chainIdentifier[0]) {
			++ k;
			molSasa->chainSasa[k-1].last = i-1;
			molSasa->chainSasa[k].first = i;
		}

		/*___________________________________________________________________________*/
		/* sum atomic SASA to residue, chain and molecule SASA */
		if (constant_sasa[0].atomDataSasa[type->residueType[i]][type->atomType[i]].polarity == 0) {
			molSasa->resSasa[j].phobicSasa += molSasa->atomSasa[i].sasa;
			molSasa->chainSasa[k].phobicSasa += molSasa->atomSasa[i].sasa;
			molSasa->phobicSasa += molSasa->atomSasa[i].sasa;
		}
		else {
			molSasa->resSasa[j].philicSasa += molSasa->atomSasa[i].sasa;
			molSasa->chainSasa[k].philicSasa += molSasa->atomSasa[i].sasa;
			molSasa->philicSasa += molSasa->atomSasa[i].sasa;
		}

		molSasa->resSasa[j].sasa += molSasa->atomSasa[i].sasa;
		molSasa->chainSasa[k].sasa += molSasa->atomSasa[i].sasa;
		molSasa->sasa += molSasa->atomSasa[i].sasa;

		/* sum number of atom overlaps to residuic number of overlaps */
		molSasa->resSasa[j].nOverlap += molSasa->atomSasa[i].nOverlap;


		/*___________________________________________________________________________*/
		/* sum atomic bSASA to residue, chain and molecule SASA */
		molSasa->resSasa[j].phobicbSasa += molSasa->atomSasa[i].phobicbSasa;
		molSasa->resSasa[j].philicbSasa += molSasa->atomSasa[i].philicbSasa;
		molSasa->resSasa[j].bSasa += molSasa->atomSasa[i].phobicbSasa + molSasa->atomSasa[i].philicbSasa;

		molSasa->chainSasa[k].phobicbSasa += molSasa->atomSasa[i].phobicbSasa;
		molSasa->chainSasa[k].philicbSasa += molSasa->atomSasa[i].philicbSasa;
		molSasa->chainSasa[k].bSasa += molSasa->atomSasa[i].phobicbSasa + molSasa->atomSasa[i].philicbSasa;

		molSasa->phobicbSasa += molSasa->atomSasa[i].phobicbSasa;
		molSasa->philicbSasa += molSasa->atomSasa[i].philicbSasa;
		molSasa->bSasa += molSasa->atomSasa[i].phobicbSasa + molSasa->atomSasa[i].philicbSasa;
	}
	molSasa->chainSasa[k].last = i - 1;

	return(0);
}

/*___________________________________________________________________________*/
/** compute SASAs */
void compute_sasa(Str *pdb, Topol *topol, Type *type, MolSasa *molSasa, \
	ConstantSasa *constant_sasa, ConstantSasa *res_sasa, Arg *arg)
{
    compute_atom_sasa(pdb, topol, type, molSasa, constant_sasa, arg); /* compute SASA per atom */
	compute_res_chain_mol_sasa(pdb, type, molSasa, constant_sasa, res_sasa); /* sum up atom SASAs */
}

