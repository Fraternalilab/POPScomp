/*==============================================================================
putsasa.c : routines for printing SASA
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "putsasa.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/** print molecule composition */
static void print_composition(FILE *sasaOutFile, Arg *arg, Argpdb *argpdb, Str *pdb)
{
	if (! arg->noHeaderOut) fprintf(sasaOutFile, "\n=== COMPOSITION ===\n");
	if (! argpdb->coarse)
		fprintf(sasaOutFile,\
			"\nProtein %8s\n%8d chains\n%8d standard residues\n"
			"%8d total residues (standard residues + HETATM residues lacking CA or P))\n"
			"%8d atoms (excluding hydrogen atoms)\n",
			arg->pdbInFileName, pdb->nChain, pdb->nResidue, pdb->nAllResidue, pdb->nAtom);
	else
		fprintf(sasaOutFile,\
			"\nProtein %8s\n%8d chains\n%8d standard residues\n"
			"%8d total residues (standard residues + HETATM residues lacking CA or P)\n"
			"%8d atoms (C-alpha and P atoms)\n",
			arg->pdbInFileName, pdb->nChain, pdb->nResidue, pdb->nAllResidue, pdb->nAtom);
}

/*____________________________________________________________________________*/
/** print molecule topology */
static void print_topology(FILE *sasaOutFile, Arg *arg, Topol *topol)
{
	if (! arg->noHeaderOut) fprintf(sasaOutFile, "\n=== TOPOLOGY ===\n");
	fprintf(sasaOutFile, "\n%8d bonds (1,2-interactions)\n"
						 "%8d angles (1,3-interactions)\n"
						 "%8d torsions (1,4-interactions excluding rings)\n"
						 "%8d non-bonded (1,>4-interactions)\n",
		topol->nBond, topol->nAngle, topol->nTorsion, topol->nNonBonded);

}

/*____________________________________________________________________________*/
/** print residue and atom types */
static void print_types(FILE *sasaOutFile, Arg *arg, Type *type, ConstantSasa *constant_sasa)
{
	unsigned int i, j;

	if (! arg->noHeaderOut) {
		fprintf(sasaOutFile, "\n=== TYPES ===\n");
		fprintf(sasaOutFile, "\nResidue\tAtom\tRadius\tParameter\tPolarity\n");
	}

	for (i = 0; i < constant_sasa->nResidueType; ++ i)
		for (j = 0; j < constant_sasa->nAtomResidue[i]; ++ j)
			fprintf(sasaOutFile, "%s\t%s\t%8.6f\t%8.6f\t%1d\n",
				constant_sasa->atomDataSasa[i][j].residueName,
				constant_sasa->atomDataSasa[i][j].atomName,
				constant_sasa->atomDataSasa[i][j].radius,
				constant_sasa->atomDataSasa[i][j].parameter,
				constant_sasa->atomDataSasa[i][j].polarity);
}

/*___________________________________________________________________________*/
/** print atom SASA  */
static void print_atom_sasa(FILE *sasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
	unsigned int i, j;
	float surface_ratio = 0.;
	char *atomName;
	char *residueName;

	if (! arg->noHeaderOut && ! arg->rout) {
		fprintf(sasaOutFile, "\n=== ATOM SASAs ===\n");
	}
	if (! arg->noHeaderOut || arg->rout) {
		fprintf(sasaOutFile, "AtomNr\tAtomNe\tResidNe\tChain\tResidNr\tiCode\tSASA/A^2\tQ(SASA)\tN(overl)\tAtomTp\tAtomGp\tSurf/A^2\n");
	}

	/* before the first line:
		print dummy lines to substitute missing hydrogen atom lines */
	if (arg->padding)  {
		j = 0;
		while (++ j < pdb->atom[0].atomNumber) {
			fprintf(sasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.2f\t%8d\t\t%2d\t\t%2d\t%10.2f\n",
				j,
				"XXX",
				pdb->atom[0].residueName,
				" ",
				pdb->atom[0].residueNumber,
				" ",
				0.,
				0.,
				0,
				0,
				0,
				0.);
		}
	}

	for (i = 0; i < pdb->nAtom; ++ i) {
		if (molSasa->atomSasa[i].surface > 0.) {
			surface_ratio = molSasa->atomSasa[i].sasa / molSasa->atomSasa[i].surface;
			/*assert(surface_ratio >= 0. && surface_ratio <= 1.);*/
		} else {
#ifdef NAN
			surface_ratio = NAN;
#else
			surface_ratio = 0./0.;
#endif
		}

		/* use original residue number of heteroresidues */
		if (pdb->atom[i].het) {
			atomName = &(pdb->atom[i].atomNameHet[0]);
			residueName = &(pdb->atom[i].residueNameHet[0]);
		} else {
			atomName = &(pdb->atom[i].atomName[0]);
			residueName = &(pdb->atom[i].residueName[0]);
		}

		fprintf(sasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.4f\t%8d\t\t%2d\t\t%2d\t%10.2f\n",
			pdb->atom[i].atomNumber,
			atomName,
			residueName,
			pdb->atom[i].chainIdentifier,
			pdb->atom[i].residueNumber,
			pdb->atom[i].icode,
			molSasa->atomSasa[i].sasa,
			surface_ratio,
			molSasa->atomSasa[i].nOverlap,
			pdb->atom[i].atomType,
			pdb->atom[i].groupID,
			molSasa->atomSasa[i].surface);
		/* after the first atom :
			print dummy lines to substitute missing hydrogen atom lines */
		if (arg->padding && ((i + 1) < pdb->nAtom)) {
			j = pdb->atom[i].atomNumber;
			while (++ j < pdb->atom[i+1].atomNumber) {
				fprintf(sasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.4f\t%8d\t\t%2d\t\t%2d%10.2f\n",
					j,
					"HXX",
					pdb->atom[i].residueName,
					" ",
					pdb->atom[i].residueNumber,
					" ",
					0.,
					0.,
					0,
					0,
					0,
					0.);
			}
		}
	}
}

/*___________________________________________________________________________*/
/** print residue SASA */
static void print_residue_sasa(FILE *sasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
    unsigned int i;
	float surface_ratio = 0.;

	if (! arg->noHeaderOut && ! arg->rout) {
		fprintf(sasaOutFile, "\n=== RESIDUE SASAs ===\n");
	}
	if (! arg->noHeaderOut || arg->rout) {
		fprintf(sasaOutFile, "ResidNe\tChain\tResidNr\tiCode\tPhob/A^2\t\tPhil/A^2\tTotal/A^2\t\tQ(SASA)\tN(overl)\tSurf/A^2\n");
	}

    for (i = 0; i < pdb->nAllResidue; ++ i) { 
		if (molSasa->resSasa[i].surface > 0.) {
			surface_ratio = molSasa->resSasa[i].sasa / molSasa->resSasa[i].surface;
			/* upper limit: there can be configurational effects leading to >1. */
			/*assert(surface_ratio >= 0 && surface_ratio <= 1.5);*/
		} else {
#ifdef NAN
			surface_ratio = NAN;
#else
			surface_ratio = 0./0.;
#endif
		}
			
		fprintf(sasaOutFile, "%8s\t%3s\t%8d\t%1s\t%10.2f\t%10.2f\t%10.2f\t%10.4f\t%8d\t%10.2f\n",
			pdb->atom[molSasa->resSasa[i].atomRef].residueName,
			pdb->atom[molSasa->resSasa[i].atomRef].chainIdentifier,
			pdb->atom[molSasa->resSasa[i].atomRef].residueNumber,
			pdb->atom[molSasa->resSasa[i].atomRef].icode,
			molSasa->resSasa[i].phobicSasa,
			molSasa->resSasa[i].philicSasa,
			molSasa->resSasa[i].sasa,
			surface_ratio,
			molSasa->resSasa[i].nOverlap,
			molSasa->resSasa[i].surface);
    }
}

/*___________________________________________________________________________*/
/** print chain SASA */
static void print_chain_sasa(FILE *sasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
    unsigned int i;

	if (! arg->noHeaderOut && ! arg->rout) {
		fprintf(sasaOutFile, "\n=== CHAIN SASAs ===\n(Atom Range excluding hydrogen atoms)\n");
	}
	if (! arg->noHeaderOut || arg->rout) {
		fprintf(sasaOutFile, "Chain\tId\tAtomRange\tResidRange\t\tPhob/A^2\t\tPhil/A^2\t\tTotal/A^2\n");
	}

    for (i = 0; i < pdb->nChain; ++ i)
		fprintf(sasaOutFile, "%3d\t%3s\t%6d->%-6d\t%5d->%-5d\t%10.2f\t%10.2f\t%10.2f\n",
			i,
			pdb->atom[molSasa->chainSasa[i].first].chainIdentifier,
			pdb->atom[molSasa->chainSasa[i].first].atomNumber,
			pdb->atom[molSasa->chainSasa[i].last].atomNumber,
			pdb->atom[molSasa->chainSasa[i].first].residueNumber,
			pdb->atom[molSasa->chainSasa[i].last].residueNumber,
			molSasa->chainSasa[i].phobicSasa,
			molSasa->chainSasa[i].philicSasa,
			molSasa->chainSasa[i].sasa);
}

/*____________________________________________________________________________*/
/** total (molecule) SASA */
void print_mol_sasa(FILE *sasaOutFile, Arg *arg, MolSasa *molSasa)
{
	if (! arg->noHeaderOut && ! arg->rout) {
		fprintf(sasaOutFile, "\n=== MOLECULE SASAs ===\n");
	}
	if (! arg->noHeaderOut || arg->rout) {
		fprintf(sasaOutFile, "Phob/A^2\t\tPhil/A^2\t\tTotal/A^2\n");
	}

    fprintf(sasaOutFile, "%10.2f\t%10.2f\t%10.2f\n",
			molSasa->phobicSasa,
			molSasa->philicSasa,
			molSasa->sasa);
}

/*____________________________________________________________________________*/
/** print atom neighbour list */
void print_neighbour_list(FILE *neighbourOutFile, Arg *arg, Str *pdb, Topol *topol)
{
	unsigned int i, j;
	if (! arg->noHeaderOut) fprintf(neighbourOutFile, "\n=== ATOM NEIGHBOUR LIST ===\n\n");
	for (i = 0; i < pdb->nAtom; ++ i) {
		/* old format */
		/*
		fprintf(neighbourOutFile, "atom %d, N(n) %d, ", 
			i, topol->neighbourState[i][0]);
		for (j = 1; j <= topol->neighbourState[i][0]; ++ j)
			fprintf(neighbourOutFile, "%d ", topol->neighbourState[i][j]);
		fprintf(neighbourOutFile, "\n");
		*/
		/* new format */
		fprintf(neighbourOutFile, "%d:%s\t%d\t", 
			pdb->atom[i].atomNumber, pdb->atom[i].chainIdentifier, topol->neighbourState[i][0]);
		for (j = 1; j <= topol->neighbourState[i][0]; ++ j)
			fprintf(neighbourOutFile, "%d:%s ",
				pdb->atom[topol->neighbourState[i][j]].atomNumber,
				pdb->atom[topol->neighbourState[i][j]].chainIdentifier);
		fprintf(neighbourOutFile, "\n");
	}
}

/*____________________________________________________________________________*/
/** print neighbour parameter */
void print_neighbour_parameter(FILE *parameterOutFile, Str *pdb, Type *type, \
	Topol *topol, MolSasa *molSasa)
{
	unsigned int i, j, k;
	float dummy = 0.;
	for (i = 0; i < pdb->nAtom; ++ i) {
		 /* central atom */
		fprintf(parameterOutFile, "%6d %02d%02d %10.4f %6.0f ", 
			pdb->atom[i].atomNumber, type->residueType[i], type->atomType[i],
			molSasa->atomSasa[i].surface, topol->neighbourPar[i][0]);
		/* neighbours */
		for (j = 1; j <= topol->neighbourPar[i][0]; ++ j)
			fprintf(parameterOutFile, "%10.4f ", topol->neighbourPar[i][j]);
		/* make all rows equally long */
		for (k = j; k < 256; ++ k)
			fprintf(parameterOutFile, "%2.0f ", dummy);
		fprintf(parameterOutFile, "\n");
	}
}

/*____________________________________________________________________________*/
/** print interface residue pairs */
/* nearest neighbour atom pairs on different chains */
void print_interface(FILE *interfaceOutFile, Arg *arg, Str *pdb, Type *type, Topol *topol)
{
	unsigned int i, j;

	for (i = 0; i < pdb->nAtom; ++ i) {
		j = topol->interfaceNn[i];

		if (strcmp(pdb->atom[i].chainIdentifier, pdb->atom[j].chainIdentifier) != 0) {
			fprintf(interfaceOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.4f\n",
				pdb->atom[i].atomNumber,
				pdb->atom[i].atomName,
				pdb->atom[i].residueName,
				pdb->atom[i].chainIdentifier,
				pdb->atom[i].residueNumber,
				pdb->atom[i].icode,
				pdb->atom[j].atomNumber,
				pdb->atom[j].atomName,
				pdb->atom[j].residueName,
				pdb->atom[j].chainIdentifier,
				pdb->atom[j].residueNumber,
				pdb->atom[j].icode,
				topol->interfaceNnDist[i]);
		}
	}
}

/*____________________________________________________________________________*/
/** print SASAs */
void print_sasa(Arg *arg, Argpdb *argpdb, Str *pdb, Type *type, Topol *topol, \
				MolSasa *molSasa, ConstantSasa *constant_sasa, int frame)
{
	char sasatrajOutFileName[256];
	char rpopsOutFileName[256];
	FILE *rpopsOutFile;

	/* for single (reference) molecule */
	if (frame < 0) {
		if (! arg->silent) {
			fprintf(stdout, "\tSASA of reference molecule: %s\n", arg->sasaOutFileName);
			if (arg->rout) {
				arg->sasaOutFile = NULL;
			} else {
				arg->sasaOutFile = safe_open(arg->sasaOutFileName, "w");
			}
		}
	} else {
			sprintf(&(sasatrajOutFileName[0]), "%s.%d.%s", arg->sasatrajOutFileName, frame, "out");
			if (! arg->rout) {
				arg->sasaOutFile = NULL;
			} else {
				arg->sasaOutFile = safe_open(sasatrajOutFileName, "w");
			}
	}
	/* composition */
	if (arg->compositionOut) {
		print_composition(arg->sasaOutFile, arg, argpdb, pdb);
	}

	/* topology */
	if (arg->topologyOut)
		print_topology(arg->sasaOutFile, arg, topol);

	/* residue/atom types */
	if (arg->typeOut)
		print_types(arg->sasaOutFile, arg, type, constant_sasa);

	/* atom SASA */
	if (arg->atomOut && ! argpdb->coarse) {
		if (arg->rout) {
			sprintf(rpopsOutFileName, "%s.%s", arg->sasaOutFileName, "rpopsAtom");
			rpopsOutFile = safe_open(rpopsOutFileName, "w");
			print_atom_sasa(rpopsOutFile, arg, pdb, molSasa);
			fclose(rpopsOutFile);
		} else {
			print_atom_sasa(arg->sasaOutFile, arg, pdb, molSasa);
		}
	}

	/* residue SASA */
	if (arg->residueOut) {
		if (arg->rout) {
			sprintf(rpopsOutFileName, "%s.%s", arg->sasaOutFileName, "rpopsResidue");
			rpopsOutFile = safe_open(rpopsOutFileName, "w");
			print_residue_sasa(rpopsOutFile, arg, pdb, molSasa);
			fclose(rpopsOutFile);
		} else {
			print_residue_sasa(arg->sasaOutFile, arg, pdb, molSasa);
		}
	}

	/* chain SASA */
	if (arg->chainOut) {
		if (arg->rout) {
			sprintf(rpopsOutFileName, "%s.%s", arg->sasaOutFileName, "rpopsChain");
			rpopsOutFile = safe_open(rpopsOutFileName, "w");
			print_chain_sasa(rpopsOutFile, arg, pdb, molSasa);
			fclose(rpopsOutFile);
		} else {
			print_chain_sasa(arg->sasaOutFile, arg, pdb, molSasa);
		}
	}

	/* molecule SASA */
	if (! arg->noTotalOut) {
		if (arg->rout) {
			sprintf(rpopsOutFileName, "%s.%s", arg->sasaOutFileName, "rpopsMolecule");
			rpopsOutFile = safe_open(rpopsOutFileName, "w");
			print_mol_sasa(rpopsOutFile, arg, molSasa);
			fclose(rpopsOutFile);
		} else {
			print_mol_sasa(arg->sasaOutFile, arg, molSasa);
		}
	}

	if (! arg->rout) {
		fclose(arg->sasaOutFile);
	}

	/* neighbour list */
	if (arg->neighbourOut) {
		if (frame < 0)
			arg->neighbourOutFile = safe_open(arg->neighbourOutFileName, "w");
		else
			arg->neighbourOutFile = safe_open(arg->neighbourOutFileName, "a");

		print_neighbour_list(arg->neighbourOutFile, arg, pdb, topol);
		fclose(arg->neighbourOutFile);
	}

	/* neighbour parameters (for benchmarking) */
	if (arg->parameterOut) {
		if (frame < 0)
			arg->parameterOutFile = safe_open(arg->parameterOutFileName, "w");
		else
			arg->parameterOutFile = safe_open(arg->parameterOutFileName, "a");

		print_neighbour_parameter(arg->parameterOutFile, pdb, type, topol, molSasa);
		fclose(arg->parameterOutFile);
	}

	/* interface residue pairs */
	if (arg->interfaceOut) {
		if (frame < 0)
			arg->interfaceOutFile = safe_open(arg->interfaceOutFileName, "w");
		else
			arg->interfaceOutFile = safe_open(arg->interfaceOutFileName, "a");

		print_interface(arg->interfaceOutFile, arg, pdb, type, topol);
		fclose(arg->interfaceOutFile);
	}
}

