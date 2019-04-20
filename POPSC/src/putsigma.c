/*==============================================================================
putsigma.c : routines for printing sigma and SFE valus
Copyright (C) 2011 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "putsigma.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*___________________________________________________________________________*/
/** print atom SFE */
static void print_atom_sfe(FILE *sigmaOutFile, Arg *arg, Str *pdb, MolSFE *molSFE)
{
	unsigned int i, j;

	if (! arg->noHeaderOut) {
		fprintf(sigmaOutFile, "\n=== ATOM Solvation Free Energy ===\n");
		fprintf(sigmaOutFile, "\nAtomNr\tAtomNe\tResiNe\tChain\tResidNr\tSFEt/(kJ/mol)\tSFEg/(kJ/mol)\tAtom Type\tAtom Group\n");
	}

	/* before the first line:
		print dummy lines to substitute missing hydrogen atom lines */
	if (arg->padding)  {
		j = 0;
		while (++ j < pdb->atom[0].atomNumber) {
			fprintf(sigmaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%10.2f\t\t%10.2f\t\t\t%2d\t\t%2d\n",
				j, "XXX", pdb->atom[0].residueName, " ",
				pdb->atom[0].residueNumber,
				0., 0., 0, 0);
		}
	}

	for (i = 0; i < pdb->nAtom; ++ i) {
		fprintf(sigmaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t\t%10.2f\t\t\t%2d\t\t%2d\n",
			pdb->atom[i].atomNumber, pdb->atom[i].atomName,
			pdb->atom[i].residueName, pdb->atom[i].chainIdentifier,
			pdb->atom[i].residueNumber,
			pdb->atom[i].icode,
			molSFE->atomSFE[i].sfe_type, molSFE->atomSFE[i].sfe_group,
			pdb->atom[i].atomType, pdb->atom[i].groupID);
		/* after the first atom :
			print dummy lines to substitute missing hydrogen atom lines */
		if (arg->padding && ((i + 1) < pdb->nAtom)) {
			j = pdb->atom[i].atomNumber;
			while (++ j < pdb->atom[i+1].atomNumber) {
				fprintf(sigmaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t\t%10.2f\t\t\t%2d\t\t%2d\n",
					j, "HXX", pdb->atom[i].residueName, " ",
					pdb->atom[i].residueNumber,
					pdb->atom[i].icode,
					0., 0., 0, 0);
			}
		}
	}
}

/*___________________________________________________________________________*/
/** print residue SFE */
static void print_residue_sfe(FILE *sigmaOutFile, Arg *arg, Str *pdb, MolSFE *molSFE)
{
    unsigned int i;

	if (! arg->noHeaderOut) {
		fprintf(sigmaOutFile, "\n=== RESIDUE Solvation Free Energy ===\n");
		fprintf(sigmaOutFile, "\nResid\tChain\tResidNr\tSFEt/(kJ/mol)\tSFEg/(kJ/mol)\n");
	}

    for (i = 0; i < pdb->nAllResidue; ++ i) { 
		fprintf(sigmaOutFile, "%3s\t%3s\t%8d\t%1s\t%10.2f\t\t%10.2f\n",
			pdb->atom[molSFE->resSFE[i].atomRef].residueName,
			pdb->atom[molSFE->resSFE[i].atomRef].chainIdentifier,
			pdb->atom[molSFE->resSFE[i].atomRef].residueNumber,
			pdb->atom[molSFE->resSFE[i].atomRef].icode,
			molSFE->resSFE[i].sfe_type,
			molSFE->resSFE[i].sfe_group);
    }
}

/*___________________________________________________________________________*/
/** print chain SFE */
static void print_chain_sfe(FILE *sigmaOutFile, Arg *arg, Str *pdb, MolSFE *molSFE)
{
    unsigned int i;

	if (! arg->noHeaderOut) {
		fprintf(sigmaOutFile, "\n=== CHAIN Solvation Free Energy ===\n(Atom Range excluding hydrogen atoms)\n");
		fprintf(sigmaOutFile, "\nChain\tId\tAtom Range\tResidue Range\t\tSFEt/(kJ/mol)\tSFEg/(kJ/mol)\n");
	}

    for (i = 0; i < pdb->nChain; ++ i)
		fprintf(sigmaOutFile, "%3d\t%3s\t%6d->%-6d\t%5d->%-5d\t%10.2f\t\t%10.2f\n",
			i,
			pdb->atom[molSFE->chainSFE[i].first].chainIdentifier,
			pdb->atom[molSFE->chainSFE[i].first].atomNumber,
			pdb->atom[molSFE->chainSFE[i].last].atomNumber,
			pdb->atom[molSFE->chainSFE[i].first].residueNumber,
			pdb->atom[molSFE->chainSFE[i].last].residueNumber,
			molSFE->chainSFE[i].sfe_type,
			molSFE->chainSFE[i].sfe_group);
}

/*____________________________________________________________________________*/
/** total (molecule) SFE */
void print_mol_sfe(FILE *sigmaOutFile, Arg *arg, MolSFE *molSFE)
{
	if (! arg->noHeaderOut) fprintf(sigmaOutFile, "\n=== MOLECULE Solvation Free Energy ===\n");
    fprintf(sigmaOutFile, "\nSFEt: %10.2f\n"
			"SFEg: %10.2f\n\n",
			molSFE->sfe_type,
			molSFE->sfe_group);
}

/*____________________________________________________________________________*/
/** print SFEs */
void print_sfe(Arg *arg, Argpdb *argpdb, Str *pdb, Type *type, Topol *topol, \
				MolSFE *molSFE, ConstantSigma *constant_sigma, int frame)
{
	/* for single (reference) molecule */
	if (frame < 0) {
		if (! arg->silent)
			fprintf(stdout, "\tSFE of input molecule: %s\n", arg->sigmaOutFileName);
		arg->sigmaOutFile = safe_open(arg->sigmaOutFileName, "w");

		/* atom SFE */
		if (arg->atomOut)
			print_atom_sfe(arg->sigmaOutFile, arg, pdb, molSFE);

		/* residue SFE */
		if (arg->residueOut)
			print_residue_sfe(arg->sigmaOutFile, arg, pdb, molSFE);

		/* chain SFE */
		if (arg->chainOut)
			print_chain_sfe(arg->sigmaOutFile, arg, pdb, molSFE);

		/* molecule SFE */
		if (! arg->noTotalOut)
			print_mol_sfe(arg->sigmaOutFile, arg, molSFE);

		fclose(arg->sigmaOutFile);

	/* for trajectory */
	} else {
		if (! arg->silent)
		fprintf(arg->sigmatrajOutFile, "%d\t%10.2f\t%10.2f\n",
				(frame + 1), molSFE->sfe_type, molSFE->sfe_group);
	}
}

