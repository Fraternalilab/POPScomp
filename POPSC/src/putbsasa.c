/*==============================================================================
putbsasa.c : routines for printing bSASA
Copyright (C) 2016-2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "putbsasa.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*___________________________________________________________________________*/
/** print atom bSASA  */
static void print_atom_bsasa(FILE *bsasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
	unsigned int i, j;

	if (! arg->noHeaderOut) {
		fprintf(bsasaOutFile, "\n=== ATOM bSASAs ===\n");
		fprintf(bsasaOutFile, "\nAtomNr\tAtomNe\tResidNe\tChain\tResidNr\tiCode\tcPhob/A^2\tQ(cPhob)\tcPhil/A^2\tQ(cPhil)\tcTotal/A^2\n");
	}

	/* before the first line:
		print dummy lines to substitute missing hydrogen atom lines */
	if (arg->padding)  {
		j = 0;
		while (++ j < pdb->atom[0].atomNumber) {
			fprintf(bsasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.8f\t%10.2f\t%10.8f\t%10.2f\n",
				j,
				"XXX",
				pdb->atom[0].residueName,
				" ",
				pdb->atom[0].residueNumber,
				" ",
				0.,
				0.,
				0.,
				0.,
				0.);
		}
	}

	for (i = 0; i < pdb->nAtom; ++ i) {
		fprintf(bsasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.8f\t%10.2f\t%10.8f\t%10.2f\n",
			pdb->atom[i].atomNumber,
			pdb->atom[i].atomName,
			pdb->atom[i].residueName,
			pdb->atom[i].chainIdentifier,
			pdb->atom[i].residueNumber,
			pdb->atom[i].icode,
			molSasa->atomSasa[i].phobicbSasa,
			molSasa->atomSasa[i].bSasa > 0. ? (molSasa->atomSasa[i].phobicbSasa / molSasa->atomSasa[i].bSasa) : 0.,
			molSasa->atomSasa[i].philicbSasa,
			molSasa->atomSasa[i].bSasa > 0. ? (molSasa->atomSasa[i].philicbSasa / molSasa->atomSasa[i].bSasa) : 0.,
			molSasa->atomSasa[i].bSasa);
		/* after the first atom :
			print dummy lines to substitute missing hydrogen atom lines */
		if (arg->padding && ((i + 1) < pdb->nAtom)) {
			j = pdb->atom[i].atomNumber;
			while (++ j < pdb->atom[i+1].atomNumber) {
				fprintf(bsasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.8f\t%10.2f\t%10.8f\t%10.2f\n",
					j,
					"HXX",
					pdb->atom[i].residueName,
					" ",
					pdb->atom[i].residueNumber,
					" ",
					0.,
					0.,
					0.,
					0.,
					0.);
			}
		}
	}
}

/*___________________________________________________________________________*/
/** print residue bSASA */
static void print_residue_bsasa(FILE *bsasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
    unsigned int i;

	if (! arg->noHeaderOut) {
		fprintf(bsasaOutFile, "\n=== RESIDUE bSASAs ===\n");
		fprintf(bsasaOutFile, "\nResidNe\tChain\tResidNr\tiCode\tcPhob/A^2\tQ(cPhob)\tcPhil/A^2\tQ(cPhil)\tcTotal/A^2\n");
	}

    for (i = 0; i < pdb->nAllResidue; ++ i) { 
		fprintf(bsasaOutFile, "%8s\t%3s\t%8d\t%1s\t%10.2f\t%10.8f\t%10.2f\t%10.8f\t%10.2f\n",
			pdb->atom[molSasa->resSasa[i].atomRef].residueName,
			pdb->atom[molSasa->resSasa[i].atomRef].chainIdentifier,
			pdb->atom[molSasa->resSasa[i].atomRef].residueNumber,
			pdb->atom[molSasa->resSasa[i].atomRef].icode,
			molSasa->resSasa[i].phobicbSasa,
			molSasa->resSasa[i].bSasa > 0. ? (molSasa->resSasa[i].phobicbSasa / molSasa->resSasa[i].bSasa) : 0.,
			molSasa->resSasa[i].philicbSasa,
			molSasa->resSasa[i].bSasa > 0. ? (molSasa->resSasa[i].philicbSasa / molSasa->resSasa[i].bSasa) : 0.,
			molSasa->resSasa[i].bSasa);
    }
}

/*___________________________________________________________________________*/
/** print chain bSASA */
static void print_chain_bsasa(FILE *bsasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
    unsigned int i;

	if (! arg->noHeaderOut) {
		fprintf(bsasaOutFile, "\n=== CHAIN bSASAs ===\n(Atom Range excluding hydrogen atoms)\n");
		fprintf(bsasaOutFile, "\nChain\tId\tAtomRange\tResidRange\tcPhob/A^2\tcPhil/A^2\t\tcTotal/A^2\n");
	}

    for (i = 0; i < pdb->nChain; ++ i)
		fprintf(bsasaOutFile, "%3d\t%3s\t%6d->%-6d\t%5d->%-5d\t%10.2f\t%10.2f\%10.2f\n",
			i,
			pdb->atom[molSasa->chainSasa[i].first].chainIdentifier,
			pdb->atom[molSasa->chainSasa[i].first].atomNumber,
			pdb->atom[molSasa->chainSasa[i].last].atomNumber,
			pdb->atom[molSasa->chainSasa[i].first].residueNumber,
			pdb->atom[molSasa->chainSasa[i].last].residueNumber,
			molSasa->chainSasa[i].phobicbSasa,
			molSasa->chainSasa[i].philicbSasa,
			molSasa->chainSasa[i].bSasa);
}

/*____________________________________________________________________________*/
/** total (molecule) bSASA */
void print_mol_bsasa(FILE *bsasaOutFile, Arg *arg, MolSasa *molSasa)
{
	if (! arg->noHeaderOut) fprintf(bsasaOutFile, "\n=== MOLECULE bSASAs ===\n");
		fprintf(bsasaOutFile, "\ncPhob/A^2\tcPhil/A^2\t\tcTotal/A^2\n");
    fprintf(bsasaOutFile, "%10.2f\t%10.2f\t%10.2f\n\n",
			molSasa->phobicbSasa,
			molSasa->philicbSasa,
			molSasa->bSasa);
}

/*____________________________________________________________________________*/
/** print bSASAs */
void print_bsasa(Arg *arg, Argpdb *argpdb, Str *pdb, Type *type, Topol *topol, \
				MolSasa *molSasa, ConstantSasa *constant_sasa, int frame)
{
	char bsasatrajOutFileName[256];

	/* for single (reference) molecule */
	if (frame < 0) {
		if (! arg->silent) {
			fprintf(stdout, "\tbSASA of reference molecule: %s\n", arg->bsasaOutFileName);
			arg->bsasaOutFile = safe_open(arg->bsasaOutFileName, "w");
		}
	} else {
			sprintf(&(bsasatrajOutFileName[0]), "%s.%d.%s", arg->bsasatrajOutFileName, frame, "out");
			arg->bsasaOutFile = safe_open(bsasatrajOutFileName, "w");
	}

	/* atom bSASA */
	if (arg->atomOut && ! argpdb->coarse)
		print_atom_bsasa(arg->bsasaOutFile, arg, pdb, molSasa);

	/* residue bSASA */
	if (arg->residueOut)
		print_residue_bsasa(arg->bsasaOutFile, arg, pdb, molSasa);

	/* chain bSASA */
	if (arg->chainOut)
		print_chain_bsasa(arg->bsasaOutFile, arg, pdb, molSasa);

	/* molecule bSASA */
	if (! arg->noTotalOut)
		print_mol_bsasa(arg->bsasaOutFile, arg, molSasa);

	fclose(arg->bsasaOutFile);
}

