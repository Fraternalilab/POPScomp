/*==============================================================================
putsasa.c : routines for printing SASA
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "putsasa.h"

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
			pdb->pdbID, pdb->nChain, pdb->nResidue, pdb->nAllResidue, pdb->nAtom);
	else
		fprintf(sasaOutFile,\
			"\nProtein %8s\n%8d chains\n%8d standard residues\n"
			"%8d total residues (standard residues + HETATM residues lacking CA or P)\n"
			"%8d atoms (C-alpha and P atoms)\n",
			pdb->pdbID, pdb->nChain, pdb->nResidue, pdb->nAllResidue, pdb->nAtom);
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
/** print a surface ratio; 'NA' when there is no reference surface */
static void print_ratio(FILE *outFile, float sasa, float surface)
{
	if (surface > 0.)
		fprintf(outFile, "%10.4f", sasa / surface);
	else
		fprintf(outFile, "%10s", "NA");
}

/*___________________________________________________________________________*/
/** print atom SASA  */
static void print_atom_sasa(FILE *sasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
	unsigned int i, j;

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
			fprintf(sasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.4f\t%8d\t%2d\t%2d\t%10.2f\n",
				j, "XXX", pdb->atom[0].residueName, "-", pdb->atom[0].residueNumber, "-",
				0., 0., 0, 0, 0, 0.);
		}
	}

	for (i = 0; i < pdb->nAtom; ++ i) {
		fprintf(sasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t",
			pdb->atom[i].atomNumber,
			pdb->atom[i].atomName,
			pdb->atom[i].residueName,
			chain_label(&(pdb->atom[i])),
			pdb->atom[i].residueNumber,
			pdb->atom[i].icode,
			molSasa->atomSasa[i].sasa);
		print_ratio(sasaOutFile, molSasa->atomSasa[i].sasa, molSasa->atomSasa[i].surface);
		fprintf(sasaOutFile, "\t%8d\t%2d\t%2d\t%10.2f\n",
			molSasa->atomSasa[i].nOverlap,
			pdb->atom[i].atomType,
			pdb->atom[i].groupID,
			molSasa->atomSasa[i].surface);
		/* after the first atom :
			print dummy lines to substitute missing hydrogen atom lines */
		if (arg->padding && ((i + 1) < pdb->nAtom)) {
			j = pdb->atom[i].atomNumber;
			while (++ j < pdb->atom[i+1].atomNumber) {
				fprintf(sasaOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.2f\t%10.4f\t%8d\t%2d\t%2d\t%10.2f\n",
					j, "HXX", pdb->atom[i].residueName, "-", pdb->atom[i].residueNumber, "-",
					0., 0., 0, 0, 0, 0.);
			}
		}
	}
}

/*___________________________________________________________________________*/
/** print residue SASA */
static void print_residue_sasa(FILE *sasaOutFile, Arg *arg, Str *pdb, MolSasa *molSasa)
{
    unsigned int i;
	Atom *ref;

	if (! arg->noHeaderOut && ! arg->rout) {
		fprintf(sasaOutFile, "\n=== RESIDUE SASAs ===\n");
	}
	if (! arg->noHeaderOut || arg->rout) {
		fprintf(sasaOutFile, "ResidNe\tChain\tResidNr\tiCode\tPhob/A^2\tPhil/A^2\tSASA/A^2\tQ(SASA)\tN(overl)\tSurf/A^2\n");
	}

    for (i = 0; i < pdb->nAllResidue; ++ i) { 
		ref = &(pdb->atom[molSasa->resSasa[i].atomRef]);
		fprintf(sasaOutFile, "%8s\t%3s\t%8d\t%1s\t%10.2f\t%10.2f\t%10.2f\t",
			ref->residueName,
			chain_label(ref),
			ref->residueNumber,
			ref->icode,
			molSasa->resSasa[i].phobicSasa,
			molSasa->resSasa[i].philicSasa,
			molSasa->resSasa[i].sasa);
		/* upper limit: there can be configurational effects leading to >1. */
		print_ratio(sasaOutFile, molSasa->resSasa[i].sasa, molSasa->resSasa[i].surface);
		fprintf(sasaOutFile, "\t%8d\t%10.2f\n",
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
		fprintf(sasaOutFile, "Chain\tId\tAtomRange\tResidRange\tPhob/A^2\tPhil/A^2\tSASA/A^2\n");
	}

    for (i = 0; i < pdb->nChain; ++ i) {
		fprintf(sasaOutFile, "%3d\t%3s\t%6d->%-6d\t%5d->%-5d\t%10.2f\t%10.2f\t%10.2f\n",
			i,
			chain_label(&(pdb->atom[molSasa->chainSasa[i].first])),
			pdb->atom[molSasa->chainSasa[i].first].atomNumber,
			pdb->atom[molSasa->chainSasa[i].last].atomNumber,
			pdb->atom[molSasa->chainSasa[i].first].residueNumber,
			pdb->atom[molSasa->chainSasa[i].last].residueNumber,
			molSasa->chainSasa[i].phobicSasa,
			molSasa->chainSasa[i].philicSasa,
			molSasa->chainSasa[i].sasa);
	}
}

/*____________________________________________________________________________*/
/** total (molecule) SASA */
void print_mol_sasa(FILE *sasaOutFile, Arg *arg, MolSasa *molSasa)
{
	if (! arg->noHeaderOut && ! arg->rout) {
		fprintf(sasaOutFile, "\n=== MOLECULE SASAs ===\n");
	}
	if (! arg->noHeaderOut || arg->rout) {
		fprintf(sasaOutFile, "Phob/A^2\tPhil/A^2\tSASA/A^2\n");
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
		fprintf(neighbourOutFile, "%d:%s\t%d\t", 
			pdb->atom[i].atomNumber, chain_label(&(pdb->atom[i])), topol->neighbourState[i][0]);
		for (j = 1; j <= topol->neighbourState[i][0]; ++ j)
			fprintf(neighbourOutFile, "%d:%s ",
				pdb->atom[topol->neighbourState[i][j]].atomNumber,
				chain_label(&(pdb->atom[topol->neighbourState[i][j]])));
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
/** print interface (residue pairs) */
/* nearest neighbour atom pairs on different chains */
void print_interface(FILE *interfaceOutFile, Arg *arg, Str *pdb, Type *type, Topol *topol)
{
	unsigned int i;
	int j;

	for (i = 0; i < pdb->nAtom; ++ i) {
		/* nearest overlapping atom on a different chain; -1 if there is none */
		if ((j = topol->interfaceNn[i]) < 0)
			continue;

		fprintf(interfaceOutFile, "%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%8d\t%3s\t%3s\t%1s\t%6d\t%1s\t%10.4f\n",
			pdb->atom[i].atomNumber,
			pdb->atom[i].atomName,
			pdb->atom[i].residueName,
			chain_label(&(pdb->atom[i])),
			pdb->atom[i].residueNumber,
			pdb->atom[i].icode,
			pdb->atom[j].atomNumber,
			pdb->atom[j].atomName,
			pdb->atom[j].residueName,
			chain_label(&(pdb->atom[j])),
			pdb->atom[j].residueNumber,
			pdb->atom[j].icode,
			topol->interfaceNnDist[i]);
	}
}

/*____________________________________________________________________________*/
/** open an R-format (POPScomp) output file '<prefix>[.<frame>].<suffix>' */
static FILE *open_rpops(Arg *arg, int frame, const char *suffix, const char *mode)
{
	char name[1024];
	int n;

	if (frame < 0)
		n = snprintf(name, sizeof(name), "%s.%s", arg->routPrefix, suffix);
	else
		n = snprintf(name, sizeof(name), "%s.%d.%s", arg->routPrefix, frame, suffix);
	if (n < 0 || (size_t)n >= sizeof(name))
		ErrorSpec("Output file name too long", arg->routPrefix);

	return open_output(arg, name, mode);
}

/*____________________________________________________________________________*/
/** open the SASA output file of the reference molecule (frame < 0) or of a frame */
static FILE *open_sasa(Arg *arg, int frame)
{
	char name[1024];
	int n;

	if (frame < 0)
		return open_output(arg, arg->sasaOutFileName, "w");

	n = snprintf(name, sizeof(name), "%s.%d.%s", arg->sasatrajOutFileName, frame, "out");
	if (n < 0 || (size_t)n >= sizeof(name))
		ErrorSpec("Output file name too long", arg->sasatrajOutFileName);
	return open_output(arg, name, "w");
}

/*____________________________________________________________________________*/
/** print SASAs */
void print_sasa(Arg *arg, Argpdb *argpdb, Str *pdb, Type *type, Topol *topol, \
				MolSasa *molSasa, ConstantSasa *constant_sasa, int frame)
{
	FILE *rpopsOutFile;
	/* files that collect all frames of a trajectory are appended to */
	const char *collectMode = (frame < 0) ? "w" : "a";

	/* SASA output file: one per molecule or frame; R format writes one file per table */
	if (frame < 0 && ! arg->silent)
		fprintf(stdout, "\tSASA of reference molecule: %s\n", arg->sasaOutFileName);
	arg->sasaOutFile = arg->rout ? NULL : open_sasa(arg, frame);

	/* composition */
	if (arg->compositionOut) {
		rpopsOutFile = arg->rout ? open_rpops(arg, frame, "rpopsComposition", "w") : arg->sasaOutFile;
		print_composition(rpopsOutFile, arg, argpdb, pdb);
		if (arg->rout) fclose(rpopsOutFile);
	}

	/* topology */
	if (arg->topologyOut) {
		rpopsOutFile = arg->rout ? open_rpops(arg, frame, "rpopsTopology", "w") : arg->sasaOutFile;
		print_topology(rpopsOutFile, arg, topol);
		if (arg->rout) fclose(rpopsOutFile);
	}

	/* residue/atom types */
	if (arg->typeOut) {
		rpopsOutFile = arg->rout ? open_rpops(arg, frame, "rpopsTypes", "w") : arg->sasaOutFile;
		print_types(rpopsOutFile, arg, type, constant_sasa);
		if (arg->rout) fclose(rpopsOutFile);
	}

	/* atom SASA */
	if (arg->atomOut && ! argpdb->coarse) {
		rpopsOutFile = arg->rout ? open_rpops(arg, frame, "rpopsAtom", "w") : arg->sasaOutFile;
		print_atom_sasa(rpopsOutFile, arg, pdb, molSasa);
		if (arg->rout) fclose(rpopsOutFile);
	}
	
	/* print '0' to pops.out.rpopsAtom for Shiny reactive file reader
	     when no atom SASAs are being computed under '--coarse' */
	if (argpdb->coarse && arg->rout) {
		rpopsOutFile = open_rpops(arg, frame, "rpopsAtom", "w");
		fprintf(rpopsOutFile, "%d\n", 0);
		fclose(rpopsOutFile);
	}

	/* residue SASA */
	if (arg->residueOut) {
		rpopsOutFile = arg->rout ? open_rpops(arg, frame, "rpopsResidue", "w") : arg->sasaOutFile;
		print_residue_sasa(rpopsOutFile, arg, pdb, molSasa);
		if (arg->rout) fclose(rpopsOutFile);
	}

	/* chain SASA */
	if (arg->chainOut) {
		rpopsOutFile = arg->rout ? open_rpops(arg, frame, "rpopsChain", "w") : arg->sasaOutFile;
		print_chain_sasa(rpopsOutFile, arg, pdb, molSasa);
		if (arg->rout) fclose(rpopsOutFile);
	}

	/* molecule SASA */
	if (! arg->noTotalOut) {
		rpopsOutFile = arg->rout ? open_rpops(arg, frame, "rpopsMolecule", "w") : arg->sasaOutFile;
		print_mol_sasa(rpopsOutFile, arg, molSasa);
		if (arg->rout) fclose(rpopsOutFile);
	}

	if (! arg->rout) {
		fclose(arg->sasaOutFile);
		arg->sasaOutFile = NULL;
	}

	/* neighbour list */
	if (arg->neighbourOut) {
		arg->neighbourOutFile = arg->rout ? open_rpops(arg, -1, "neighbours.out", collectMode) :
								open_output(arg, arg->neighbourOutFileName, collectMode);
		print_neighbour_list(arg->neighbourOutFile, arg, pdb, topol);
		fclose(arg->neighbourOutFile);
	}

	/* neighbour parameters (for benchmarking) */
	if (arg->parameterOut) {
		arg->parameterOutFile = open_output(arg, arg->parameterOutFileName, collectMode);
		print_neighbour_parameter(arg->parameterOutFile, pdb, type, topol, molSasa);
		fclose(arg->parameterOutFile);
	}

	/* interface residue pairs */
	if (arg->interfaceOut) {
		arg->interfaceOutFile = arg->rout ? open_rpops(arg, -1, "interface.out", collectMode) :
								open_output(arg, arg->interfaceOutFileName, collectMode);
		print_interface(arg->interfaceOutFile, arg, pdb, type, topol);
		fclose(arg->interfaceOutFile);
	}
}
