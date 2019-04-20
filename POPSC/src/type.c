/*==============================================================================
type.c : determine atom and residue types
Copyright (C) 2002-2018 Franca Fraternali
Copyright (C) 2008-2018 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin and Valerie Hindie
Read the COPYING file for license information.
==============================================================================*/

#include "type.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/** get atom and residue types */
int get_types(Str *pdb, Type *type, ConstantSasa *constant_sasa)
{
	unsigned int i, j;
	char p_residue[8] = "";
	char p_atom[4] = "";
	int found_residue = 0;
	int found_atom = 0;

	type->atomType = safe_malloc(pdb->nAtom * sizeof(int));
	type->residueType = safe_malloc(pdb->nAtom * sizeof(int));

	/*____________________________________________________________________________*/
	for (i = 0; i < pdb->nAtom; ++ i) {
		found_residue = 0;
		found_atom = 0;

		/*____________________________________________________________________________*/
		/* residue */
		strip_char(pdb->atom[i].residueName, p_residue); /* remove spaces */
		/* ... match its residue type to the ones in constant SASA */
		for (j = 0, found_residue = 0; j < constant_sasa->nResidueType; ++ j) {
			/* if matched, assign residue type 'j' to atom 'i' */
			if (strcmp(p_residue, constant_sasa->atomDataSasa[j][0].residueName) == 0) {
				type->residueType[i] = j;
				++ found_residue; /* flag up match */
				break;
			}
		}

		/* if no match is found for this residue name */
		if (! found_residue) {
			strip_char(pdb->atom[i].residueName, p_residue); /* remove spaces */
			fprintf(stderr, "Warning: Unknown residue name '%s' -> '%s' of residue number %d\n"
							"Setting residue type to UNK (for unknown polymer residues)\n",
				pdb->atom[i].residueName, p_residue, pdb->atom[i].residueNumber);
			type->residueType[i] = 37; /* UNK is residue number 37 in 'sasa_data.h' */
			++ found_residue; /* flag up match */
			break;
		}

		/*____________________________________________________________________________*/
		/* atom */
		/* ... match its atom type to the ones in constant SASA */
		strip_char(pdb->atom[i].atomName, p_atom); /* remove spaces */
		for (j = 0, found_atom = 0; j < constant_sasa->nAtomResidue[type->residueType[i]]; ++ j) {
			/* if matched, assign atom type 'j' to atom 'i' */
			if (strcmp(p_atom, constant_sasa->atomDataSasa[type->residueType[i]][j].atomName) == 0) {
				type->atomType[i] = j;
				++ found_atom; /* flag up match */
				break;
			}
		}

		/* OXT is special because it is not residue-specific */
		if (! found_atom) {
			/* 'OXT' has variable residue names, but is assigned to 
				residue 'ANY' = 20 and atom 'OXT' = 0 */
			if (strcmp(p_atom, constant_sasa->atomDataSasa[20][0].atomName) == 0) {
				type->residueType[i] = 20;
				type->atomType[i] = 0;
				++ found_atom; /* flag up match */
				continue;
			}
		}

		/* O1 is special because it is not residue-specific;
			not PDB format, but frequently used in MD and graphics programs */
		if (! found_atom) {
			/* 'O1' has variable residue names, but is assigned to 
				residue 'ANY' = 20 and atom 'O1' = 1 */
			if (strcmp(p_atom, constant_sasa->atomDataSasa[20][1].atomName) == 0) {
				type->residueType[i] = 20;
				type->atomType[i] = 1;
				++ found_atom; /* flag up match */
				continue;
			}
		}

		/* O2 is special because it is not residue-specific;
			not PDB format, but frequently used in MD and graphics programs */
		if (! found_atom) {
			/* 'O2' has variable residue names, but is assigned to 
				residue 'ANY' = 20 and atom 'O2' = 2 */
			if (strcmp(p_atom, constant_sasa->atomDataSasa[20][2].atomName) == 0) {
				type->residueType[i] = 20;
				type->atomType[i] = 2;
				++ found_atom; /* flag up match */
				continue;
			}
		}


		/*____________________________________________________________________________*/
		/* unknown atom means unknown parametrisation: exit */
		if (! found_atom)
			ErrorSpec("Unknown type of atom", pdb->atom[i].atomName);
	}

	return 0;
}

