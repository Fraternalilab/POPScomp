/*===============================================================================
atomgroup.c : definition of atom groups
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#include "atomgroup.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/** assign atom group ID number */
/* assign atom group number to atom type *
 * the structure 'constAtomGroup' in 'atomgroup.h' contains the definition 
 * of atom types and their group IDs */
int get_atomgroup(Str *str, Atomgroup *atomGroup)
{
	unsigned int i, j;
	int residue_type; /* position in residue array */
	int atom_type; /* position in atom array */
	char p_residue[4] = "";
	char p_atom[8] = "";
	int found_residue = 0;
	int found_atom = 0;

	for (i = 0; i < str->nAtom; ++ i) {
		/* initialise atom group ID: unmatched=0 */
		str->atom[i].atomType = 0;
		str->atom[i].groupID = 0;
		found_residue = 0;
		found_atom = 0;

		/*____________________________________________________________________________*/
		/* residue */
		/* ... match residue type to the ones in constant atom group data */
		strip_char(str->atom[i].residueName, p_residue); /* remove spaces */
		for (j = 0, residue_type = -1, found_residue = 0; j < atomGroup->nResidueType; ++ j) {
			/* if matched, assign residue type 'j' to atom 'i' */
			if (strcmp(p_residue, atomGroup->atomDataGroup[j][0].residueName) == 0) {
				residue_type = j;
				++ found_residue; /* flag up match */
				break;
			}
		}

		/* if no match is found for this residue name */
		if (! found_residue) {
			strip_char(str->atom[i].residueName, p_residue); /* remove spaces */
			fprintf(stderr, "Warning: Unknown residue name '%s' -> '%s' of residue number %d\n"
							"Setting residue type to UNK (for unknown polymer residues)\n",
				str->atom[i].residueName, p_residue, str->atom[i].residueNumber);
			residue_type = 37; /* UNK is residue number 37 in 'sasa_data.h' */
			++ found_residue; /* flag up match */
			break;
		}

		/* atom */
		/* ... match the atom type to the ones belonging to 'residue_type' */
		strip_char(str->atom[i].atomName, p_atom); /* remove spaces */
		for (j = 0, atom_type = -1, found_atom = 0; j < atomGroup->nAtomResidue[residue_type]; ++ j) {
			/* if matched, assign atom type 'j' to atom 'i' */
			if (strcmp(p_atom, atomGroup->atomDataGroup[residue_type][j].atomName) == 0) {
				atom_type = j;
				++ found_atom; /* flag up match */
				break;
			}
		}

		/* OXT is special because it is not residue-specific */
		if (! found_atom) {
			/* 'OXT' has variable residue names, but is assigned to 
				residue 'ANY' = 20 and atom 'OXT' = 0 */
			if (strcmp(p_atom, atomGroup->atomDataGroup[20][0].atomName) == 0) {
				residue_type = 20;
				atom_type = 0;
				++ found_atom; /* flag up match */
				continue;
			}
		}

		/* O1 is special because it is not residue-specific;
			not PDB format, but frequently used in MD and graphics programs */
		if (! found_atom) {
			/* 'O1' has variable residue names, but is assigned to 
				residue 'ANY' = 20 and atom 'O1' = 1 */
			if (strcmp(p_atom, atomGroup->atomDataGroup[20][1].atomName) == 0) {
				residue_type = 20;
				atom_type = 1;
				++ found_atom; /* flag up match */
				continue;
			}
		}

		/* O2 is special because it is not residue-specific;
			not PDB format, but frequently used in MD and graphics programs */
		if (! found_atom) {
			/* 'O2' has variable residue names, but is assigned to 
				residue 'ANY' = 20 and atom 'O2' = 2 */
			if (strcmp(p_atom, atomGroup->atomDataGroup[20][2].atomName) == 0) {
				residue_type = 20;
				atom_type = 2;
				++ found_atom; /* flag up match */
				continue;
			}
		}

		/* assign atom type and group ID number to atom */
		if (found_atom) {
			str->atom[i].atomType = atomGroup->atomDataGroup[residue_type][atom_type].atomType;
			str->atom[i].groupID = atomGroup->atomDataGroup[residue_type][atom_type].groupID;
		}
		else {
			WarningSpec("Unknown group of atom", str->atom[i].atomName);
			str->atom[i].atomType = 0; /* default atom group ID */
			str->atom[i].groupID = 0; /* default atom group ID */
		}
	}

	return 0;
}
