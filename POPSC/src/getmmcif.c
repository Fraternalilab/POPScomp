/*==============================================================================
getmmcif.c : routines for reading MMCIF structures
Copyright (C) 2004 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "cif_reader.h"
#include "getmmcif.h"

/*____________________________________________________________________________*/
/* map MMCIF structure */
int map_structure_mmcif(Arg *arg, Argpdb *argpdb, Str *str, Structure *s) {

	int i;
	int min_res = INT_MAX;
	int max_res = INT_MIN;

	/*____________________________________________________________________________*/
	/* initialise/allocate memory for set of (64) selected (CA) atom entries */
	str->nAtom = 0;
	str->nAllAtom = 0;
	str->nResidue = 0;
	str->nAllResidue = 0;
	str->nChain = 0;

	/*____________________________________________________________________________*/
	/* number of atoms */
	str->nAtom = s->natom;
	str->nAllAtom = str->nAtom;

	/* number of residues */
	min_res = INT_MAX;
	max_res = INT_MIN;
	for (i = 0; i < str->nAtom; ++ i) {
		min_res = s->res_number[i] < min_res ? s->res_number[i] : min_res;
		max_res = s->res_number[i] > max_res ? s->res_number[i] : max_res;
	}
	str->nResidue = max_res - min_res + 1;
	str->nAllResidue = str->nResidue;

	/* number of chains */
	str->nChain = s->chain_number;

	/*____________________________________________________________________________*/
	/* allocate memory for structure */
	str->atom = safe_malloc(str->nAtom * sizeof(Atom));
	str->atomMap = safe_malloc(str->nAtom * sizeof(int));

	/* array of residue-centric atom indices */
	str->resAtom = safe_malloc(str->nResidue * sizeof(int));

	/* allocate memory for sequence residues */
	str->sequence.res = safe_malloc(str->nResidue * sizeof(char));

	/*____________________________________________________________________________*/
	/* map entries from MMCIF structure 's' to PDB structure 'str' */
	/* that could be done directly, but I wanted to have a clean separation
		between the C++ reader and the 'getpdb'-type assignment that is used
		in the PDB and XML reading functions */

	for (i = 0; i < str->nAtom; ++ i) {
		/* atoms */
		str->atom[i].atomNumber = s->atom_number[i];
		strcpy(str->atom[i].atomName, s->atom_name[i]);

		/* residues */
		str->atom[i].residueNumber = s->res_number[i];
		strcpy(str->atom[i].residueName, s->res_name[i]);
		str->atom[i].icode[0] = s->ins_code[i];
		str->atom[i].icode[1] = '\0';

		/* chains */
		strcpy(str->atom[i].chainIdentifier, s->chain_name[i]);

		/* coordinates */
		str->atom[i].pos.x = s->xyz[3*i + 0];
		str->atom[i].pos.y = s->xyz[3*i + 1];
		str->atom[i].pos.z = s->xyz[3*i + 2];
	}


	/*____________________________________________________________________________*/
    free_structure(s);

    return 0;
}

