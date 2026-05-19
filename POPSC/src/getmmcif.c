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
	int min_chain = INT_MAX;
	int max_chain = INT_MIN;

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
	/*
	min_chain = INT_MAX;
	max_chain = INT_MIN;
	for (i = 0; i < str->nAtom; ++ i) {
		min_chain = s->chain_number[i] < min_chain ? s->chain_number[i] : min_chain;
		max_chain = s->chain_number[i] > max_chain ? s->chain_number[i] : max_chain;
	}
	str->nChain = max_chain - min_chain + 1;
	*/
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
    free_structure(s);

    return 0;
}

