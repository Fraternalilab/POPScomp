/*===============================================================================
atomgroup.h : definition of atom groups
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef ATOMGROUP_H
#define ATOMGROUP_H

#include "getpdb.h"
#include "modstring.h"

/*___________________________________________________________________________*/
/* data structures */
/* holds atom-specific data for each residue */
typedef struct {
	char residueName[4]; /* for PDB entries */
	char atomName[8];
	int atomType; /* GROMOS atom type number */
	int groupID; /* group ID number (sequential index) */
} AtomdataGroup;

/* holds information about residues */
typedef struct {
	int nResidueType; /* number of residue types */
	int nAtomResidue[64]; /* number of atoms per residue type (in the order given by 'atomgroup') */
	AtomdataGroup atomDataGroup[64][64]; /* atom-specific data */
} Atomgroup;

/*____________________________________________________________________________*/
/* prototypes */
int get_atomgroup(Str *str, Atomgroup *atomgroup);

#endif

