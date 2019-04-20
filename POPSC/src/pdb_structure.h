/*===============================================================================
pdb_structure.h : PDB structure definition
Copyright (C) 2004-2008 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef PDB_STRUCTURE_H
#define PDB_STRUCTURE_H

#include "seq.h"
#include "vector.h"

/*____________________________________________________________________________*/
/* structures */

/* atom : definition of PDB atom format, numbers indicate columns */
typedef struct atom
{
	struct atom *next;
	char recordName[8]; /* Record type; 1 - 6*/
	int  entityId; /* Entity ID; in XML files */
	int atomNumber; /* Atom serial number;  7 - 11 */
	char atomName[8]; /* Atom name; 13 - 16 */
	char atomName_raw[8]; /* Atom name; 13 - 14 */
	char alternativeLocation[2]; /* Alternate location indicator; 17 */
	char residueName[4]; /* Residue name; 18 - 20 */
	char chainIdentifier[2]; /* Chain identifier; 22 */
	int residueNumber; /* Residue sequence number; 23 - 26 */
	char icode[2]; /* Code for insertion of residues; 27 */
	Vec pos; /* position vector (x, y, z) */
	float occupancy; /* Occupancy; 55 - 60 */
	float temperatureFactor; /* Temperature factor; 61 - 66 */
	char segmentIdentifier[5]; /* Segment identifier; 73 - 76 */
	char element[3]; /* Element symbol; 77 - 78 */
	char charge[3]; /* Charge on the atom; 79 - 80 */
	int formalCharge;
	float partialCharge;
	char altpos; /* Alternate position indicator */
	char secstr; /* Secondary structure */
	char description[32]; /* everything before coordinates */
	Vec tpos; /* transformed position vector */
	int atomType; /* GROMOS atom type */
	int groupID; /* atom group ID */
	int modelNumber; /* model number */
	int het; /* heteroatom flag */
	char atomNameHet[8]; /* original atom name of heteroatom; 13 - 16 */
	char residueNameHet[8]; /* original residue name of heteroresidue */
} Atom;

/* residue */
typedef struct residue 
{ 
	struct residue *next, *prev; 
	Atom *start, *stop; 
	int residuenumber; 
	char chainIdentifier[2]; 
	char insert[8]; 
	char residueName[4]; 
	char residueID[8]; 
	int het; /* heterotamp flag */
	char residueNameHet[4]; /* original residue name of heteroresidue */
} Residue; 

/* chain */
typedef struct chain 
{ 
	struct chain *next, *prev; 
	Atom *start, *stop; 
	Residue *residues; 
	char chainIdentifier[2]; 
} Chain; 

/* molecular structure */
typedef struct str
{
	struct str *next, *prev;
	Atom *atom; /* array of atoms constituting structure */
	int *resAtom; /* atom indices of CA and P atoms */
	int *atomMap; /* map of the selected atom count to the original atom count */
	int nAtom; /* number of selected (CA) atoms */
	int nAllAtom; /* number of all atoms */
	int nResidue; /* number of residues (CA and P atoms) */
	int nAllResidue; /* number of all residues (including HETATM) */
	int nChain; /* number of chains */
	int multiModel; /* multiple models */
	int modelNumber; /* number of PDB model */
	Seq sequence; /* amino acid sequence of structure */
	Seq strSequence; /* sequence of string-encoded structure */
	Chain *chains;
} Str;

/* secondary structure */
typedef struct secstr
{
	struct secstr *next, *prev;
	char chain1[8];
	char insert1[8];
	char chain2[8];
	char insert2[8];
	int  res1;
	int  res2;
	char type;
}  Secstr;

/* ensemble of structures; multiple models */
typedef struct
{
	Str *str;
	int nStr;
} Ensmbl;

#endif

