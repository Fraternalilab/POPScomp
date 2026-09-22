/*==============================================================================
getmmcif.c : routines for reading MMCIF structures
Copyright (C) 2004 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "cif_reader.h"
#include "getmmcif.h"

/* aacode() is shared with the PDB reader (getpdb.c) */

/*____________________________________________________________________________*/
/** standardise non-standard atom names */
__inline__ static int standardise_name(char *residueName, char *atomName)
{
	/* GRO 'ILE CD' to PDB 'ILE CD1' */
	/* atom names are unpadded in mmCIF */
	if ((strcmp(residueName, "ILE") == 0) && (strcmp(atomName, "CD") == 0))
		strcpy(atomName, "CD1");

	return 0;
}

/*____________________________________________________________________________*/
/** process HET residues */
__inline__ static int process_het(Str *str, char *line, regex_t *regexPattern, char (*hetAtomNewname)[32], int nHetAtom)
{
	int hetAtomNr = -1;

	/* atom name: assign only allowed atom elements, otherwise atom is skipped */
	if ((hetAtomNr = match_patterns(regexPattern, nHetAtom, &(str->atom[str->nAtom].atomName[0]))) >= 0) {
		/* store original atom name in 'Het' and overwrite with new name
			that is a generic name for the SASA parameters */
		sprintf(str->atom[str->nAtom].atomNameHet, "%s", &(str->atom[str->nAtom].atomName[0]));
		sprintf(str->atom[str->nAtom].atomName, "%s", &(hetAtomNewname[hetAtomNr][0]));

		sprintf(str->atom[str->nAtom].residueNameHet, "%s", &(str->atom[str->nAtom].residueName[0]));
		sprintf(str->atom[str->nAtom].residueName, "%s", "HET");
		/* set heteroatom flag */
		str->atom[str->nAtom].het = 1;
	} else {
		WarningSpec("Skipping HETATM", str->atom[str->nAtom].atomName);
		return 1;
	}

	return 0;
}

/*____________________________________________________________________________*/
/* map MMCIF structure */
int map_structure_mmcif(Arg *arg, Argpdb *argpdb, Str *str, Structure *s) {

	int i = 0; /* counter for PDB entries (skipping some MMCIF entries) */
	int x = -1; /* counter for MMCIF entries */
	unsigned int k = 0;
	char residueAltloc = 0; /* altloc selected for the current residue */

	/*____________________________________________________________________________*/
	/* initialise/allocate memory for set of selected atom entries */
	str->nAtom = 0;
	str->nAllAtom = s->natom; /* all atom records of the model: trajectory frame size */
	str->nResidue = 0;
	str->nAllResidue = 0;
	str->nChain = 0;

	/*____________________________________________________________________________*/
	/* allocate memory for structure; at most one entry per mmCIF atom */
	str->atom = safe_malloc((s->natom + 1) * sizeof(Atom));
	str->atomMap = safe_malloc((s->natom + 1) * sizeof(int));

	/* array of residue-centric atom indices: at most one per atom */
	str->resAtom = safe_malloc((s->natom + 1) * sizeof(int));

	/* allocate memory for sequence residues */
	str->sequence.res = safe_malloc((s->natom + 1) * sizeof(char));

	/*____________________________________________________________________________*/
	/* Map entries from MMCIF structure 's' to PDB structure 'str'. */
	/* That could be mapped directly in the MMCIF reader,
		but for consistency and potential debugging purpose,
		I wanted to have a clean separation between
		the C++ reader and the 'getpdb'-type assignment
		that is also used in the PDB and XML reading functions. */

	for (x = 0; x < s->natom; ++ x) {

		/* select one alternative location per residue: the first one encountered */
		if (x == 0 || s->res_number[x] != s->res_number[x - 1] ||
			s->ins_code[x] != s->ins_code[x - 1] ||
			strcmp(s->chain_name[x], s->chain_name[x - 1]) != 0) {
			residueAltloc = 0;
		}
		if (s->altloc[x] != ' ') {
			if (residueAltloc == 0)
				residueAltloc = s->altloc[x];
			if (s->altloc[x] != residueAltloc)
				continue;
		}

		/* skip hydrogens/deuteriums */
		if (! argpdb->hydrogens && s->element[x] != NULL &&
			(strcmp(s->element[x], "H") == 0 || strcmp(s->element[x], "D") == 0)) {
			continue;
		}

		/* only polymer (ATOM) records are processed, HETATM entries are skipped */
		if (s->record_type[x] != 'A') {
			continue;
		}

		/* water is not part of the solute surface */
		if (is_water(s->res_name[x])) {
			continue;
		}

		memset(&(str->atom[i]), 0, sizeof(Atom));
		str->atom[i].recordIndex = x;
		strcpy(str->atom[i].recordName, "ATOM  ");

		/* atoms */
		str->atom[i].atomNumber = s->atom_number[x];
		snprintf(str->atom[i].atomName, sizeof(str->atom[i].atomName), "%s", s->atom_name[x]);
		str->atom[i].alternativeLocation[0] = s->altloc[x];
		str->atom[i].alternativeLocation[1] = '\0';

		/* residues */
		str->atom[i].residueNumber = s->res_number[x];
		snprintf(str->atom[i].residueName, sizeof(str->atom[i].residueName), "%s", s->res_name[x]);
		/* no insertion code is written as '-', as in the PDB reader */
		str->atom[i].icode[0] = (s->ins_code[x] == ' ' || s->ins_code[x] == '\0') ? '-' : s->ins_code[x];
		str->atom[i].icode[1] = '\0';

		/* chains */
		if (strlen(s->chain_name[x]) >= sizeof(str->atom[i].chainIdentifier))
			WarningSpec("Truncating chain identifier", s->chain_name[x]);
		snprintf(str->atom[i].chainIdentifier, sizeof(str->atom[i].chainIdentifier), "%s", s->chain_name[x]);

		/* element */
		snprintf(str->atom[i].element, sizeof(str->atom[i].element), "%s", s->element[x]);

		/* coordinates */
		str->atom[i].pos.x = s->xyz[3*x + 0];
		str->atom[i].pos.y = s->xyz[3*x + 1];
		str->atom[i].pos.z = s->xyz[3*x + 2];

		/* detect CA/N3 atoms */
		if ((strcmp(str->atom[i].atomName, "CA") == 0) ||
			(strcmp(str->atom[i].atomName, "N3") == 0)) {
			str->resAtom[k] = i;
			str->sequence.res[k++] = aacode(str->atom[i].residueName);
		}

		standardise_name(str->atom[i].residueName, str->atom[i].atomName);

		/* in coarse mode record only the representative atoms: CA (amino acids), P (nucleotides) */
		if (argpdb->coarse && ! is_coarse_atom(str->atom[i].residueName, str->atom[i].atomName)) {
			continue;
		}

		/* map copied atom index i to original mmCIF atom index */
		str->atomMap[i] = x;

		++i;
	}

	str->nAtom = i;
	str->nResidue = k;
	str->sequence.res[k] = '\0';

	/* residue and chain indices, residue and chain counts */
	index_structure(str);
	/* identifier for output file names */
	set_pdbID_from_filename(str, arg->mmcifInFileName);

	/*____________________________________________________________________________*/
    free_structure(s);

    return 0;
}
