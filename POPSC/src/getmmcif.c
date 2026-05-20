/*==============================================================================
getmmcif.c : routines for reading MMCIF structures
Copyright (C) 2004 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "cif_reader.h"
#include "getmmcif.h"

/*____________________________________________________________________________*/
/* match PDB residue name against constant residue name array */
__inline__ static char scan_array(char *code3, char *residue_array[], int shift)
{
	unsigned int i;
	char residue = ' ';

	for (i = 0; i < 26; ++ i)
		if (strncmp(code3, residue_array[i], 3) == 0) {
			residue = i + shift; /* shift=65 for UPPER, shift=97 for lower */
			break;
		}

	return residue;
}

/*____________________________________________________________________________*/
/** amino acid 3-letter to 1-letter code conversion */
__inline__ static char aacode(char *code3)
{
	char residue = ' '; /* 1-letter residue name */

	/* three-letter code of amino acid residues, exception HET (X) */
	char *aa3[] = {"ALA","---","CYS","ASP","GLU","PHE","GLY","HIS","ILE","---","LYS","LEU","MET","ASN","---","PRO","GLN","ARG","SER","THR","UNL","VAL","TRP","HET","TYR","UNK"};
	/* nucleotide residues */
	char *nuc[] = {"  A"," DA","  C"," DC","---","---","  G"," DG","  I"," DI","---","---","---","  N"," DN","---","---","---"," DT","  T","  U"," DU","---","---","---","---"};

	/* match against amino acid residues */
	residue = scan_array(code3, aa3, 65);

	/* match against nucleotide residues */
	if (residue == ' ')
		residue = scan_array(code3, nuc, 97);

	/* residue not found */
	if (residue == ' ') {
		Warning("Non-standard residue.");
		residue = 'X';
	}
	
	return residue;
}

/*____________________________________________________________________________*/
/** standardise non-standard atom names */
__inline__ static int standardise_name(char *residueName, char *atomName)
{
	/* GRO 'ILE CD' to PDB 'ILE CD1' */
	if ((strcmp(residueName, "ILE") == 0) && (strcmp(atomName, " CD ") == 0))
		strcpy(atomName, " CD1");

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

	char resbuf;
	int ca_p = 0;
	/* for HETATM entries */
	regex_t *regexPattern = 0; /* regular atom patterns */
	/* allowed HETATM atom types (standard N,CA,C,O) and elements (any N,C,O,P,S) */
	const int nHetAtom = 9;
	char hetAtomPattern[9][32] = {{" N  "},{" CA "},{" C  "},{" O  "},{".{1}C[[:print:]]{1,3}"},{".{1}N[[:print:]]{1,3}"},{".{1}O[[:print:]]{1,3}"},{".{1}P[[:print:]]{1,3}"},{".{1}S[[:print:]]{1,3}"}};

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

	/* number of residues */
	str->nResidue = s->nresidue; 

	/* number of chains */
	str->nChain = s->chain_number;

	/* for debigging only
	printf("nAtom %d , nResidue %d , nChain %d\n", str->nAtom, str->nResidue, str->nChain);
	*/

	/*____________________________________________________________________________*/
	/* allocate memory for structure */
	str->atom = safe_malloc(str->nAtom * sizeof(Atom));
	str->atomMap = safe_malloc(str->nAtom * sizeof(int));

	/* array of residue-centric atom indices */
	str->resAtom = safe_malloc((str->nResidue + 1) * sizeof(int));

	/* allocate memory for sequence residues */
	str->sequence.res = safe_malloc((str->nResidue + 1) * sizeof(char));

	/* compile allowed HETATM element patterns */
	regexPattern = safe_malloc(nHetAtom * sizeof(regex_t));
	compile_patterns(regexPattern, &(hetAtomPattern[0]), nHetAtom);

	/*____________________________________________________________________________*/
	/* Map entries from MMCIF structure 's' to PDB structure 'str'. */
	/* That could be mapped directly in the MMCIF reader,
		but for consistency and potential debugging purpose,
		I wanted to have a clean separation between
		the C++ reader and the 'getpdb'-type assignment
		that is also used in the PDB and XML reading functions. */

	for (x = 0; x < s->natom; ++ x) {

		ca_p = 0;

		/* Skip hydrogens/deuteriums before copying into str->atom[i] */
		if (!argpdb->hydrogens) {
			if (s->element[x][0] == 'H' || s->element[x][0] == 'D') {
				++str->nAllAtom;
				continue;
			}
		}

		if (!argpdb->hydrogens) {
		    if (s->element[x] != NULL &&
      		  (s->element[x][0] == 'H' || s->element[x][0] == 'D')) {
      		  ++str->nAllAtom;
       		 continue;
    		}
		}

		/* atoms */
		str->atom[i].atomNumber = s->atom_number[x];
		strcpy(str->atom[i].atomName, s->atom_name[x]);
		str->atom[i].alternativeLocation[0] = s->altloc[x];
		str->atom[i].alternativeLocation[1] = '\0';

		/* residues */
		str->atom[i].residueNumber = s->res_number[x];
		strcpy(str->atom[i].residueName, s->res_name[x]);
		str->atom[i].icode[0] = s->ins_code[x];
		str->atom[i].icode[1] = '\0';

		/* chains */
		strcpy(str->atom[i].chainIdentifier, s->chain_name[x]);

		/* coordinates */
		str->atom[i].pos.x = s->xyz[3*x + 0];
		str->atom[i].pos.y = s->xyz[3*x + 1];
		str->atom[i].pos.z = s->xyz[3*x + 2];

		/* record type */
		if (s->record_type[x] == 'A') {
			str->atom[i].het = 0;
		} else if (s->record_type[x] == 'H') {
			str->atom[i].het = 1;
		} else {
			++str->nAllAtom;
			continue;
		}

		/* aa code */
		if (str->atom[i].het == 0) {
			assert((resbuf = aacode(str->atom[i].residueName)) != ' ');
		}

		/* Skip HETATM entries */
		if (str->atom[i].het == 1) {
			++str->nAllAtom;
			continue;
		}

		/* detect CA/N3 atoms */
		if ((strncmp(str->atom[i].atomName, " CA ", 4) == 0) ||
			(strncmp(str->atom[i].atomName, " N3 ", 4) == 0)) {

			str->resAtom[k] = i;
			str->sequence.res[k++] = aacode(str->atom[i].residueName);
			ca_p = 1;
		}

		standardise_name(str->atom[i].residueName, str->atom[i].atomName);

		/* in coarse mode record only CA and P entries */
		if (!ca_p && argpdb->coarse) {
			++str->nAllAtom;
			continue;
		}

		/* count residues among copied heavy ATOM entries */
		if (i == 0 ||
			str->atom[i].residueNumber != str->atom[i - 1].residueNumber ||
			strcmp(str->atom[i].icode, str->atom[i - 1].icode) != 0 ||
			strcmp(str->atom[i].chainIdentifier, str->atom[i - 1].chainIdentifier) != 0) {

			++str->nAllResidue;
		}

		/* for debugging only
		printf("x %d : i %d : %d %s %s %s %d  %5.3f %5.3f %5.3f\n",
			   x, i,
			   str->atom[i].atomNumber,
			   str->atom[i].atomName,
			   str->atom[i].residueName,
			   str->atom[i].chainIdentifier,
			   str->atom[i].residueNumber,
			   str->atom[i].pos.x,
			   str->atom[i].pos.y,
			   str->atom[i].pos.z);
		*/

		/* map copied atom index i to original mmCIF atom index/count */
		str->atomMap[i] = str->nAllAtom;

		++i;
		++str->nAllAtom;
	}

	str->nAtom = i;
	str->nResidue = k;
	str->sequence.res[k] = '\0';

	/*____________________________________________________________________________*/
	/* free the compiled regular expressions */
	free_patterns(regexPattern, nHetAtom);
	free(regexPattern);
    free_structure(s);


    return 0;
}


