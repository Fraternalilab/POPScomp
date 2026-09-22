/*==============================================================================
getpdb.c : routines for reading PDB structures
Copyright (C) 2004-2026 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "getpdb.h"

/*____________________________________________________________________________*/
/* remove spaces from atom and residue names */
void remove_spaces(char *s)
{
    char *src = s;
    char *dst = s;

    while (*src) {
        if (*src != ' ')
            *dst++ = *src;
        src++;
    }

    *dst = '\0';
}

/*____________________________________________________________________________*/
/* match PDB residue name against constant residue name array */
__inline__ static char scan_array(char *code3, char *residue_array[], int shift)
{
	unsigned int i;
	char residue = ' ';

	for (i = 0; i < 26; ++ i)
		if (strcmp(code3, residue_array[i]) == 0) {
			residue = i + shift; /* shift=65 for UPPER, shift=97 for lower */
			break;
		}

	return residue;
}

/*____________________________________________________________________________*/
/** amino acid 3-letter to 1-letter code conversion */
char aacode(char *code3)
{
	char residue = ' '; /* 1-letter residue name */

	/* three-letter code of amino acid residues, exception HET (X) */
	char *aa3[] = {"ALA","---","CYS","ASP","GLU","PHE","GLY","HIS","ILE","---","LYS","LEU","MET","ASN","---","PRO","GLN","ARG","SER","THR","UNL","VAL","TRP","HET","TYR","UNK"};
	/* nucleotide residues */
	/* residue names are space-stripped on input: compare unpadded names */
	char *nuc[] = {"A","DA","C","DC","---","---","G","DG","I","DI","---","---","---","N","DN","---","---","---","DT","T","U","DU","---","---","---","---"};

	/* match against amino acid residues */
	residue = scan_array(code3, aa3, 65);

	/* match against nucleotide residues */
	if (residue == ' ')
		residue = scan_array(code3, nuc, 97);

	/* residue not found */
	if (residue == ' ')
		residue = 'X';
	
	return residue;
}

/*____________________________________________________________________________*/
/** water residue names; water is not part of the solute surface */
int is_water(const char *residueName)
{
	return (strcmp(residueName, "HOH") == 0 || strcmp(residueName, "WAT") == 0 ||
			strcmp(residueName, "DOD") == 0);
}

/*____________________________________________________________________________*/
/** nucleotide residue names */
static int is_nucleotide(const char *residueName)
{
	const char *nuc[] = {"A","C","G","I","N","T","U","DA","DC","DG","DI","DN","DT","DU"};
	unsigned int i;

	for (i = 0; i < sizeof(nuc) / sizeof(nuc[0]); ++ i)
		if (strcmp(residueName, nuc[i]) == 0)
			return 1;
	return 0;
}

/*____________________________________________________________________________*/
/** representative atom of a residue in coarse-grained mode:
	P for nucleotides, CA for all other residues (see coarse table in sasa_data.h) */
int is_coarse_atom(const char *residueName, const char *atomName)
{
	if (is_nucleotide(residueName))
		return (strcmp(atomName, "P") == 0);
	return (strcmp(atomName, "CA") == 0);
}

/*____________________________________________________________________________*/
/** map common force-field residue names onto the parametrised standard residue */
static void standardise_residue(Atom *atom)
{
	unsigned int i;
	static int warned = 0;
	const char *alias[][2] = {
		{"HIE","HIS"},{"HID","HIS"},{"HIP","HIS"},{"HSD","HIS"},{"HSE","HIS"},{"HSP","HIS"},
		{"CYX","CYS"},{"CYM","CYS"},{"ASH","ASP"},{"GLH","GLU"},{"LYN","LYS"},{"MSE","MET"}};

	/* GROMOS/GROMACS water oxygen: parametrised as 'O2' of residue SOL */
	if (strcmp(atom->residueName, "SOL") == 0 && strcmp(atom->atomName, "OW") == 0) {
		strcpy(atom->atomName, "O2");
		return;
	}

	for (i = 0; i < sizeof(alias) / sizeof(alias[0]); ++ i) {
		if (strcmp(atom->residueName, alias[i][0]) == 0) {
			if (! warned) {
				fprintf(stderr, "Warning: mapping force-field residue names onto standard residues"
								" (first: %s -> %s)\n", alias[i][0], alias[i][1]);
				warned = 1;
			}
			strcpy(atom->residueName, alias[i][1]);
			/* selenomethionine selenium takes the place of the methionine sulphur */
			if (strcmp(alias[i][0], "MSE") == 0 && strcmp(atom->atomName, "SE") == 0)
				strcpy(atom->atomName, "SD");
			return;
		}
	}
}

/*____________________________________________________________________________*/
/** assign running residue and chain indices and count residues and chains;
	one definition for all input formats: a new residue starts when the chain,
	the residue number or the insertion code changes, a new chain starts
	when the chain identifier changes */
void index_structure(Str *str)
{
	int i;

	str->nAllResidue = 0;
	str->nChain = 0;

	for (i = 0; i < str->nAtom; ++ i) {
		standardise_residue(&(str->atom[i]));

		if (i == 0 ||
			strcmp(str->atom[i].chainIdentifier, str->atom[i - 1].chainIdentifier) != 0) {
			++ str->nChain;
		}
		str->atom[i].chainIndex = str->nChain - 1;

		if (i == 0 ||
			str->atom[i].chainIndex != str->atom[i - 1].chainIndex ||
			str->atom[i].residueNumber != str->atom[i - 1].residueNumber ||
			strcmp(str->atom[i].icode, str->atom[i - 1].icode) != 0) {
			++ str->nAllResidue;
		}
		str->atom[i].residueIndex = str->nAllResidue - 1;
	}
}

/*____________________________________________________________________________*/
/** generic atom name of the UNK residue parameters for an atom of an unknown residue,
	chosen by element: N, O, S and P map onto themselves, carbon onto the side-chain 'CG';
	returns NULL for other elements */
const char *unk_atom_name(Atom *atom)
{
	char e = '\0';
	unsigned int i;

	/* element column if present, else the first letter of the atom name */
	if (strlen(atom->element) == 1)
		e = atom->element[0];
	else if (strlen(atom->element) == 0)
		for (i = 0; i < strlen(atom->atomName); ++ i)
			if (isalpha((unsigned char)atom->atomName[i])) {
				e = atom->atomName[i];
				break;
			}

	switch (e) {
		case 'C': return "CG";
		case 'N': return "N";
		case 'O': return "O";
		case 'S': return "S";
		case 'P': return "P";
		default: return NULL;
	}
}

/*___________________________________________________________________________*/
/** chain identifier for output: an empty (blank) chain identifier is written as '-'
	for compatibility with POPSR Shiny and whitespace-separated parsers */
const char *chain_label(Atom *atom)
{
	return (strlen(atom->chainIdentifier) == 0) ? "-" : atom->chainIdentifier;
}

/*____________________________________________________________________________*/
/** exit on a structure that cannot be processed;
	in JSON mode an empty JSON file is written and the exit status is 0
	(server behaviour), otherwise the exit status is 1 */
void exit_structure_error(Arg *arg, Str *str, const char *message)
{
	char name[512];
	FILE *jsonFile;

	fprintf(stderr, "Error: %s\n", message);

	if (arg->jsonOut) {
		snprintf(name, sizeof(name), "%s.json", str->pdbID);
		jsonFile = open_output(arg, name, "w");
		fclose(jsonFile);
		exit(0);
	}
	exit(1);
}

/*____________________________________________________________________________*/
/** set the structure identifier from the input file name, without extensions */
void set_pdbID_from_filename(Str *str, const char *fileName)
{
	const char *base = strrchr(fileName, '/');
	char *dot;

	base = base ? base + 1 : fileName;
	snprintf(str->pdbID, sizeof(str->pdbID), "%s", base);
	/* strip compression and format extensions */
	while ((dot = strrchr(str->pdbID, '.')) != NULL && dot != str->pdbID &&
		(strcmp(dot, ".gz") == 0 || strcmp(dot, ".pdb") == 0 || strcmp(dot, ".ent") == 0 ||
		 strcmp(dot, ".cif") == 0 || strcmp(dot, ".mmcif") == 0 || strcmp(dot, ".xml") == 0))
		*dot = '\0';
}

/*____________________________________________________________________________*/
/** standardise non-standard atom names */
__inline__ static int standardise_name(char *residueName, char *atomName)
{
	/* GRO 'ILE CD' to PDB 'ILE CD1' */
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
/** process HET residues */
/*
__inline__ static int process_het(Str *str, char *line, regex_t *regexPattern, char (*hetAtomNewname)[32], int nHetAtom)
{
	int hetAtomNr = -1;

	if ((hetAtomNr = match_patterns(regexPattern, nHetAtom, &(str->atom[str->nAtom].atomName[0]))) >= 0) {
		sprintf(str->atom[str->nAtom].atomName, "%s", &(hetAtomNewname[hetAtomNr][0]));
		sprintf(str->atom[str->nAtom].residueName, "%s", "HET");
		fprintf(stderr, "Setting atom %d name %s to %s of residue HET\n",
					str->nAtom, str->atom[str->nAtom].atomName, 
					&(hetAtomNewname[hetAtomNr][0]));
	} else {
		WarningSpec("Skipping HETATM", str->atom[str->nAtom].atomName);
		return 1;
	}

	return 0;
}
*/

/*____________________________________________________________________________*/
/** read PDB file */

/**
http://www.wwpdb.org/docs.html
The ATOM record:

COLUMNS        DATA TYPE       FIELD         DEFINITION
________________________________________________________________________________
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X
                                               in Angstroms.
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
                                               in Angstroms.
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
                                               in Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier,
                                               left-justified.
77 - 78        LString(2)      element       Element symbol,
                                               right-justified.
79 - 80        LString(2)      charge        Charge on the atom.
*/

int read_pdb(FILE *pdbInFile, gzFile *pdbgzInFile, Arg *arg, Argpdb *argpdb, Str *str)
{
	unsigned int i, j;
	unsigned int k = 0;
	char line[256];
	unsigned int allocated_atom = 64;
	unsigned int allocated_residue = 64;
	char atomName[8] = "";
	int nModel = 0; /* number of MODEL records seen */
	int recordIndex = 0; /* index of ATOM/HETATM record within the model */
	int prevResidueNumber = 0; /* residue key of the previous record, for altloc selection */
	char prevIcode = 0;
	char prevChain = 0;
	char residueAltloc = 0; /* altloc selected for the current residue */
	/* for HETATM entries */
	regex_t *regexPattern = 0; /* regular atom patterns */
	/* allowed HETATM atom types (standard N,CA,C,O) and elements (any N,C,O,P,S) */
	const int nHetAtom = 9;
	char hetAtomPattern[9][32] = {{"N"},{"CA"},{"C"},{"O"},{".{1}C[[:print:]]{1,3}"},{".{1}N[[:print:]]{1,3}"},{".{1}O[[:print:]]{1,3}"},{".{1}P[[:print:]]{1,3}"},{".{1}S[[:print:]]{1,3}"}};

	/*____________________________________________________________________________*/
	/* initialise/allocate memory for set of (64) selected (CA) atom entries */
	str->nAtom = 0;
	str->nAllAtom = 0;
	str->nResidue = 0;
	str->nAllResidue = 0;
	str->nChain = 0;
	
	str->atom = safe_malloc(allocated_atom * sizeof(Atom));
	str->atomMap = safe_malloc(allocated_atom * sizeof(int));

	/* array of residue-centric atom indices */
	str->resAtom = safe_malloc(allocated_residue * sizeof(int));

	/* allocate memory for sequence residues */
	str->sequence.res = safe_malloc(allocated_residue * sizeof(char));

	/* compile allowed HETATM element patterns */
	regexPattern = safe_malloc(nHetAtom * sizeof(regex_t));
	compile_patterns(regexPattern, &(hetAtomPattern[0]), nHetAtom);

	/*____________________________________________________________________________*/
	/* not all PDB data types are used in this program to save resources */
    while (1) {
		if (arg->zipped) {
			if (gzgets(*pdbgzInFile, line, sizeof(line)) == 0) {
				break;
			}
		} else {
			if (fgets(line, sizeof(line), pdbInFile) == 0) {
				break;
			}
		}

		/* discard the remainder of lines longer than the buffer */
		if (strchr(line, '\n') == NULL) {
			int ch;
			if (arg->zipped) {
				while ((ch = gzgetc(*pdbgzInFile)) != -1 && ch != '\n');
			} else {
				while ((ch = fgetc(pdbInFile)) != EOF && ch != '\n');
			}
		}
	

		/*____________________________________________________________________________*/
		/* check conditions to start assigning this entry */
		/* read only the first model: stop at its ENDMDL or at the second MODEL record */
		if (strncmp(line, "MODEL", 5) == 0) {
			if (++ nModel > 1)
				break;
			continue;
		}
		if (strncmp(line, "ENDMDL", 6) == 0) {
			break;
		}

		/* read only ATOM/HETATM records */
		if((strncmp(line, "ATOM  ", 6) != 0) && (strncmp(line, "HETATM", 6) != 0)) {
			continue;
		}

		/* pad short lines with spaces up to the charge column */
		for (i = strcspn(line, "\r\n"); i < 80; ++ i)
			line[i] = ' ';
		line[80] = '\0';

		/* index of this record in the model (trajectories list all records) */
		++ recordIndex;

		/* select one alternative location per residue: the first one encountered */
		if (line[21] != prevChain || atoi(&line[22]) != prevResidueNumber || line[26] != prevIcode) {
			prevChain = line[21];
			prevResidueNumber = atoi(&line[22]);
			prevIcode = line[26];
			residueAltloc = 0;
		}
		if (line[16] != ' ') {
			if (residueAltloc == 0)
				residueAltloc = line[16];
			if (line[16] != residueAltloc)
				continue;
		}

		/*____________________________________________________________________________*/
		/* read this entry */
		memset(&(str->atom[str->nAtom]), 0, sizeof(Atom));
		str->atom[str->nAtom].recordIndex = recordIndex - 1;

		/* atom number */
		str->atom[str->nAtom].atomNumber = atoi(&line[6]);

		/* record type */
		for (i = 0, j = 0; i < 6; ) {
			str->atom[str->nAtom].recordName[j++] = line[i++];
		}
		str->atom[str->nAtom].recordName[j] = '\0';
		
		/* atom name */
		for (i = 12, j = 0; i < 16; ) {
			str->atom[str->nAtom].atomName[j++] = line[i++];
		}
		str->atom[str->nAtom].atomName[j] = '\0';
		remove_spaces(str->atom[str->nAtom].atomName);

		/* alternative location */
		str->atom[str->nAtom].alternativeLocation[0] = line[16];
		str->atom[str->nAtom].alternativeLocation[1] = '\0';

		/* residue name */
		for (i = 17, j = 0; i < 20; ) {
			str->atom[str->nAtom].residueName[j++] = line[i++];
		}
		str->atom[str->nAtom].residueName[j] = '\0';
		remove_spaces(str->atom[str->nAtom].residueName);

		/* chain identifier */
		str->atom[str->nAtom].chainIdentifier[0] = line[21];
		str->atom[str->nAtom].chainIdentifier[1] = '\0';
		remove_spaces(str->atom[str->nAtom].chainIdentifier);

		/* residue number */
		str->atom[str->nAtom].residueNumber = atoi(&line[22]);

		/* code for insertion of residues */
		str->atom[str->nAtom].icode[0] = line[26];
		str->atom[str->nAtom].icode[1] = '\0';

		/* avoid space character in icode */
		if (isspace(str->atom[str->nAtom].icode[0]))
			str->atom[str->nAtom].icode[0] = '-';

		/* coordinates */
		str->atom[str->nAtom].pos.x = atof(&line[30]);
		str->atom[str->nAtom].pos.y = atof(&line[38]);
		str->atom[str->nAtom].pos.z = atof(&line[46]);

		/* element */
		for (i = 76, j = 0; i < 78; ) {
			str->atom[str->nAtom].element[j++] = line[i++];
		}
		str->atom[str->nAtom].element[j] = '\0';
		remove_spaces(str->atom[str->nAtom].element);

		/* description: everything before coordinates */
		for (i = 0, j = 0; i < 30; ) {
			str->atom[str->nAtom].description[j++] = line[i++];
		}
		str->atom[str->nAtom].description[j] = '\0';

		/*____________________________________________________________________________*/
		/* check conditions to record this entry */
		/* if no hydrogens set, skip hydrogen lines, including deuterium */
		if (! argpdb->hydrogens) {
			strip_char(str->atom[str->nAtom].atomName, &(atomName[0]));
			/* skip patterns 'H...' and '?H..', where '?' is a digit */
			if ((atomName[0] == 'H') || \
				((atomName[0] >= 48) && (atomName[0] <= 57) && (atomName[1] == 'H'))) {
				continue;
			}
			/* same for D (deuterium)*/
			if ((atomName[0] == 'D') || \
				((atomName[0] >= 48) && (atomName[0] <= 57) && (atomName[1] == 'D'))) {
				continue;
			}
		}

		/* HETATM entries are not processed */
		if (strncmp(line, "HETATM", 6) == 0) {
			continue;
		}

		/* water written as ATOM records is not processed either */
		if (is_water(str->atom[str->nAtom].residueName)) {
			continue;
		}

		/* detect CA and N3 atoms of standard residues for residue allocation */
		/* atomName has been space-stripped above: compare unpadded names */
		if ((strcmp(str->atom[str->nAtom].atomName, "CA") == 0) ||
		(strcmp(str->atom[str->nAtom].atomName, "N3") == 0)) {
			str->resAtom[k] = str->nAtom;
			str->sequence.res[k ++] = aacode(str->atom[str->nAtom].residueName);
			if (k == allocated_residue) {
				allocated_residue += 64;
				str->resAtom = safe_realloc(str->resAtom, allocated_residue * sizeof(int));
				str->sequence.res = safe_realloc(str->sequence.res, allocated_residue * sizeof(char));
			}
		}

		/* standardise non-standard atom names */
		standardise_name(str->atom[str->nAtom].residueName, str->atom[str->nAtom].atomName);

		/* in coarse mode record only the representative atoms: CA (amino acids), P (nucleotides) */
		if (argpdb->coarse && ! is_coarse_atom(str->atom[str->nAtom].residueName, str->atom[str->nAtom].atomName))
			continue;

		/*____________________________________________________________________________*/
		/* records original atom order (count) */
		str->atomMap[str->nAtom] = str->atom[str->nAtom].recordIndex;
		/* increment to next atom entry */
		++ str->nAtom;

		/*____________________________________________________________________________*/
		/* allocate more memory if needed */
		if (str->nAtom == allocated_atom) {
			allocated_atom += 64;
			str->atom = safe_realloc(str->atom, allocated_atom * sizeof(Atom));
			str->atomMap = safe_realloc(str->atomMap, allocated_atom * sizeof(int));
		}
	}
	str->sequence.res[k] = '\0';
	str->nResidue = k;
	/* all ATOM/HETATM records of the model: the atom count of trajectory frames */
	str->nAllAtom = recordIndex;

	/*____________________________________________________________________________*/
	/* free the compiled regular expressions */
	free_patterns(regexPattern, nHetAtom);
	free(regexPattern);

	/*____________________________________________________________________________*/
	return 0;
}

/*____________________________________________________________________________*/
/** read CONECT record */
int read_conect(FILE *pdbInFile)
{
	char line[80];

	if (fseek(pdbInFile, 0L, SEEK_SET) != 0) {
		fprintf(stderr, "File pointer repositioning error\n");
		exit(1);
	}
	
    while(fgets(line, 80, pdbInFile) != 0) {
        if (strncmp(line, "CONECT", 6) == 0) {
			;
		}
	}

	return 0;
}

/*_____________________________________________________________________________*/
/** read PDB structure */
void read_structure(Arg *arg, Argpdb *argpdb, Str *pdb)
{
	FILE *pdbInFile = 0;
	gzFile pdbgzInFile;

/*
	typedef voidp gzFile;

	gzFile gzopen (const char *path, const char *mode);
    Opens a gzip (.gz) file for reading or writing. The mode parameter is as in fopen ("rb" or "wb") but can also include a compression level ("wb9") or a strategy: 'f' for filtered data as in "wb6f", 'h' for Huffman only compression as in "wb1h". (See the description of deflateInit2 for more information about the strategy parameter.)

    gzopen can be used to read a file which is not in gzip format ; in this case gzread will directly read from the file without decompression.

    gzopen returns NULL if the file could not be opened or if there was insufficient memory to allocate the (de)compression state ; errno can be checked to distinguish the two cases (if errno is zero, the zlib error is Z_MEM_ERROR).
*/

    pdb->sequence.name = safe_malloc((strlen(basename(arg->pdbInFileName)) + 1) * sizeof(char));
    strcpy(pdb->sequence.name, basename(arg->pdbInFileName));

	/* gzipped or raw input file */
	/* passing both types of file pointers to read_pdb, but only one will be used */
	if (arg->zipped) {
		pdbgzInFile = gzopen(arg->pdbInFileName, "r");
		read_pdb(pdbInFile, &pdbgzInFile, arg, argpdb, pdb);
		gzclose(pdbgzInFile);
	} else {
		pdbInFile = safe_open(arg->pdbInFileName, "r");
		read_pdb(pdbInFile, &pdbgzInFile, arg, argpdb, pdb);
		fclose(pdbInFile);
	}

    /* check for empty pdb structure and exit */
    if (pdb->nAtom == 0) {
        ErrorSpecNoexit("Could not find atoms in input file. Please check input parameters:\nformat (--pdb or --pdbml), compression (--zipped or not) and file name",
			arg->pdbInFileName);
        free(pdb->atom);
        free(pdb->sequence.res);
        free(pdb->sequence.name);
        exit(1);
    }

	/* residue and chain indices, residue and chain counts */
	index_structure(pdb);
	/* identifier for output file names */
	set_pdbID_from_filename(pdb, arg->pdbInFileName);

    if (! arg->silent) fprintf(stdout, "\tPDB file: %s\n"
										"\tPDB file content:\n"
										"\t\tall atoms = %d\n"
										"\t\tprocessed atoms (C,N,O,S,P) = %d\n"
										"\t\tresidues (CA||N3) = %d\n"
										"\t\tchains = %d\n",
							arg->pdbInFileName, pdb->nAllAtom,
							pdb->nAtom, pdb->nResidue, pdb->nChain);
}

