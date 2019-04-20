/*==============================================================================
getpdb.c : routines for reading PDB structures
Copyright (C) 2004 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "getpdb.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

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
	char line[80];
	char stopline[80] = "";
    int stopflag = 0;
	unsigned int allocated_atom = 64;
	unsigned int allocated_residue = 64;
	char atomName[] = "    ";
	char resbuf;
	int ca_p = 0;
	/* for HETATM entries */
	regex_t *regexPattern = 0; /* regular atom patterns */
	/* allowed HETATM atom types (standard N,CA,C,O) and elements (any N,C,O,P,S) */
	const int nHetAtom = 9;
	char hetAtomPattern[9][32] = {{" N  "},{" CA "},{" C  "},{" O  "},{".{1}C[[:print:]]{1,3}"},{".{1}N[[:print:]]{1,3}"},{".{1}O[[:print:]]{1,3}"},{".{1}P[[:print:]]{1,3}"},{".{1}S[[:print:]]{1,3}"}};
	char hetAtomNewname[9][32] = {{" N  "},{" CA "},{" C  "},{" O  "},{" C_ "},{" N_ "},{" O_ "},{" P_ "},{" S_ "}};

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
    /* count the number of models */
	if (arg->zipped) {
		while(gzgets(*pdbgzInFile, line, 80) != 0) {
			if (strncmp(line, "MODEL ", 6) == 0) {
				if (stopflag == 0) {
					stopflag = 1;
					continue;
				} else {
					strcpy(stopline, line);
					break;
				}
			}
		}
	} else {
		while(fgets(line, 80, pdbInFile) != 0) {
			if (strncmp(line, "MODEL ", 6) == 0) {
				if (stopflag == 0) {
					stopflag = 1;
					continue;
				} else {
					strcpy(stopline, line);
					break;
				}
			}
		}
	}

    /* rewind the file handle to the start */
	if (arg->zipped) {
		if (gzseek(*pdbgzInFile, 0L, SEEK_SET) != 0) {
			/* handle repositioning error */
		}
	} else {
		if (fseek(pdbInFile, 0L, SEEK_SET) != 0) {
			/* handle repositioning error */
		}
	}

	/*____________________________________________________________________________*/
	/* not all PDB data types are used in this program to save resources */
    while (1) {
		if (arg->zipped) {
			if (gzgets(*pdbgzInFile, line, 80) == 0) {
				break;
			}
		} else {
			if (fgets(line, 80, pdbInFile) == 0) {
				break;
			}
		}
	
		ca_p = 0; /* CA or P flag */

		/*____________________________________________________________________________*/
		/* check conditions to start assigning this entry */
		/* skip other models */
		if((strcmp(line, stopline) == 0) && (stopflag == 1)) {
			break;
		}

		/* read only ATOM/HETATM records */
		if((strncmp(line, "ATOM  ", 6) != 0) && (strncmp(line, "HETATM", 6) != 0)) {
			continue;
		}

        /* skip alternative locations except for location 'A' */ 
		if (line[16] != 32 && line[16] != 65) {
			/*fprintf(stderr, "Warning: Skipping atom %d in alternative location %c\n",
				atoi(&line[6]), line[16]);*/
			continue;
		}

		/*____________________________________________________________________________*/
		/* read this entry */
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

		/* alternative location */
		/*str->atom[str->nAtom].alternativeLocation[0] = line[16];	
		str->atom[str->nAtom].alternativeLocation[1] = '\0';*/

		/* residue name */
		for (i = 17, j = 0; i < 20; ) {
			str->atom[str->nAtom].residueName[j++] = line[i++];
		}
		str->atom[str->nAtom].residueName[j] = '\0';

		/* chain identifier */
		str->atom[str->nAtom].chainIdentifier[0] = line[21];
		str->atom[str->nAtom].chainIdentifier[1] = '\0';

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

		/*printf("x %6.4f, y %6.4f, z %6.4f\n", str->atom[str->nAtom].x,
			str->atom[str->nAtom].y, str->atom[str->nAtom].z);*/

		/* occupancy */
		/*str->atom[str->nAtom].occupancy = atof(&line[54]);*/

		/* temperature factor */
		/*str->atom[str->nAtom].temp_f = atof(&line[60]);*/

		/* segment identifier */
		/*for (i = 72, j = 0; i < 76; )
			str->atom[str->nAtom].segmentIdentifier[j++] = line[i++];
		str->atom[str->nAtom].segmentIdentifier[j] = '\0';*/

		/* element */
		for (i = 76, j = 0; i < 78; ) {
			str->atom[str->nAtom].element[j++] = line[i++];
		}
		str->atom[str->nAtom].element[j] = '\0';

		/* charge */
		/*for (i = 78, j = 0; i < 80; )
			str->atom[str->nAtom].charge[j++] = line[i++];
		str->atom[str->nAtom].charge[j] = '\0';*/

		/* description: everything before coordinates */
		for (i = 0, j = 0; i < 30; ) {
			str->atom[str->nAtom].description[j++] = line[i++];
		}
		str->atom[str->nAtom].description[j] = '\0';

		/*____________________________________________________________________________*/
		/* check conditions to record this entry */
		/* if no hydrogens set, skip hydrogen lines */
		if (! argpdb->hydrogens) {
			strip_char(str->atom[str->nAtom].atomName, &(atomName[0]));
			/* skip patterns 'H...' and '?H..', where '?' is a digit */
			if ((atomName[0] == 'H') || \
				((atomName[0] >= 48) && (atomName[0] <= 57) && (atomName[1] == 'H'))) {
				++ str->nAllAtom;
				continue;
			}
		}

		/* aa code */
		if (strncmp(line, "ATOM  ", 6) == 0) {
			assert((resbuf = aacode(str->atom[str->nAtom].residueName)) != ' ');
		}
	
		/* process HETATM entries */
		if (strncmp(line, "HETATM", 6) == 0) {
			if (process_het(str, &(line[0]), regexPattern, &(hetAtomNewname[0]), nHetAtom) != 0) {
				continue;
			}
		}

		/* detect CA and N3 atoms of standard residues for residue allocation */
		if ((strncmp(str->atom[str->nAtom].atomName, " CA ", 4) == 0) ||
		(strncmp(str->atom[str->nAtom].atomName, " N3 ", 4) == 0)) {
			str->resAtom[k] = str->nAtom;
			str->sequence.res[k ++] = aacode(str->atom[str->nAtom].residueName);
			if (k == allocated_residue) {
				allocated_residue += 64;
				str->resAtom = safe_realloc(str->resAtom, allocated_residue * sizeof(int));
				str->sequence.res = safe_realloc(str->sequence.res, allocated_residue * sizeof(char));
			}
			++ ca_p;
		}

		/* standardise non-standard atom names */
		standardise_name(str->atom[str->nAtom].residueName, str->atom[str->nAtom].atomName);

		/* in coarse mode record only CA and P entries */
		if (!ca_p && argpdb->coarse)
			continue;

		/*____________________________________________________________________________*/
		/* count number of allResidues (including HETATM residues) */
        if (str->nAtom == 0 ||
			str->atom[str->nAtom].residueNumber != str->atom[str->nAtom - 1].residueNumber ||
			strcmp(str->atom[str->nAtom].icode, str->atom[str->nAtom - 1].icode) != 0) {
			++ str->nAllResidue;
		}

		/*____________________________________________________________________________*/
		/* count number of chains */
        if (str->nAtom == 0 ||
			str->atom[str->nAtom].chainIdentifier[0] != str->atom[str->nAtom - 1].chainIdentifier[0]) {
			++ str->nChain;
		}

		/*____________________________________________________________________________*/
		/* records original atom order (count) */
		str->atomMap[str->nAtom] = str->nAllAtom;
		/* increment to next atom entry */
		++ str->nAtom;
		++ str->nAllAtom;

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
    if (pdb->nAtom == 0)
    {
        ErrorSpecNoexit("Could not find atoms in input file. Please check input parameters:\nformat (--pdb or --pdbml), compression (--zipped or not) and file name",
			arg->pdbInFileName);
        free(pdb->atom);
        free(pdb->sequence.res);
        free(pdb->sequence.name);
        exit(1);
    }

    if (! arg->silent) fprintf(stdout, "\tPDB file: %s\n"
										"\tPDB file content:\n"
										"\t\tall atoms = %d\n"
										"\t\tprocessed atoms (C,N,O,S,P) = %d\n"
										"\t\tresidues (CA||N3) = %d\n"
										"\t\tchains = %d\n",
							arg->pdbInFileName, pdb->nAllAtom,
							pdb->nAtom, pdb->nResidue, pdb->nChain);
}

