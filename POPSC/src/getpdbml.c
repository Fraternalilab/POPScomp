/*==============================================================================
getpdbml.c : routines for reading PDBML structures
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/
#include "getpdbml.h"
#include "pdb_structure.h"

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

	/* three-letter code of amino acid residues, exception HET -> X */
	char *aa3[] = {"ALA","---","CYS","ASP","GLU","PHE","GLY","HIS","ILE","---","LYS","LEU","MET","ASN","---","PRO","GLN","ARG","SER","THR","UNL","VAL","TRP","HET","TYR","UNK"};
	/* nucleotide residues */
	char *nuc[] = {"  A"," DA","  C"," DC","---","---","  G"," DG","  I"," DI","---","---","---","  N"," DN","---","---","---"," DT","  T","  U"," DU","---","---","---","---"};

	/* match against amino acid residues */
	residue = scan_array(code3, aa3, 65);

	/* match against nucleotide residues */
	if (residue == ' ') {
		residue = scan_array(code3, nuc, 97);
	} else if (residue == ' ') {
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
/* initialise all entries of a given atom */
__inline__ static void init_atom(Str *pdb)
{
	pdb->atom[pdb->nAtom].temperatureFactor = 0.;
	pdb->atom[pdb->nAtom].pos.x = 0.;
	pdb->atom[pdb->nAtom].pos.y = 0.;
	pdb->atom[pdb->nAtom].pos.z = 0.;
	strcpy(pdb->atom[pdb->nAtom].chainIdentifier, "");
	strcpy(pdb->atom[pdb->nAtom].atomName, "");
	strcpy(pdb->atom[pdb->nAtom].atomNameHet, "");
	strcpy(pdb->atom[pdb->nAtom].residueName, "");
	pdb->atom[pdb->nAtom].residueNumber = 0;
	strcpy(pdb->atom[pdb->nAtom].recordName, "");
	strcpy(pdb->atom[pdb->nAtom].icode, "");
	pdb->atom[pdb->nAtom].occupancy = 0.;
	pdb->atom[pdb->nAtom].modelNumber = 0;
	strcpy(pdb->atom[pdb->nAtom].element, "");
	strcpy(pdb->atom[pdb->nAtom].charge, "");
	pdb->atom[pdb->nAtom].formalCharge = 0;
	pdb->atom[pdb->nAtom].partialCharge = 0.;
	pdb->atom[pdb->nAtom].het = 0;
}

/*____________________________________________________________________________*/
/* The XML library transparently handles compression when doing
     file-based accesses. That is different from the 'read_structure' routine
     in 'getpdb', where the 'gz' library is being invoked explicitly. */
int parseXML(const char *filename, Str *pdb) {
    xmlDoc *doc; /* the resulting document tree */
    xmlNode *root_node = 0;
	xmlNode *cur_node = 0;
	xmlNode *site_node = 0;
	xmlNode *atom_node = 0;
	unsigned int allocated_atom = 64;
	unsigned int allocated_residue = 64;
	xmlChar *content = 0;
	unsigned int k = 0;
	int ca_p = 0;
	char resbuf;
	char line[80];
	regex_t *regexPattern = 0; /* regular atom patterns */
	/* allowed HETATM atom types (standard N,CA,C,O) and elements (any N,C,O,P,S) */
	const int nHetAtom = 6;
	char hetAtomPattern[6][32] = {{"N"},{"CA"},{"C"},{"O"},{"P"},{"S"}};
	char hetAtomNewname[6][32] = {{"N_"},{"CA"},{"C_"},{"O_"},{"P_"},{"S_"}};

	/*____________________________________________________________________________*/
    /* parse the file and get the document (DOM) */
	if ((doc = xmlReadFile(filename, NULL, 0)) == NULL) {
        fprintf(stderr, "XML Parser: Failed to read %s\n", filename);
		exit(1);
	}

	/*____________________________________________________________________________*/
	/* parse document tree */
	/* set root node */
	root_node = xmlDocGetRootElement(doc);

	/* atom sites are set as child nodes */
	for (cur_node = root_node->children; cur_node; cur_node = cur_node->next) {
		if(strcmp("atom_siteCategory", (char *)cur_node->name) == 0) {
			site_node = cur_node;
			break;
		}
	}

	/*____________________________________________________________________________*/
	/* allocate PDB structure */
	pdb->atom = safe_malloc(allocated_atom * sizeof(Atom));
	pdb->atomMap = safe_malloc(allocated_atom * sizeof(int));
	/* array of residue-centric atom indices */
	pdb->resAtom = safe_malloc(allocated_residue * sizeof(int));
	/* allocate memory for sequence residues */
	pdb->sequence.res = safe_malloc(allocated_residue * sizeof(char));

	/* initialise indices */
	pdb->nAtom = 0;
	pdb->nAllAtom = 0;
	pdb->nResidue = 0;
	pdb->nAllResidue = 0;
	pdb->nChain = 0;

	/* compile allowed HETATM element patterns */
	regexPattern = safe_malloc(nHetAtom * sizeof(regex_t));
	compile_patterns(regexPattern, &(hetAtomPattern[0]), nHetAtom);

	/*____________________________________________________________________________*/
	/* traverse XML tree (atom sites) */
	for (atom_node = site_node->children; atom_node; atom_node = atom_node->next) {
		if (strcmp("atom_site", (char *)atom_node->name) == 0) {
			/* initialise all entries of this atom */ 
			init_atom(pdb);

			/* children (= entries) of this atom site */
			for (cur_node = atom_node->children; cur_node; cur_node = cur_node->next) {
				/* assign node content to string */
				content = xmlNodeGetContent(cur_node);

				/* copy string content to PDB data structure */
				/* temperature factor */
				if (strcmp((char *)cur_node->name, "B_iso_or_equiv") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].temperatureFactor));
				}
				/* x coordinate */
				if (strcmp((char *)cur_node->name, "Cartn_x") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].pos.x));
				}
				/* y coordinate */
				if (strcmp((char *)cur_node->name, "Cartn_y") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].pos.y));
				}
				/* z coordinate */
				if (strcmp((char *)cur_node->name, "Cartn_z") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].pos.z));
				}
				/* chain identifier */
				if (strcmp((char *)cur_node->name, "auth_asym_id") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].chainIdentifier);
				}
				/* atom name */
				if (strcmp((char *)cur_node->name, "auth_atom_id") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].atomName);
				}
				/* residue name */
				if (strcmp((char *)cur_node->name, "auth_comp_id") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].residueName);
				}
				/* residue number */
				if (strcmp((char *)cur_node->name, "auth_seq_id") == 0) {
					sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].residueNumber));
				}
				/* ATOM or HETATM */
				if (strcmp((char *)cur_node->name, "group_PDB") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].recordName);
				}
				/* insert code */
				if (strcmp((char *)cur_node->name, "pdbx_PDB_ins_code") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].icode);
				}
				/* occupancy */
				if (strcmp((char *)cur_node->name, "occupancy") == 0) {
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].occupancy));
				}
				/* model number */
				if (strcmp((char *)cur_node->name, "pdbx_PDB_model_num") == 0) {
					sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].modelNumber));
				}
				/* atom element */
				if (strcmp((char *)cur_node->name, "type_symbol") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].element);
				}
				/* charge */
				if (strcmp((char *)cur_node->name, "pdbx_formal_charge") == 0) {
					sscanf((char *)content, "%s", pdb->atom[pdb->nAtom].charge);
					sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].formalCharge));
					sscanf((char *)content, "%f", &(pdb->atom[pdb->nAtom].partialCharge));
				}
			}

			/*____________________________________________________________________________*/
			/* select entries to record */
			/* only first MODEL if several are persent in PDB entry */
			if (pdb->atom[pdb->nAtom].modelNumber > 1) {
				goto ENDPARSE;
			}

			/* skip hydrogen atoms */
			if (strcmp(pdb->atom[pdb->nAtom].element, "H") == 0) {
				++ pdb->nAllAtom;
				continue;
			}

			/* aa code */
			if (strcmp(pdb->atom[pdb->nAtom].recordName, "ATOM") == 0) {
				assert((resbuf = aacode(pdb->atom[pdb->nAtom].residueName)) != ' ');
			}
	
			/* process HETATM entries */
			if (strcmp(pdb->atom[pdb->nAtom].recordName, "HETATM") == 0) {
				if (process_het(pdb, &(line[0]), regexPattern,
						&(hetAtomNewname[0]), nHetAtom) != 0) {
					continue;
				}
			}

			/* detect CA and N3 atoms for residue allocation */
			if ((strcmp(pdb->atom[pdb->nAtom].atomName, "CA") == 0) ||
				(strcmp(pdb->atom[pdb->nAtom].atomName, "N3") == 0)) {
				pdb->resAtom[k] = pdb->nAtom;
				pdb->sequence.res[k ++] = aacode(pdb->atom[pdb->nAtom].residueName);
				if (k == allocated_residue) {
					allocated_residue += 64;
					pdb->resAtom = safe_realloc(pdb->resAtom, allocated_residue * sizeof(int));
					pdb->sequence.res = safe_realloc(pdb->sequence.res, allocated_residue * sizeof(char));
				}
				++ ca_p;
			}

			/* standardise non-standard atom names (here GRO ILE_CD) */
			standardise_name(pdb->atom[pdb->nAtom].residueName,
								pdb->atom[pdb->nAtom].atomName);
			

			/*____________________________________________________________________________*/
			/* count number of allResidues (including HETATM residues) */
			if (pdb->nAtom == 0 ||
				pdb->atom[pdb->nAtom].residueNumber != pdb->atom[pdb->nAtom - 1].residueNumber ||
				strcmp(pdb->atom[pdb->nAtom].icode, pdb->atom[pdb->nAtom - 1].icode) != 0) {
				++ pdb->nAllResidue;
			}

			/*____________________________________________________________________________*/
			/* count number of chains */
			if (pdb->nAtom == 0 ||
				pdb->atom[pdb->nAtom].chainIdentifier[0] != pdb->atom[pdb->nAtom - 1].chainIdentifier[0]) {
				++ pdb->nChain;
			}

			/* records original atom order (count) */
			pdb->atomMap[pdb->nAtom] = pdb->nAllAtom;

			/* increment atom number */
			++ pdb->nAtom;
			++ pdb->nAllAtom;

			/* allocate more memory if needed */
			if (pdb->nAtom == allocated_atom) {
				allocated_atom += 64;
				pdb->atom = safe_realloc(pdb->atom, allocated_atom * sizeof(Atom));
				pdb->atomMap = safe_realloc(pdb->atomMap, allocated_atom * sizeof(int));
			}
		}
	}

	ENDPARSE:
	/* fprintf(stderr, "\tUsing only the first model of the PDB entry\n"); */

	pdb->nResidue = k;

	/*____________________________________________________________________________*/
	/* free global variables */
    xmlFreeDoc(doc);

	/*____________________________________________________________________________*/
	/* free the compiled regular expressions */
	free_patterns(regexPattern, nHetAtom);
	free(regexPattern);

	return 0;
}

/*____________________________________________________________________________*/
/* parse PDB file in XML format */
void read_structure_xml(Arg *arg, Argpdb *argpdb, Str *pdb)
{
    LIBXML_TEST_VERSION;

    pdb->sequence.name = safe_malloc((strlen(basename(arg->pdbmlInFileName)) + 1) * sizeof(char));
    strcpy(pdb->sequence.name, basename(arg->pdbmlInFileName));

    parseXML(arg->pdbmlInFileName, pdb);
    xmlCleanupParser();
    xmlMemoryDump();

    /* check for empty pdb structure and exit */
    if (pdb->nAtom == 0)
    {
        ErrorSpecNoexit("Invalid PDB file", arg->pdbmlInFileName);
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
							arg->pdbmlInFileName, pdb->nAllAtom,
							pdb->nAtom, pdb->nResidue, pdb->nChain);
}

