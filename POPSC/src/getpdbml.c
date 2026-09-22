/*==============================================================================
getpdbml.c : routines for reading PDBML structures
https://doi.org/10.1093/bioinformatics/bti082
http://xmlsoft.org/index.html
Copyright (C) 2018-2019 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/
#include "getpdbml.h"
#include "pdb_structure.h"

/* aacode() is shared with the PDB reader (getpdb.c) */

/*____________________________________________________________________________*/
/** standardise non-standard atom names */
__inline__ static int standardise_name(char *residueName, char *atomName)
{
	/* GRO 'ILE CD' to PDB 'ILE CD1' */
	/* atom names are unpadded in PDBML */
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
/* copy node content into a fixed-size field: leading/trailing blanks removed,
	truncated to the field size; empty or missing content leaves the field unchanged */
static void copy_content(char *dst, size_t size, const xmlChar *content)
{
	const char *c = (const char *)content;
	size_t n;

	if (c == NULL)
		return;
	while (*c == ' ' || *c == '\t' || *c == '\n' || *c == '\r')
		++ c;
	n = strcspn(c, " \t\n\r");
	if (n == 0)
		return;
	if (n >= size)
		n = size - 1;
	memcpy(dst, c, n);
	dst[n] = '\0';
}

/*____________________________________________________________________________*/
/* initialise all entries of a given atom */
__inline__ static void init_atom(Str *pdb)
{
	memset(&(pdb->atom[pdb->nAtom]), 0, sizeof(Atom));
	/* setting this to 'space' (ASCII 32) to check later whether
	   it is still that or another alternative location */
	strcpy(pdb->atom[pdb->nAtom].alternativeLocation, " ");
	/* no insertion code is written as '-', as in the PDB reader */
	strcpy(pdb->atom[pdb->nAtom].icode, "-");
}

/*____________________________________________________________________________*/
/* The XML library transparently handles compression when doing
     file-based accesses. That is different from the 'read_structure' routine
     in 'getpdb', where the 'gz' library is being invoked explicitly. */
int parseXML(const char *filename, Argpdb *argpdb, Str *pdb) {
    xmlDoc *doc; /* the resulting document tree */
    xmlNode *root_node = 0;
	xmlNode *cur_node = 0;
	xmlNode *site_node = 0;
	xmlNode *atom_node = 0;
	unsigned int allocated_atom = 64;
	unsigned int allocated_residue = 64;
	xmlChar *data = (xmlChar*)"datablockName";
	xmlChar *id = (xmlChar*)"id";
	xmlChar *content = 0;
	unsigned int k = 0;
	int recordIndex = 0; /* index of the atom_site record within the model */
	char residueAltloc = 0; /* altloc selected for the current residue */
	Atom *prev = 0; /* previous atom_site record, for altloc selection */
	Atom prevRecord;
	regex_t *regexPattern = 0; /* regular atom patterns */
	/* allowed HETATM atom types (standard N,CA,C,O) and elements (any N,C,O,P,S) */
	const int nHetAtom = 6;
	char hetAtomPattern[6][32] = {{"N"},{"CA"},{"C"},{"O"},{"P"},{"S"}};

	/*____________________________________________________________________________*/
    /* parse the file and get the document (DOM) */
	if ((doc = xmlReadFile(filename, NULL, 0)) == NULL) {
        fprintf(stderr, "XML Parser: Failed to read %s\n", filename);
		exit(1);
	}

	/*____________________________________________________________________________*/
	/* parse document tree */
	/* set root node */
	if ((root_node = xmlDocGetRootElement(doc)) == NULL) {
        fprintf(stderr, "XML Parser: No root element in %s\n", filename);
		exit(1);
	}
	/* extract pdbID via root node property "datablockName" */
	strcpy(pdb->pdbID, "");
	if ((content = xmlGetProp(root_node, data)) != NULL) {
		copy_content(pdb->pdbID, sizeof(pdb->pdbID), content);
		/* PDBML datablock names carry the suffix '-noatom' etc.: keep them as they are */
		xmlFree(content);
	}

	/* traverse XML tree: read datablock content and halt at atom_siteCategory */
	/* traverse atom sites as child nodes below */
	for (cur_node = root_node->children; cur_node; cur_node = cur_node->next) {
		if (cur_node->name && strcmp("atom_siteCategory", (char *)cur_node->name) == 0) {
			site_node = cur_node;
			break;
		}
	}
	if (site_node == NULL) {
        fprintf(stderr, "XML Parser: No atom_siteCategory in %s\n", filename);
		exit(1);
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
	/* traverse XML tree: atom sites */
	for (atom_node = site_node->children; atom_node; atom_node = atom_node->next) {
		if (atom_node->name && strcmp("atom_site", (char *)atom_node->name) == 0) {
			/* initialise all entries of this atom */ 
			init_atom(pdb);

			/* extract atom number via atom_site "id" */
			if ((content = xmlGetProp(atom_node, id)) != NULL) {
				sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].atomNumber));
				xmlFree(content);
			}
			
			/* children (= entries) of this atom site */
			for (cur_node = atom_node->children; cur_node; cur_node = cur_node->next) {
				if (cur_node->name == NULL)
					continue;
				/* assign node content to string */
				if ((content = xmlNodeGetContent(cur_node)) == NULL)
					continue;

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
					copy_content(pdb->atom[pdb->nAtom].chainIdentifier,
						sizeof(pdb->atom[pdb->nAtom].chainIdentifier), content);
				}
				/* atom name */
				if (strcmp((char *)cur_node->name, "auth_atom_id") == 0) {
					copy_content(pdb->atom[pdb->nAtom].atomName,
						sizeof(pdb->atom[pdb->nAtom].atomName), content);
				}
				/* alternative location */
				if (strcmp((char *)cur_node->name, "label_alt_id") == 0) {
					copy_content(pdb->atom[pdb->nAtom].alternativeLocation,
						sizeof(pdb->atom[pdb->nAtom].alternativeLocation), content);
				}
				/* residue name */
				if (strcmp((char *)cur_node->name, "auth_comp_id") == 0) {
					copy_content(pdb->atom[pdb->nAtom].residueName,
						sizeof(pdb->atom[pdb->nAtom].residueName), content);
				}
				/* residue number */
				if (strcmp((char *)cur_node->name, "auth_seq_id") == 0) {
					sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].residueNumber));
				}
				/* ATOM or HETATM */
				if (strcmp((char *)cur_node->name, "group_PDB") == 0) {
					copy_content(pdb->atom[pdb->nAtom].recordName,
						sizeof(pdb->atom[pdb->nAtom].recordName), content);
				}
				/* insert code */
				if (strcmp((char *)cur_node->name, "pdbx_PDB_ins_code") == 0) {
					copy_content(pdb->atom[pdb->nAtom].icode,
						sizeof(pdb->atom[pdb->nAtom].icode), content);
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
					copy_content(pdb->atom[pdb->nAtom].element,
						sizeof(pdb->atom[pdb->nAtom].element), content);
				}
				/* charge */
				if (strcmp((char *)cur_node->name, "pdbx_formal_charge") == 0) {
					copy_content(pdb->atom[pdb->nAtom].charge,
						sizeof(pdb->atom[pdb->nAtom].charge), content);
					sscanf((char *)content, "%d", &(pdb->atom[pdb->nAtom].formalCharge));
				}

				xmlFree(content);
			}

			/*____________________________________________________________________________*/
			/* select entries to record */
			/* only first MODEL if several are present in PDB entry */
			if (pdb->atom[pdb->nAtom].modelNumber > 1) {
				break;
			}

			/* index of this record in the model (trajectories list all records) */
			pdb->atom[pdb->nAtom].recordIndex = recordIndex ++;

			/* select one alternative location per residue: the first one encountered */
			if (prev == 0 ||
				strcmp(pdb->atom[pdb->nAtom].chainIdentifier, prev->chainIdentifier) != 0 ||
				pdb->atom[pdb->nAtom].residueNumber != prev->residueNumber ||
				strcmp(pdb->atom[pdb->nAtom].icode, prev->icode) != 0) {
				residueAltloc = 0;
			}
			prevRecord = pdb->atom[pdb->nAtom];
			prev = &prevRecord;
			if (pdb->atom[pdb->nAtom].alternativeLocation[0] != ' ' &&
				pdb->atom[pdb->nAtom].alternativeLocation[0] != '.' &&
				pdb->atom[pdb->nAtom].alternativeLocation[0] != '?') {
				if (residueAltloc == 0)
					residueAltloc = pdb->atom[pdb->nAtom].alternativeLocation[0];
				if (pdb->atom[pdb->nAtom].alternativeLocation[0] != residueAltloc)
					continue;
			}

			/* skip hydrogen and deuterium atoms unless requested */
			if (! argpdb->hydrogens &&
				(strcmp(pdb->atom[pdb->nAtom].element, "H") == 0 ||
				 strcmp(pdb->atom[pdb->nAtom].element, "D") == 0)) {
				continue;
			}

			/* HETATM entries are not processed */
			if (strcmp(pdb->atom[pdb->nAtom].recordName, "ATOM") != 0) {
				continue;
			}

			/* water is not part of the solute surface */
			if (is_water(pdb->atom[pdb->nAtom].residueName)) {
				continue;
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
			}

			/* standardise non-standard atom names (here GRO ILE_CD) */
			standardise_name(pdb->atom[pdb->nAtom].residueName,
								pdb->atom[pdb->nAtom].atomName);

			/* in coarse mode record only the representative atoms: CA (amino acids), P (nucleotides) */
			if (argpdb->coarse &&
				! is_coarse_atom(pdb->atom[pdb->nAtom].residueName, pdb->atom[pdb->nAtom].atomName)) {
				continue;
			}

			/* records original atom order (count) */
			pdb->atomMap[pdb->nAtom] = pdb->atom[pdb->nAtom].recordIndex;

			/* increment atom number */
			++ pdb->nAtom;

			/* allocate more memory if needed */
			if (pdb->nAtom == allocated_atom) {
				allocated_atom += 64;
				pdb->atom = safe_realloc(pdb->atom, allocated_atom * sizeof(Atom));
				pdb->atomMap = safe_realloc(pdb->atomMap, allocated_atom * sizeof(int));
			}
		}
	}

	pdb->nResidue = k;
	pdb->sequence.res[k] = '\0';
	/* all atom_site records of the model: the atom count of trajectory frames */
	pdb->nAllAtom = recordIndex;

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

    parseXML(arg->pdbmlInFileName, argpdb, pdb);
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

	/* residue and chain indices, residue and chain counts */
	index_structure(pdb);
	/* identifier for output file names: datablock name, else file name */
	if (strlen(pdb->pdbID) == 0)
		set_pdbID_from_filename(pdb, arg->pdbmlInFileName);

    if (! arg->silent) fprintf(stdout, "\tPDB file: %s\n"
										"\tPDB identifier: %s\n"
										"\tPDB file content:\n"
										"\t\tall atoms = %d\n"
										"\t\tprocessed atoms (C,N,O,S,P) = %d\n"
										"\t\tresidues (CA||N3) = %d\n"
										"\t\tchains = %d\n",
							arg->pdbmlInFileName, pdb->pdbID, pdb->nAllAtom,
							pdb->nAtom, pdb->nResidue, pdb->nChain);
}

