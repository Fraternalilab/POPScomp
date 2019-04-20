/*==============================================================================
json.c : JSON routines, using the cJSON library
To validate JSON output format:
../funpdbe-client/funpdbe_client.py --path=pops.json --mode=validate
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "json.h"
#include "cJSON.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
void print_json(Arg *arg, cJSON *json)
{
	char outpath[256];
	/* print JSON object to string */
	char *popsOutJson = cJSON_Print(json);

	if (! arg->silent)
		fprintf(stdout, "\tSASA of input molecule: %s\n", arg->jsonOutFileName);

	/* print string to file */
	sprintf(outpath, "%s/%s", arg->outDirName, arg->jsonOutFileName);
	arg->jsonOutFile = safe_open(outpath, "w");
	fprintf(arg->jsonOutFile, "%s", popsOutJson);
	fclose(arg->jsonOutFile);

	free(popsOutJson);
}

/*____________________________________________________________________________*/
/* The natural way to make residues dependent on chains would be a
	double loop iterating over chains and residues. The POPS PDB data structure
	is such that it is easier here to iterate over residues
	and to create chains on the fly, to which new residue array are attached. */
void make_resSasaJson(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *json)
{
	char reslab[8];

	/* 'json' is the root object to which everything else will be attached */

	/* indices to iterate through arrays defined below */
	unsigned int c = 0; /* chain index */
	char isChainLabel[2] = "@"; /* dummy chain name, never in structure */
	unsigned int r = 0; /* residue index */

	/* header: attached to 'json' */
	cJSON_AddStringToObject(json, "data_resource", "popscomp");
	cJSON_AddStringToObject(json, "resource_version", "2.3.0");
	cJSON_AddStringToObject(json, "software_version", "2.3.0");
	cJSON_AddStringToObject(json, "resource_entry_url", "https://github.com/Fraternalilab/POPSCOMP");
	cJSON_AddStringToObject(json, "release_date", "04/11/2018");
	cJSON_AddStringToObject(json, "pdb_id", "1f3r");

	/* add Chain array */
	cJSON *chains = cJSON_AddArrayToObject(json, "chains");
	strcpy(isChainLabel, pdb->atom[pdb->resAtom[0]].chainIdentifier);

	/* iterate over all Chains */
	for (c = 0; c < pdb->nChain; ++ c) {
		cJSON *chain = cJSON_CreateObject();
		cJSON_AddItemToArray(chains, chain);
		cJSON_AddStringToObject(chain, "chain_label", isChainLabel); 
		/*cJSON_AddObjectToObject(chain, "additional_chain_annotations");*/

		/* add Residue array for new Chain */
		cJSON *residues = cJSON_AddArrayToObject(chain, "residues");

		/* iterate over all Residues */
		for ( ; r < pdb->nResidue; ++ r) {
			/* add Residue */
			cJSON *residue = cJSON_CreateObject();
			cJSON_AddItemToArray(residues, residue);
			sprintf(reslab, "%d", pdb->atom[pdb->resAtom[r]].residueNumber);
			cJSON_AddStringToObject(residue, "pdb_res_label", reslab);
			/* use original residue number of heteroresidues */
			if (pdb->atom[pdb->resAtom[r]].het) {
				cJSON_AddStringToObject(residue, "aa_type",
									pdb->atom[pdb->resAtom[r]].residueNameHet);
			} else {
				cJSON_AddStringToObject(residue, "aa_type",
									pdb->atom[pdb->resAtom[r]].residueName);
			}
			/*cJSON_AddObjectToObject(residue, "additional_residue_annotations");*/

			/* add Site_Data array */
			cJSON *site_data = cJSON_AddArrayToObject(residue, "site_data");

			/* add Sites */
			cJSON *phil = cJSON_CreateObject();
			cJSON_AddItemToArray(site_data, phil);
			cJSON_AddNumberToObject(phil, "site_id_ref", 1);
			/* truncate SASA area to 4 digits */
			/* rounded to 1 postcomma digit would suffice,
				but the format validator trips over identical entries */
			cJSON_AddNumberToObject(phil, "raw_score",
								trunc(resSasa[r].philicSasa * 1e4) / 1e4);
			cJSON_AddNumberToObject(phil, "confidence_score", 0.9);
			cJSON_AddStringToObject(phil, "confidence_classification", "high");

			cJSON *phob = cJSON_CreateObject();
			cJSON_AddItemToArray(site_data, phob);
			cJSON_AddNumberToObject(phob, "site_id_ref", 2);
			cJSON_AddNumberToObject(phob, "raw_score",
								trunc(resSasa[r].phobicSasa * 1e4) / 1e4);
			cJSON_AddNumberToObject(phob, "confidence_score", 0.9);
			cJSON_AddStringToObject(phob, "confidence_classification", "high");

			cJSON *total = cJSON_CreateObject();
			cJSON_AddItemToArray(site_data, total);
			cJSON_AddNumberToObject(total, "site_id_ref", 3);
			cJSON_AddNumberToObject(total, "raw_score",
								trunc(resSasa[r].sasa * 1e4) / 1e4);
			cJSON_AddNumberToObject(total, "confidence_score", 0.9);
			cJSON_AddStringToObject(total, "confidence_classification", "high");

			/* if the next residue has a different chain identifier */
			if (r < (pdb->nResidue - 1) &&
				strcmp(pdb->atom[pdb->resAtom[r+1]].chainIdentifier, isChainLabel) != 0) {
				/* increment Residue, update Chain identifier and break out
				   into enclosing loop for creating a new Chain,
				   including its cognate Residue array */
				r ++; 
				strcpy(isChainLabel, pdb->atom[pdb->resAtom[r]].chainIdentifier);
				break;
			}
		}
	}

	/* add Sites array */
	cJSON *sites = cJSON_AddArrayToObject(json, "sites");

	/* add Site information */
	/* hydrophilic SASA */
	cJSON *phobsite = cJSON_CreateObject();
	cJSON_AddItemToArray(sites, phobsite);
	cJSON_AddNumberToObject(phobsite, "site_id", 1);
	cJSON_AddStringToObject(phobsite, "label", "hydrophobic SASA [A^2]");
	/*
	cJSON_AddStringToObject(site, "source_database", "pdb");
	cJSON_AddStringToObject(site, "source_accession", "1f3r");
	cJSON_AddStringToObject(site, "source_release_date", "01/01/2017");
	*/
	/* hydrophilic SASA */
	cJSON *philsite = cJSON_CreateObject();
	cJSON_AddItemToArray(sites, philsite);
	cJSON_AddNumberToObject(philsite, "site_id", 2);
	cJSON_AddStringToObject(philsite, "label", "hydrophilic SASA [A^2]");
	/*
	cJSON_AddStringToObject(site, "source_database", "pdb");
	cJSON_AddStringToObject(site, "source_accession", "1f3r");
	cJSON_AddStringToObject(site, "source_release_date", "01/01/2017");
	*/
	/* total SASA */
	cJSON *totalsite = cJSON_CreateObject();
	cJSON_AddItemToArray(sites, totalsite);
	cJSON_AddNumberToObject(totalsite, "site_id", 3);
	cJSON_AddStringToObject(totalsite, "label", "total SASA [A^2]");
	/*
	cJSON_AddStringToObject(site, "source_database", "pdb");
	cJSON_AddStringToObject(site, "source_accession", "1f3r");
	cJSON_AddStringToObject(site, "source_release_date", "01/01/2017");
	*/
	/* entry annotations */
	/*cJSON_AddObjectToObject(json, "additional_entry_annotations");*/

	/* add Evidence array */
	cJSON *evidence_code_ontology = cJSON_AddArrayToObject(json, "evidence_code_ontology");
	cJSON *evidence_code_onto = cJSON_CreateObject();
	cJSON_AddItemToArray(evidence_code_ontology, evidence_code_onto);
    cJSON_AddStringToObject(evidence_code_onto, "eco_term", "computational combinatorial evidence used in automatic assertion");
    cJSON_AddStringToObject(evidence_code_onto, "eco_code", "ECO_0000246");
}

/*____________________________________________________________________________*/
void print_jsonb(Arg *arg, cJSON *jsonb)
{
	/* print JSON object to string */
	char *popsbOutJson = cJSON_Print(jsonb);

	if (! arg->silent)
		fprintf(stdout, "\tbSASA of input molecule: %s\n", arg->jsonbOutFileName);

	/* print string to file */
	arg->jsonbOutFile = safe_open(arg->jsonbOutFileName, "w");
	fprintf(arg->jsonbOutFile, "%s", popsbOutJson);
	fclose(arg->jsonbOutFile);

	free(popsbOutJson);
}

/*____________________________________________________________________________*/
/* the same as the previous routine for bSasa */
void make_resbSasaJson(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *jsonb)
{
	char reslab[8];

	/* 'jsonb' is the root object to which everything else will be attached */

	/* indices to iterate through arrays defined below */
	unsigned int c = 0; /* chain index */
	char isChainLabel[2] = "@"; /* dummy chain name, never in structure */
	unsigned int r = 0; /* residue index */

	/* header: attached to 'jsonb' */
	cJSON_AddStringToObject(jsonb, "data_resource", "popscomp");
	cJSON_AddStringToObject(jsonb, "resource_version", "3.0.0");
	cJSON_AddStringToObject(jsonb, "software_version", "3.0.0");
	cJSON_AddStringToObject(jsonb, "resource_entry_url", "https://github.com/Fraternalilab/POPSCOMP");
	cJSON_AddStringToObject(jsonb, "release_date", "21/02/2018");
	cJSON_AddStringToObject(jsonb, "pdb_id", "1f3r");

	/* add Chain array */
	cJSON *chains = cJSON_AddArrayToObject(jsonb, "chains");
	strcpy(isChainLabel, pdb->atom[pdb->resAtom[0]].chainIdentifier);

	/* iterate over all Chains */
	for (c = 0; c < pdb->nChain; ++ c) {
		cJSON *chain = cJSON_CreateObject();
		cJSON_AddItemToArray(chains, chain);
		cJSON_AddStringToObject(chain, "chain_label", isChainLabel); 
		/*cJSON_AddObjectToObject(chain, "additional_chain_annotations");*/

		/* add Residue array for new Chain */
		cJSON *residues = cJSON_AddArrayToObject(chain, "residues");

		/* iterate over all Residues */
		for ( ; r < pdb->nResidue; ++ r) {
			/* add Residue */
			cJSON *residue = cJSON_CreateObject();
			cJSON_AddItemToArray(residues, residue);
			sprintf(reslab, "%d", pdb->atom[pdb->resAtom[r]].residueNumber);
			cJSON_AddStringToObject(residue, "pdb_res_label", reslab);
			cJSON_AddStringToObject(residue, "aa_type",
									pdb->atom[pdb->resAtom[r]].residueName);
			/*cJSON_AddObjectToObject(residue, "additional_residue_annotations");*/

			/* add Site_Data array */
			cJSON *site_data = cJSON_AddArrayToObject(residue, "site_data");

			/* add Sites */
			cJSON *phil = cJSON_CreateObject();
			cJSON_AddItemToArray(site_data, phil);
			cJSON_AddNumberToObject(phil, "site_id_ref", 1);
			/* truncate bSASA area to 4 digits */
			/* rounded to 1 postcomma digit would suffice,
				but the format validator trips over identical entries */
			cJSON_AddNumberToObject(phil, "raw_score",
								trunc(resSasa[r].philicbSasa * 1e4) / 1e4);
			cJSON_AddNumberToObject(phil, "confidence_score", 0.9);
			cJSON_AddStringToObject(phil, "confidence_classification", "high");

			cJSON *phob = cJSON_CreateObject();
			cJSON_AddItemToArray(site_data, phob);
			cJSON_AddNumberToObject(phob, "site_id_ref", 2);
			cJSON_AddNumberToObject(phob, "raw_score",
								trunc(resSasa[r].phobicbSasa * 1e4) / 1e4);
			cJSON_AddNumberToObject(phob, "confidence_score", 0.9);
			cJSON_AddStringToObject(phob, "confidence_classification", "high");

			cJSON *total = cJSON_CreateObject();
			cJSON_AddItemToArray(site_data, total);
			cJSON_AddNumberToObject(total, "site_id_ref", 3);
			cJSON_AddNumberToObject(total, "raw_score",
								trunc(resSasa[r].bSasa * 1e4) / 1e4);
			cJSON_AddNumberToObject(total, "confidence_score", 0.9);
			cJSON_AddStringToObject(total, "confidence_classification", "high");

			/* if the next residue has a different chain identifier */
			if (r < (pdb->nResidue - 1) &&
				strcmp(pdb->atom[pdb->resAtom[r+1]].chainIdentifier, isChainLabel) != 0) {
				/* increment Residue, update Chain identifier and break out
				   into enclosing loop for creating a new Chain,
				   including its cognate Residue array */
				r ++; 
				strcpy(isChainLabel, pdb->atom[pdb->resAtom[r]].chainIdentifier);
				break;
			}
		}
	}

	/* add Sites array */
	cJSON *sites = cJSON_AddArrayToObject(jsonb, "sites");
	cJSON *site = cJSON_CreateObject();
	cJSON_AddItemToArray(sites, site);

	/* add Site information */
	/* hydrophilic bSASA */
	cJSON_AddNumberToObject(site, "site_id", 1);
	cJSON_AddStringToObject(site, "label", "hydrophobic bSASA [A^2]");
	cJSON_AddStringToObject(site, "source_database", "pdb");
	cJSON_AddStringToObject(site, "source_accession", "1f3r");
	cJSON_AddStringToObject(site, "source_release_date", "01/01/2017");

	/* hydrophilic bSASA */
	cJSON_AddNumberToObject(site, "site_id", 2);
	cJSON_AddStringToObject(site, "label", "hydrophilic bSASA [A^2]");
	cJSON_AddStringToObject(site, "source_database", "pdb");
	cJSON_AddStringToObject(site, "source_accession", "1f3r");
	cJSON_AddStringToObject(site, "source_release_date", "01/01/2017");

	/* total bSASA */
	cJSON_AddNumberToObject(site, "site_id", 3);
	cJSON_AddStringToObject(site, "label", "total bSASA [A^2]");
	cJSON_AddStringToObject(site, "source_database", "pdb");
	cJSON_AddStringToObject(site, "source_accession", "1f3r");
	cJSON_AddStringToObject(site, "source_release_date", "01/01/2017");

	/* entry annotations */
	/*cJSON_AddObjectToObject(jsonb, "additional_entry_annotations");*/

	/* add Evidence array */
	cJSON *evidence_code_ontology = cJSON_AddArrayToObject(jsonb, "evidence_code_ontology");
	cJSON *evidence_code_onto = cJSON_CreateObject();
	cJSON_AddItemToArray(evidence_code_ontology, evidence_code_onto);
    cJSON_AddStringToObject(evidence_code_onto, "eco_term", "computational combinatorial evidence used in automatic assertion");
    cJSON_AddStringToObject(evidence_code_onto, "eco_code", "ECO_0000246");
}

