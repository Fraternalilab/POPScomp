/*==============================================================================
json.c : JSON routines, using the cJSON library
To validate JSON output format:
../funpdbe-client/funpdbe_client.py --path=pops.json --mode=validate
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "json.h"
#include "cJSON.h"

/*____________________________________________________________________________*/
/** write a JSON object to '<outDir>/<pdbID><suffix>' */
static void write_json(Arg *arg, Str *pdb, cJSON *json, const char *suffix)
{
	char name[512];
	char *popsOutJson = cJSON_Print(json);
	FILE *jsonFile;

	snprintf(name, sizeof(name), "%s%s", pdb->pdbID, suffix);

	if (! arg->silent)
		fprintf(stdout, "\tSASA of input molecule: %s\n", name);

	jsonFile = open_output(arg, name, "w");
	fprintf(jsonFile, "%s", popsOutJson);
	fclose(jsonFile);

	free(popsOutJson);
}

/*____________________________________________________________________________*/
void print_json(Arg *arg, Str *pdb, cJSON *json)
{
	write_json(arg, pdb, json, ".json");
}

/*____________________________________________________________________________*/
void print_jsonb(Arg *arg, Str *pdb, cJSON *jsonb)
{
	write_json(arg, pdb, jsonb, ".b.json");
}

/*____________________________________________________________________________*/
/** add one site value to the site data of a residue */
static void add_site_data(cJSON *site_data, int site_id, double value)
{
	cJSON *site = cJSON_CreateObject();
	cJSON_AddItemToArray(site_data, site);
	cJSON_AddNumberToObject(site, "site_id_ref", site_id);
	/* truncate SASA area to 4 digits */
	/* rounded to 1 postcomma digit would suffice,
		but the format validator trips over identical entries */
	cJSON_AddNumberToObject(site, "raw_score", trunc(value * 1e4) / 1e4);
	cJSON_AddNumberToObject(site, "confidence_score", 0.9);
	cJSON_AddStringToObject(site, "confidence_classification", "high");
}

/*____________________________________________________________________________*/
/** add one site definition */
static void add_site(cJSON *sites, int site_id, const char *label)
{
	cJSON *site = cJSON_CreateObject();
	cJSON_AddItemToArray(sites, site);
	cJSON_AddNumberToObject(site, "site_id", site_id);
	cJSON_AddStringToObject(site, "label", label);
}

/*____________________________________________________________________________*/
/* The natural way to make residues dependent on chains would be a
	double loop iterating over chains and residues. The POPS PDB data structure
	is such that it is easier here to iterate over residues
	and to create chains on the fly, to which new residue array are attached.
	Residues are the SASA residues ('resSasa', one per residue index),
	so that residue values and labels always refer to the same residue. */
static void make_json(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *json, int buried)
{
	char reslab[32];
	unsigned int r = 0; /* residue index */
	int chainIndex = -1; /* chain index of the current chain object */
	Atom *ref; /* reference atom of the residue */
	cJSON *chains = NULL;
	cJSON *residues = NULL;

	/* header: attached to 'json' */
	cJSON_AddStringToObject(json, "data_resource", "POPScomp_PDBML");
	cJSON_AddStringToObject(json, "resource_version", "weekly");
	/* VERSION is set by the build system (configure.ac), so that the version
		of the JSON output cannot drift from the version of the program */
	cJSON_AddStringToObject(json, "software_version", VERSION);
	cJSON_AddStringToObject(json, "resource_entry_url", "https://github.com/Fraternalilab/POPScomp");
	cJSON_AddStringToObject(json, "release_date", "25/10/2021");
	cJSON_AddStringToObject(json, "pdb_id", pdb->pdbID);

	if (pdb->nAllResidue == 0 || pdb->atom == NULL) {
    	fprintf(stderr, "No residues/atoms available for JSON output\n");
    	return;
	}

	/* add Chain array */
	chains = cJSON_AddArrayToObject(json, "chains");

	/* iterate over all Residues, opening a new Chain whenever the chain index changes */
	for (r = 0; r < pdb->nAllResidue; ++ r) {
		ref = &(pdb->atom[resSasa[r].atomRef]);

		if (ref->chainIndex != chainIndex) {
			cJSON *chain = cJSON_CreateObject();
			cJSON_AddItemToArray(chains, chain);
			cJSON_AddStringToObject(chain, "chain_label", ref->chainIdentifier); 
			residues = cJSON_AddArrayToObject(chain, "residues");
			chainIndex = ref->chainIndex;
		}

		/* add Residue: residue number and insertion code */
		cJSON *residue = cJSON_CreateObject();
		cJSON_AddItemToArray(residues, residue);
		if (ref->icode[0] == '-')
			snprintf(reslab, sizeof(reslab), "%d", ref->residueNumber);
		else
			snprintf(reslab, sizeof(reslab), "%d%s", ref->residueNumber, ref->icode);
		cJSON_AddStringToObject(residue, "pdb_res_label", reslab);
		cJSON_AddStringToObject(residue, "aa_type", ref->residueName);

		/* add Site_Data array: 1 hydrophilic, 2 hydrophobic, 3 total */
		cJSON *site_data = cJSON_AddArrayToObject(residue, "site_data");
		if (! buried) {
			add_site_data(site_data, 1, resSasa[r].philicSasa);
			add_site_data(site_data, 2, resSasa[r].phobicSasa);
			add_site_data(site_data, 3, resSasa[r].sasa);
		} else {
			add_site_data(site_data, 1, resSasa[r].philicbSasa);
			add_site_data(site_data, 2, resSasa[r].phobicbSasa);
			add_site_data(site_data, 3, resSasa[r].bSasa);
		}
	}

	/* add Sites array: labels match the site data above */
	cJSON *sites = cJSON_AddArrayToObject(json, "sites");
	if (! buried) {
		add_site(sites, 1, "hydrophilic SASA [A^2]");
		add_site(sites, 2, "hydrophobic SASA [A^2]");
		add_site(sites, 3, "total SASA [A^2]");
	} else {
		add_site(sites, 1, "hydrophilic bSASA [A^2]");
		add_site(sites, 2, "hydrophobic bSASA [A^2]");
		add_site(sites, 3, "total bSASA [A^2]");
	}

	/* add Evidence array */
	cJSON *evidence_code_ontology = cJSON_AddArrayToObject(json, "evidence_code_ontology");
	cJSON *evidence_code_onto = cJSON_CreateObject();
	cJSON_AddItemToArray(evidence_code_ontology, evidence_code_onto);
    cJSON_AddStringToObject(evidence_code_onto, "eco_term", "computational combinatorial evidence used in automatic assertion");
    cJSON_AddStringToObject(evidence_code_onto, "eco_code", "ECO_0000246");
}

/*____________________________________________________________________________*/
void make_resSasaJson(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *json)
{
	make_json(arg, pdb, resSasa, json, 0);
}

/*____________________________________________________________________________*/
/* the same as the previous routine for bSasa */
void make_resbSasaJson(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *jsonb)
{
	make_json(arg, pdb, resSasa, jsonb, 1);
}
