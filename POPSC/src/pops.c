/*=============================================================================
POPS* : Parameter OPtimised Surface of proteins and nucleic acids
Copyright (C) 2002-2018 Franca Fraternali (program author, parametrisation)
Copyright (C) 2008-2018 Jens Kleinjung (modular C code)
Copyright (C) 2002 Luigi Cavallo (parametrisation)
Copyright (C) 2002 Kuang Lin and Valerie Hindie (translation to C)

POPS is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Read the COPYING file for license and distribution details.

Fraternali, F. and van Gunsteren, W.F.
An efficient mean solvation force model for use in molecular dynamics simulations
  of proteins in aqueous solution.
  Journal of Molecular Biology 256 (1996) 939-948.

Fraternali, F. and Cavallo, L.
Parameter optimized surfaces (POPS*): analysis of key interactions and 
  conformational changes in the ribosome.
  Nucleic Acids Research 30 (2002) 2950-2960.

Cavallo, L., Kleinjung, J. and Fraternali, F.
POPS: A fast algorithm for solvent accessible surface areas at atomic 
  and residue level.
  Nucleic Acids Research 31 (2003) 3364-3366.

Kleinjung, J. and Fraternali, F.
  POPSCOMP: an automated interaction analysis of biomolecular complexes.
  Nucleic Acids Research 33 (2005) W342-W346.

=============================================================================*/

#include "config.h"
#include "pops.h"

/*____________________________________________________________________________*/
/* global variables */
#ifdef MPI
#include <mpi.h>
    int my_rank; /* rank of 'this' node */
    int nodes; /* number of nodes */
#else
    int my_rank = 0; 
    int nodes = 1; 
#endif

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	int i;

    /*____________________________________________________________________________*/
	/* data structures and variables */
    Arg arg; /** data structure for command line arguments */
    Argpdb argpdb; /** data structure for PDB command line arguments */
	Str pdb; /** data structure for input PDB molecule */
	Traj traj; /** data structure for trajectory */
	MolSasa molSasa; /* data structure for molecular SASA values,
					invoking structures for residuic and atomic SASA values */
	MolSFE molSFE; /* data structure for molecular SFE values,
					invoking structures for residuic and atomic SFE values */
	Topol topol; /* molecular topology */
	Type type; /* atom and residue types */
	ConstantSasa *constant_sasa; /* selected SASA constants */ 
	ConstantSasa *res_sasa; /* residue SASA constants */ 
	ConstantSigma *constant_sigma; /* SIGMA constants */ 
	Atomgroup *atomGroup; /* atom group constants */
	extern int status; /* program status from 'error' library */

	/* create JSON object */
	cJSON *resSasaJson = cJSON_CreateObject();
    if (resSasaJson == NULL) {
		cJSON_Delete(resSasaJson);
		ErrorSpec("Exiting", "JSON object returned NULL");
    }

	/* create JSONb object */
	cJSON *resSasaJsonb = cJSON_CreateObject();
    if (resSasaJsonb == NULL) {
		cJSON_Delete(resSasaJsonb);
		ErrorSpec("Exiting", "JSONb object returned NULL");
    }

    /*________________________________________________________________________*/
    /* MPI */
#ifdef MPI
    /* initialize MPI routines */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    /*____________________________________________________________________________*/
    /** parse command line arguments */
    parse_args(argc, &(argv[0]), &arg, &argpdb);

    /*____________________________________________________________________________*/
	/** select SASA parameters */
	if (! argpdb.coarse) {
		constant_sasa = &(constant_sasa_data[0]); /* data array element 0: atomic */
	} else {
		constant_sasa = &(constant_sasa_data[1]); /* data array element 1: coarse-grained */
	}
	res_sasa = &(constant_sasa_data[1]);

    /*____________________________________________________________________________*/
	/** select SIGMA parameters */
	if (! argpdb.coarse) {
		constant_sigma = &(constant_sigma_data[0]); /* data array element 0: atomic */
	} else {
		constant_sigma = &(constant_sigma_data[1]); /* data array element 1: coarse-grained */
	}

    /*____________________________________________________________________________*/
    /** read input structure, either XML format or classic PDB format */
	if (! arg.silent) fprintf(stdout, "Input structure\n");
	if (arg.pdbml) {
		read_structure_xml(&arg, &argpdb, &pdb);
	} else {
		read_structure(&arg, &argpdb, &pdb);
	}

    /*____________________________________________________________________________*/
    /** read input trajectory */
	if (arg.trajInFileName) {
		if (! arg.silent) fprintf(stdout, "Input trajectory\n");
		read_gromos_traj(&traj, &arg, pdb.nAllAtom);
	}

    /*___________________________________________________________________________*/
    /* assign atom group ID number */
    atomGroup = &(constAtomGroup[0]); /* group array element 0: POPS grouping */
	if (! arg.silent) fprintf(stdout, "Groups\n");
    get_atomgroup(&pdb, atomGroup); /* assign atom grouping */

	/*____________________________________________________________________________*/
    /** define atom and residue types */
	if (! arg.silent) fprintf(stdout, "Types\n");
    get_types(&pdb, &type, constant_sasa);

    /*____________________________________________________________________________*/
    /** compute molecular topology */
	if (! arg.silent) fprintf(stdout, "Topology\n");
	init_topology(&pdb, &topol);
    get_topology(&pdb, &type, &topol, constant_sasa, &argpdb, &arg);

    /*____________________________________________________________________________*/
    /** compute SASA: atoms, residues, chains, molecule */
	if (! arg.silent) fprintf(stdout, "Solvent Accessible Surface Area\n");
    init_sasa(&pdb, &type, &molSasa, constant_sasa, &arg);
    compute_sasa(&pdb, &topol, &type, &molSasa, constant_sasa, res_sasa, &arg);
    
    /*____________________________________________________________________________*/
	/* SASA output, bSASA is buried area */
	if (arg.jsonOut) {
		/** print JSON output */
		if (! arg.silent) fprintf(stdout, "SASA Output:\n");
		make_resSasaJson(&arg, &pdb, molSasa.resSasa, resSasaJson);
		print_json(&arg, resSasaJson);
		if (! arg.silent) fprintf(stdout, "bSASA Output:\n");
		/* disabled, because not validated against PDBe server */
		/* make_resbSasaJson(&arg, &pdb, molSasa.resSasa, resSasaJsonb);
		   print_jsonb(&arg, resSasaJsonb); */
	} else {
		/** print tabulated output */
		if (! arg.silent) fprintf(stdout, "SASA Output:\n");
		print_sasa(&arg, &argpdb, &pdb, &type, &topol, &molSasa, constant_sasa, -1);
		if (! arg.silent) fprintf(stdout, "bSASA Output:\n");
		print_bsasa(&arg, &argpdb, &pdb, &type, &topol, &molSasa, constant_sasa, -1);
	}

    /*____________________________________________________________________________*/
    /** compute Solvation Free Energy: atoms, residues, chains, molecule */
	if (! arg.silent) fprintf(stdout, "Solvation Free Energy\n");
    init_sfe(&pdb, &type, &molSFE, constant_sigma, &arg);
    compute_sfe(&pdb, &type, &molSasa, &molSFE, constant_sigma, &arg);
    
    /*____________________________________________________________________________*/
	/** print Solvation Free Energy */
	if (! arg.silent) fprintf(stdout, "SFE Output:\n");
	/* we don't have SFEs for residues yet */
	if (! argpdb.coarse && ! arg.jsonOut)
		print_sfe(&arg, &argpdb, &pdb, &type, &topol, &molSFE, constant_sigma, -1);

    /*____________________________________________________________________________*/
	/** free memory */
	free_topology(&pdb, &topol);
	free_sasa(&molSasa);
	free_sfe(&molSFE);

    /*____________________________________________________________________________*/
	/** GROMOS trajectory */
	if (arg.trajInFileName) {
		if (! arg.silent)
			fprintf(stdout, "SASA Output of trajectory frames: %s.*.out\n\t",
				arg.sasatrajOutFileName);

		for (i = 0; i < traj.nFrame; ++ i) {
			if (! arg.silent) {
				(((i+1) % 50) != 0) ? fprintf(stdout, ".") : fprintf(stdout, "%d\n\t", (i + 1));
				fflush(stdout);
			}
			assert(traj.frame[i].nAtom == pdb.nAllAtom);
			copy_coordinates(&pdb, &traj, i);
			/* topology */
			init_topology(&pdb, &topol);
			get_topology(&pdb, &type, &topol, constant_sasa, &argpdb, &arg);
			/* SASA */
			init_sasa(&pdb, &type, &molSasa, constant_sasa, &arg);
			compute_sasa(&pdb, &topol, &type, &molSasa, constant_sasa, res_sasa, &arg);
			print_sasa(&arg, &argpdb, &pdb, &type, &topol, &molSasa, constant_sasa, i);
			/* SFE */
			init_sfe(&pdb, &type, &molSFE, constant_sigma, &arg);
			compute_sfe(&pdb, &type, &molSasa, &molSFE, constant_sigma, &arg);
			/*print_sfe(&arg, &argpdb, &pdb, &type, &topol, &molSFE, constant_sigma, i);*/
			/* free memory */
			free_topology(&pdb, &topol);
			free_sasa(&molSasa);
			free_sfe(&molSFE);
		}
	}

    /*____________________________________________________________________________*/
	/** free memory */
	/* structure */
	free(pdb.sequence.name);
	free(pdb.atom);
	free(pdb.resAtom);
	free(pdb.atomMap);
	free(pdb.sequence.res);

	/* trajectory */
	if (arg.trajInFileName) {
		for (i = 0; i < (traj.nFrame + 1); ++ i)
			free(traj.frame[i].trajatom);
		free(traj.frame);
	}

	/* type */
	free(type.atomType);
	free(type.residueType);

	/* JSON object */
	cJSON_Delete(resSasaJson);
	cJSON_Delete(resSasaJsonb);

    /*________________________________________________________________________*/
    /* MPI */
#ifdef MPI
    /* stop MPI processes */
    MPI_Finalize();
#endif

    /*____________________________________________________________________________*/
	/* terminate */
	if (! arg.silent) fprintf(stdout, "\nClean termination\n\n");

	status = 0;

    return status;
}
