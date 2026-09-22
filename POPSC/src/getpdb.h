/*===============================================================================
getpdb.h : read PDB structures 
Copyright (C) 2004 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETPDB_H
#define GETPDB_H

#include <ctype.h>
#include <libgen.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>

#include "arg.h"
#include "argpdb.h"
#include "error.h"
#include "modstring.h"
#include "pdb_structure.h"
#include "pattern.h"
#include "safe.h"
#include "seq.h"
#include "vector.h"

/*____________________________________________________________________________*/
/* index of the unknown polymer residue 'UNK' in the parameter tables
	('sasa_data.h', 'sigma_data.h', 'atomgroup_data.h') */
#define UNK_RESIDUE_TYPE 36
/* index of the generic terminal-atom residue 'ANY' (OXT, O1, O2) */
#define ANY_RESIDUE_TYPE 20

/*____________________________________________________________________________*/
/* prototypes */
void remove_spaces(char *s);
char aacode(char *code3);
int is_water(const char *residueName);
int is_coarse_atom(const char *residueName, const char *atomName);
void index_structure(Str *str);
void set_pdbID_from_filename(Str *str, const char *fileName);
void exit_structure_error(Arg *arg, Str *str, const char *message);
const char *unk_atom_name(Atom *atom);
const char *chain_label(Atom *atom);
int read_pdb(FILE *pdbfile, gzFile *pdbgzInFile, Arg *arg, Argpdb *argpdb, Str *str);
void read_structure(Arg *arg, Argpdb *argpdb, Str *pdb);

#endif
