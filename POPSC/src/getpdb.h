/*===============================================================================
getpdb.h : read PDB structures 
Copyright (C) 2004 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETPDB_H
#define GETPDB_H

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
/* prototypes */
int read_pdb(FILE *pdbfile, gzFile *pdbgzInFile, Arg *arg, Argpdb *argpdb, Str *str);
void read_structure(Arg *arg, Argpdb *argpdb, Str *pdb);

#endif
