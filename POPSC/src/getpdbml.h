/*===============================================================================
getpdbml.h : read PDBML structures 
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETPDBML_H
#define GETPDBML_H

#include <libgen.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

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
void read_structure_xml(Arg *arg, Argpdb *argpdb, Str *pdb);

#endif
