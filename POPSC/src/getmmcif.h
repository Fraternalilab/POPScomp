/*===============================================================================
getmmcif.h : read MMCIF structures 

Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETMMCIF_H
#define GETMMCIF_H

#include <limits.h>

#include "arg.h"
#include "argpdb.h"
#include "getpdb.h"
#include "pdb_structure.h" 
#include "safe.h"

/*____________________________________________________________________________*/
/* prototypes */
int map_structure_mmcif(Arg *arg, Argpdb *argpdb, Str *str, Structure *s);
/* the function header for 'read_cif()' is defined
      in the C++ wrapper 'cif_header.h' */

#endif

