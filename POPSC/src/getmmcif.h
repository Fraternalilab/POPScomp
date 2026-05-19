/*===============================================================================
getmmcif.h : read MMCIF structures 

Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETMMCIF_H
#define GETMMCIF_H

#include "arg.h"
#include "argpdb.h"
#include "pdb_structure.h" 

/*____________________________________________________________________________*/
/* prototypes */
int map_structure_mmcif(Arg *arg, Argpdb *argpdb, Str *pdb, Structure *s);
/* the function header for 'read_cif()' is defined
      in the C++ wrapper 'cif_header.h' */

#endif

