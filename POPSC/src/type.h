/*===============================================================================
type.h : determine atom and residue types
Copyright (C) 2002-2018 Franca Fraternali
Copyright (C) 2008-2018 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin
Read the COPYING file for license information.
================================================================================*/

#ifndef TYPE_H
#define TYPE_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getpdb.h"
#include "modstring.h"
#include "sasa_const.h"

/*____________________________________________________________________________*/
/* structures */
typedef struct {
	int *atomType; /* identifier to denote atom type */
	int *residueType; /* identifier to denote residue type */
} Type;

/*____________________________________________________________________________*/
/* prototypes */
int get_types(Str *pdb, Type *type, ConstantSasa *constant_sasa);

#endif

