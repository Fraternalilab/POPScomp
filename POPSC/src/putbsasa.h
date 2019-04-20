/*===============================================================================
putbsasa.h : routines for printing bSASA
Copyright (C) 2016-2018 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef PUTBSASA_H
#define PUTBSASA_H

#include <stdio.h>

#include "arg.h"
#include "argpdb.h"
#include "getpdb.h"
#include "math.h"
#include "sasa.h"
#include "sasa_const.h"
#include "topol.h"
#include "type.h"

void print_bsasa(Arg *arg, Argpdb *argpdb, Str *pdb, Type *type, Topol *topol, \
				MolSasa *molSasa, ConstantSasa *constant_sasa, int frame);

#endif

