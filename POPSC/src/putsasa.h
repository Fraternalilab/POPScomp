/*===============================================================================
putsasa.h : routines for printing SASA
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef PUTSASA_H
#define PUTSASA_H

#include <stdio.h>

#include "arg.h"
#include "argpdb.h"
#include "getpdb.h"
#include "math.h"
#include "sasa.h"
#include "sasa_const.h"
#include "topol.h"
#include "type.h"

void print_sasa(Arg *arg, Argpdb *argpdb, Str *pdb, Type *type, Topol *topol, \
				MolSasa *molSasa, ConstantSasa *constant_sasa, int frame);

#endif
