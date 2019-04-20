/*==============================================================================
POPS* : Parameter OPtimised Surface of proteins and nucleic acids
Copyright (C) 2002-2018 Franca Fraternali (program author, parametrisation)
Copyright (C) 2008-2018 Jens Kleinjung (modular C code)
Copyright (C) 2002 Luigi Cavallo (parametrisation)
Copyright (C) 2002 Kuang Lin and Valerie Hindie (translation to C)

POPS is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Read the COPYING file for license and distribution details.

POPS : Parameter OPtimised Surfaces at atomic resolution.
		 F. Fraternali and L. Cavallo (2002)
		 Nucleic Acids Res. 30(13), 2950-2960.
==============================================================================*/

#if !defined POPS_H
#define POPS_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "argpdb.h"
#include "atomgroup.h"
#include "atomgroup_data.h"
#include "cJSON.h"
#include "fs.h"
#include "getpdb.h"
#include "getpdbml.h"
#include "gettraj.h"
#include "json.h"
#include "matrix.h"
#include "putsasa.h"
#include "putbsasa.h"
#include "putsigma.h"
#include "safe.h"
#include "sasa.h"
#include "sasa_const.h"
#include "sasa_data.h"
#include "sigma.h"
#include "sigma_const.h"
#include "sigma_data.h"
#include "topol.h"
#include "type.h"

#endif
