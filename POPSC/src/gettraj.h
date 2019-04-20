/*===============================================================================
getTRAJ.h : Read GROMOS96 trajectory file
Copyright (C) 2008-2011 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETTRAJ_H
#define GETTRAJ_H

#include <math.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getpdb.h"
#include "pdb_structure.h"
#include "safe.h"
#include "vector.h"

/*____________________________________________________________________________*/
/* structures */

/** trajectory atom */
typedef struct
{
	Vec pos; /* trajectory atom position */
} Trajatom;

/** trajectory frame */
typedef struct
{
	int nAtom; /* number of atoms (without ions) */
	int nIon; /* number of ions */
	int nSolv; /* number of solvent molecules */
	Trajatom *trajatom; /* trajectory atom */
} Frame;

/** trajectory */
typedef struct
{
	int nFrame; /* number of atoms (without ions) */
	Frame *frame; /* xyz coordinate vector of atom */
} Traj;


/*____________________________________________________________________________*/
/** prototypes */
void copy_coordinates(Str *pdb, Traj *traj, int frame);
int read_gromos_traj(Traj *traj, Arg *arg, int protEnd);

#endif
