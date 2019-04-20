/*==============================================================================
argpdb.h : parse command line arguments for PDB processing
Copyright (C) 2009 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARGPDB_H
#define ARGPDB_H

/*____________________________________________________________________________*/
/* structures */
typedef struct
{ 
	int coarse; /* Calpha and Pphosphate atoms only */
	int hydrogens; /* read hydrogen atoms */
	int multiModel; /* read multiple models */
	int partOcc; /* partial occupancy */
} Argpdb;

#endif
