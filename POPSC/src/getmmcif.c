/*==============================================================================
getmmcif.c : routines for reading MMCIF structures
Copyright (C) 2004 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "cif_reader.h"
#include "getmmcif.h"


/*____________________________________________________________________________*/
/* map MMCIF structure */
int map_structure_mmcif(Arg *arg, Argpdb *argpdb, Str *pdb, Structure *s) {
    Structure *struc = NULL;

    struc = read_cif("test.cif.gz");

    if(struc == NULL) {
        fprintf(stderr, "Failed to read structure\n");
        return 1;
    }

    printf("Number of atoms: %d\n", struc->natom);

    /* access coordinates */
    printf("%f %f %f\n",
           struc->xyz[0],
           struc->xyz[1],
           struc->xyz[2]);

    free_structure(struc);

    return 0;
}

