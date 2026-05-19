/*==============================================================================
cif_reader.h : wrapper header to use C++ functions from 'gemmi' 
The 'gemmi' library provides functions for the now default mmcif format.
Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif

/*___________________________________________________________________________*/

typedef struct {
    int natom;

    double* xyz;

    char** atom_name;
    int* atom_number;

    char** res_name;
    int* res_number;
	char* ins_code;

    char** chain_name;
    int chain_number;
} Structure;

/*___________________________________________________________________________*/
Structure* read_cif(const char* filename);
void free_structure(Structure* s);

/*___________________________________________________________________________*/
#ifdef __cplusplus
}
#endif

