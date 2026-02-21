/*==============================================================================
cif_header.h : wrapper header to use C++ functions from 'gemma' 
Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif

/*___________________________________________________________________________*/
typedef struct {
    int natom;
    double *xyz;   // length 3*natom
} Structure;

/*___________________________________________________________________________*/
Structure* read_cif(const char* filename);
void free_structure(Structure* s);

/*___________________________________________________________________________*/
#ifdef __cplusplus
}
#endif

