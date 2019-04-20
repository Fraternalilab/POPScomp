/*=============================================================================
matrix.h : allocate matrices (with dicontinuous memory allocation)
Copyright (C) Jens Kleinjung 2007
Read the COPYING file for license information.
=============================================================================*/

#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>

#include "safe.h"
#include "vector.h"

/*___________________________________________________________________________*/
/* prototypes */
/* integer */
int **alloc_mat2D_int(int **mat2D_int, int x, int y);
void init_mat2D_int(int **mat2D_int, int x, int y, int val);
void free_mat2D_int(int **mat2D_int, int x);
void print_mat2D_int(char *outFileName, int **mat2D_int, int x, int y);

int ***alloc_mat3D_int(int ***mat3D_int, int x, int y, int z);
void init_mat3D_int(int ***mat3D_int, int x, int y, int z, int val);
void free_mat3D_int(int ***mat3D_int, int x, int y);

int ****alloc_mat4D_int(int ****mat4D_int, int w, int x, int y, int z);
void init_mat4D_int(int ****mat4D_int, int w, int x, int y, int z, int val);
void free_mat4D_int(int ****mat4D_int, int w, int x, int y);

/* float */
float **alloc_mat2D_float(float **mat2D_float, int x, int y);
void init_mat2D_float(float **mat2D_float, int x, int y, float val);
void free_mat2D_float(float **mat2D_float, int x);
void print_mat2D_float(char *outFileName, float **mat2D_float, int x, int y);
void div_mat2D_float(float **mat2D_float, int x, int y, float a);

float ***alloc_mat3D_float(float ***mat3D_float, int x, int y, int z);
void init_mat3D_float(float ***mat3D_float, int x, int y, int z, float val);
void free_mat3D_float(float ***mat3D_float, int x, int y);

float ****alloc_mat4D_float(float ****mat4D_float, int w, int x, int y, int z);
void init_mat4D_float(float ****mat4D_float, int w, int x, int y, int z, float val);
void free_mat4D_float(float ****mat4D_float, int w, int x, int y);

/* vector */
Vec **alloc_mat2D_vec(Vec **mat2D_vec, int x, int y);
void init_mat2D_vec(Vec **mat2D_vec, int x, int y, Vec val);
void free_mat2D_vec(Vec **mat2D_vec, int x);

Vec ***alloc_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z);
void init_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z, Vec val);
void free_mat3D_vec(Vec ***mat3D_vec, int x, int y);

Vec ****alloc_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z);
void init__mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z, Vec val);
void free_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y);

#endif

