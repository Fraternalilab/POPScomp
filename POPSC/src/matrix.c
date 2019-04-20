/*=============================================================================
matrix.c : allocate matrices with dicontinuous memory allocation
Copyright (C) Jens Kleinjung 2007
Read the COPYING file for license information.
==============================================================================*/

#include "matrix.h"

/*___________________________________________________________________________*/
/** 2D integer matrix */
/** allocate */
int **alloc_mat2D_int(int **mat2D_int, int x, int y)
{
    unsigned int i;

	mat2D_int = (int **)safe_malloc(x * sizeof(int *));
	for (i = 0; i < x; ++ i)
		mat2D_int[i] = (int *)safe_malloc(y * sizeof(int));

	return mat2D_int;
}

/** initialise */
void init_mat2D_int(int **mat2D_int, int x, int y, int val)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_int[i][j] = val;
}

/** free */
void free_mat2D_int(int **mat2D_int, int x)
{
    unsigned int i;

	for (i = 0; i < x; ++ i)
		free(mat2D_int[i]);

    free(mat2D_int);
}

/** print */
void print_mat2D_int(char *outFileName, int **mat2D_int, int x, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = 0; i < x; ++ i) {
		if (i > 0) fprintf(outFile, "\n");
		for (j = 0; j < y; ++ j) {
			fprintf(outFile, "%2d ", mat2D_int[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/*___________________________________________________________________________*/
/** 3D integer matrix */
/** allocate */
int ***alloc_mat3D_int(int ***mat3D_int, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_int = (int ***)safe_malloc(x * sizeof(int **));
	for (i = 0; i < x; ++ i) {
		mat3D_int[i] = (int **)safe_malloc(y * sizeof(int *));
		for (j = 0; j < y; ++ j)
			mat3D_int[i][j] = (int *)safe_malloc(z * sizeof(int));
	}

	return mat3D_int;
}

/** initialise */
void init_mat3D_int(int ***mat3D_int, int x, int y, int z, int val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				mat3D_int[i][j][k] = val;
}

/** free */
void free_mat3D_int(int ***mat3D_int, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_int[i][j]);
        free(mat3D_int[i]);
    }
    free(mat3D_int);
}

/*___________________________________________________________________________*/
/** 4D integer matrix */
/** allocate */
int ****alloc_mat4D_int(int ****mat4D_int, int w, int x, int y, int z)
{
	unsigned int h, i, j;

	mat4D_int = (int ****)safe_malloc(w * sizeof(int ***));
	for (h = 0; h < w; ++ h) {
		mat4D_int[h] = (int ***)safe_malloc(x * sizeof(int **));
		for (i = 0; i < x; ++ i) {
			mat4D_int[h][i] = (int **)safe_malloc(y * sizeof(int *));
			for (j = 0; j < y; ++ j)
				mat4D_int[h][i][j] = (int *)safe_malloc(z * sizeof(int));
		}
	}

	return mat4D_int;
}

/** initialise */
void init_mat4D_int(int ****mat4D_int, int w, int x, int y, int z, int val)
{
    unsigned int h, i, j, k;

	for (h = 0; h < w; ++ h)
		for (i = 0; i < x; ++ i)
			for (j = 0; j < y; ++ j)
				for (k = 0; k < z; ++ k)
					mat4D_int[h][i][j][k] = val;
}

/** free */
void free_mat4D_int(int ****mat4D_int, int w, int x, int y)
{
    unsigned int h, i, j;

    for (h = 0; h < w; ++ h) {
		for (i = 0; i < x; ++ i) {
			for (j = 0; j < y; ++ j)
				free(mat4D_int[h][i][j]);
			free(mat4D_int[h][i]);
	    }
		free(mat4D_int[h]);
	}
	free(mat4D_int);
}

/*___________________________________________________________________________*/
/** 2D float matrix */
/** allocate */
float **alloc_mat2D_float(float **mat2D_float, int x, int y)
{
    unsigned int i;

	mat2D_float = (float **)safe_malloc(x * sizeof(float *));
	for (i = 0; i < x; ++ i)
		mat2D_float[i] = (float *)safe_malloc(y * sizeof(float));

	return mat2D_float;
}

/** initialise */
void init_mat2D_float(float **mat2D_float, int x, int y, float val)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_float[i][j] = val;
}

/* free */
void free_mat2D_float(float **mat2D_float, int x)
{
    unsigned int i;

	for (i = 0; i < x; ++ i)
		free(mat2D_float[i]);
    free(mat2D_float);
}

/** print */
void print_mat2D_float(char *outFileName, float **mat2D_float, int x, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = 0; i < x; ++ i) {
		if (i > 0) fprintf(outFile, "\n");
		for (j = 0; j < y; ++ j) {
			fprintf(outFile, "%3.2f ", mat2D_float[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** divide by factor */
void div_mat2D_float(float **mat2D_float, int x, int y, float a)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i) {
		for (j = 0; j < y; ++ j) {
			mat2D_float[i][j] /= a;
		}
	}
}

/*___________________________________________________________________________*/
/** 3D float matrix */
/** allocate */
float ***alloc_mat3D_float(float ***mat3D_float, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_float = (float ***)safe_malloc(x * sizeof(float **));
	for (i = 0; i < x; ++ i) {
		mat3D_float[i] = (float **)safe_malloc(y * sizeof(float *));
		for (j = 0; j < y; ++ j)
			mat3D_float[i][j] = (float *)safe_malloc(z * sizeof(float));
	}

	return mat3D_float;
}

/** initialise */
void init_mat3D_float(float ***mat3D_float, int x, int y, int z, float val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				mat3D_float[i][j][k] = val;
}


/** free */
void free_mat3D_float(float ***mat3D_float, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_float[i][j]);
        free(mat3D_float[i]);
    }
    free(mat3D_float);
}

/*___________________________________________________________________________*/
/** 4D float matrix */
/** allocate */
float ****alloc_mat4D_float(float ****mat4D_float, int w, int x, int y, int z)
{
	unsigned int h, i, j;

	mat4D_float = (float ****)safe_malloc(w * sizeof(float ***));
	for (h = 0; h < w; ++ h) {
		mat4D_float[h] = (float ***)safe_malloc(x * sizeof(float **));
		for (i = 0; i < x; ++ i) {
			mat4D_float[h][i] = (float **)safe_malloc(y * sizeof(float *));
			for (j = 0; j < y; ++ j)
				mat4D_float[h][i][j] = (float *)safe_malloc(z * sizeof(float));
		}
	}

	return mat4D_float;
}

/** initialise */
void init_mat4D_float(float ****mat4D_float, int w, int x, int y, int z, float val)
{
    unsigned int h, i, j, k;

	for (h = 0; h < w; ++ h)
		for (i = 0; i < x; ++ i)
			for (j = 0; j < y; ++ j)
				for (k = 0; k < z; ++ k)
					mat4D_float[h][i][j][k] = val;
}


/** free */
void free_mat4D_float(float ****mat4D_float, int w, int x, int y)
{
    unsigned int h, i, j;

    for (h = 0; h < w; ++ h) {
		for (i = 0; i < x; ++ i) {
			for (j = 0; j < y; ++ j)
				free(mat4D_float[h][i][j]);
			free(mat4D_float[h][i]);
	    }
		free(mat4D_float[h]);
	}
	free(mat4D_float);
}

/*___________________________________________________________________________*/
/** 2D vector matrix */
/** allocate */
Vec **alloc_mat2D_vec(Vec **mat2D_vec, int x, int y)
{
    unsigned int i;

	mat2D_vec = (Vec **)safe_malloc(x * sizeof(Vec *));
	for (i = 0; i < x; ++ i)
		mat2D_vec[i] = (Vec *)safe_malloc(y * sizeof(Vec));

	return mat2D_vec;
}

/** initialise */
void init_mat2D_vec(Vec **mat2D_vec, int x, int y, Vec val)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			v_copy(&(mat2D_vec[i][j]), &(val));
}

/** free */
void free_mat2D_vec(Vec **mat2D_vec, int x)
{
    unsigned int i;

	for (i = 0; i < x; ++ i)
		free(mat2D_vec[i]);
    free(mat2D_vec);
}

/*___________________________________________________________________________*/
/** 3D vector matrix */
/** allocate */
Vec ***alloc_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_vec = (Vec ***)safe_malloc(x * sizeof(Vec **));
	for (i = 0; i < x; ++ i) {
		mat3D_vec[i] = (Vec **)safe_malloc(y * sizeof(Vec *));
		for (j = 0; j < y; ++ j)
			mat3D_vec[i][j] = (Vec *)safe_malloc(z * sizeof(Vec));
	}

	return mat3D_vec;
}

/** initialise */
void init_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z, Vec val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				v_copy(&(mat3D_vec[i][j][k]), &(val));
}

/** free */
void free_mat3D_vec(Vec ***mat3D_vec, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_vec[i][j]);
        free(mat3D_vec[i]);
    }
    free(mat3D_vec);
}

/*___________________________________________________________________________*/
/** 4D vector matrix */
/** allocate */
Vec ****alloc_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z)
{
	unsigned int h, i, j;

	mat4D_vec = (Vec ****)safe_malloc(w * sizeof(Vec ***));
	for (h = 0; h < w; ++ h) {
		mat4D_vec[h] = (Vec ***)safe_malloc(x * sizeof(Vec **));
		for (i = 0; i < x; ++ i) {
			mat4D_vec[h][i] = (Vec **)safe_malloc(y * sizeof(Vec *));
			for (j = 0; j < y; ++ j)
				mat4D_vec[h][i][j] = (Vec *)safe_malloc(z * sizeof(Vec));
		}
	}

	return mat4D_vec;
}

/** initialise */
void init_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z, Vec val)
{
    unsigned int h, i, j, k;

	for (h = 0; h < w; ++ h)
		for (i = 0; i < x; ++ i)
			for (j = 0; j < y; ++ j)
				for (k = 0; k < z; ++ k)
					v_copy(&(mat4D_vec[h][i][j][k]), &(val));
}

/** free */
void free_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y)
{
    unsigned int h, i, j;

    for (h = 0; h < w; ++ h) {
		for (i = 0; i < x; ++ i) {
			for (j = 0; j < y; ++ j)
				free(mat4D_vec[h][i][j]);
			free(mat4D_vec[h][i]);
	    }
		free(mat4D_vec[h]);
	}
	free(mat4D_vec);
}

