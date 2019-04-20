/*==============================================================================
vector.c : Vector definitions
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "stdio.h"
#include "vector.h"

/* vector length |v| */
float v_len(Vec *v)
{
    return sqrt(pow(v->x, 2) + pow(v->y, 2) + pow(v->z, 2));
}

/*____________________________________________________________________________*/
/* vector dot product v1 * v2 */
float v_dot_pro(Vec *v1, Vec *v2)
{
    return (v1->x * v2->x + v1->y * v2->y + v1->z * v2->z);
}

/*____________________________________________________________________________*/
/* vector cross product v1 x v2 */
void v_cro_pro(Vec *v3, Vec *v1, Vec *v2)
{
    v3->x = v1->y * v2->z - v1->z * v2->y;
    v3->y = v1->z * v2->x - v1->x * v2->z;
    v3->z = v1->x * v2->y - v1->y * v2->x;
}

/*____________________________________________________________________________*/
/* vector cross product v1 x v2 */
void vector_cro_pro(Vec *v1, Vec *v2, Vec *v3)
{
    v3->x = v1->y * v2->z - v1->z * v2->y;
    v3->y = v1->z * v2->x - v1->x * v2->z;
    v3->z = v1->x * v2->y - v1->y * v2->x;
}

/*____________________________________________________________________________*/
/* vector sum */
void v_sum(Vec *v3, Vec *v1, Vec *v2)
{
    v3->x = v1->x + v2->x;
    v3->y = v1->y + v2->y;
    v3->z = v1->z + v2->z;
}

/*____________________________________________________________________________*/
/* vector sum */
void vector_sum(Vec *v1, Vec *v2, Vec *v3)
{
    v3->x = v1->x + v2->x;
    v3->y = v1->y + v2->y;
    v3->z = v1->z + v2->z;
}

/*____________________________________________________________________________*/
/* difference vector */
void v_dif(Vec *v3, Vec *v1, Vec *v2)
{
    v3->x = v1->x - v2->x;
    v3->y = v1->y - v2->y;
    v3->z = v1->z - v2->z;
}

/*____________________________________________________________________________*/
/* difference vector */
void vector_dif(Vec *v1, Vec *v2, Vec *v3)
{
    v3->x = v1->x - v2->x;
    v3->y = v1->y - v2->y;
    v3->z = v1->z - v2->z;
}

/*____________________________________________________________________________*/
/* vector angle */
float v_ang(Vec *v1, Vec *v2)
{
	float cos_angle;
	/*float c_angle, s_angle, sin_angle;*/
    /*Vec v3;*/
    /* cosinus of vector angle from dot product */
    /* cos angle = v1 * v2 / ||v1|| * ||v2|| */
    /* angle = acos (cos angle) */
    cos_angle = v_dot_pro(v1, v2) / (v_len(v1) * v_len(v2));
    /*c_angle = acos(cos_angle);*/
    /* sinus of vector angle from cross product */
    /* sin angle = length (v1 x v2) / (length (v1) * length (v2)) */
    /*v_cro_pro(&v3, v1, v2);*/
    /*sin_angle = v_len(&v3) / (v_len(v1) * v_len(v2));*/
    /*s_angle = asin(sin_angle);*/

    return (AN * acos(cos_angle));
}

/*____________________________________________________________________________*/
/* divide vector by scalar */
void v_div_sca(Vec *v2, Vec *v1, float a)
{
	v2->x = v1->x / a;
	v2->y = v1->y / a;
	v2->z = v1->z / a;
}

/*____________________________________________________________________________*/
/* divide vector by scalar */
void vector_div_sca(Vec *v1, float a, Vec *v2)
{
	v2->x = v1->x / a;
	v2->y = v1->y / a;
	v2->z = v1->z / a;
}

/*____________________________________________________________________________*/
/* multiply vector by scalar */
void v_mul_sca(Vec *v2, Vec *v1, float a)
{
	v2->x = v1->x * a;
	v2->y = v1->y * a;
	v2->z = v1->z * a;
}

/*____________________________________________________________________________*/
/* divide vector by scalar */
void vector_mul_sca(Vec *v1, float a, Vec *v2)
{
	v2->x = v1->x * a;
	v2->y = v1->y * a;
	v2->z = v1->z * a;
}

/*____________________________________________________________________________*/
/* normalise vector */
void v_norm(Vec *v2, Vec *v1)
{
	if (v_len(v1) > 0)
		v_div_sca(v2, v1, v_len(v1));
	else
		v_zero(v2);
}

/*____________________________________________________________________________*/
/* normalise vector */
void vector_norm(Vec *v1, Vec *v2)
{
	if (v_len(v1) > 0)
		vector_div_sca(v1, v_len(v1), v2);
	else
		v_zero(v2);
}

/*____________________________________________________________________________*/
/* zero vector */
void v_zero(Vec *v)
{
	v->x = 0;
	v->y = 0;
	v->z = 0;
}

/*____________________________________________________________________________*/
/* 2D rotation of vector in origin:
	x' = x cos(phi) + y sin(phi)
	y' = -x sin(phi) + y cos(phi) */
void rotate_2D(float *dim1, float *dim2, float phi)
{
	float tdim1 = 0;
	float tdim2 = 0;

	/* rotate */
	tdim1 =  *dim1 * cos(phi) + *dim2 * sin(phi);
	tdim2 = -*dim1 * sin(phi) + *dim2 * cos(phi); 

	/* assign */
	*dim1 = tdim1;
	*dim2 = tdim2;
}

/*____________________________________________________________________________*/
/* shift vector to origin, rotate in xy plane and shift back */
Vec v_shift_rotate_xy(Vec *v1, Vec *shift, float phi)
{
	Vec v2;

	v_dif(&v2, v1, shift);
	rotate_2D(&(v2.x), &(v2.y), phi);
	v_sum(&v2, shift, &v2);

	return v2;
}

/*____________________________________________________________________________*/
/* shift vector to origin, rotate in yz plane and shift back */
Vec v_shift_rotate_yz(Vec *v1, Vec *shift, float phi)
{
	Vec v2;

	v_dif(&v2, v1, shift);
	rotate_2D(&(v2.y), &(v2.z), phi);
	v_sum(&v2, shift, &v2);

	return v2;
}

/*____________________________________________________________________________*/
/* copy vector */
void v_copy(Vec *v2, Vec *v1)
{
	v2->x = v1->x;
	v2->y = v1->y;
	v2->z = v1->z;
}

/*____________________________________________________________________________*/
/* copy vector */
void vector_copy(Vec *v1, Vec *v2)
{
	v2->x = v1->x;
	v2->y = v1->y;
	v2->z = v1->z;
}

/*____________________________________________________________________________*/
/* RMSD of two vectors, the second in the transformed form! */
float v_rmsd(Vec *v1, Vec *v2)
{
	return sqrt(pow(v1->x - v2->x, 2)
              + pow(v1->y - v2->y, 2)
			  + pow(v1->z - v2->z, 2));
}

/*____________________________________________________________________________*/
/** print vector coordinates to stderr */
void v_put(Vec *v1)
{
	fprintf(stderr, "%f %f %f", v1->x, v1->y, v1->z);
}

/*____________________________________________________________________________*/
/** print vector coordinates to string */
void v_put_char(Vec *v1, char *str)
{
	sprintf(str, "%4e %4e %4e", v1->x, v1->y, v1->z);
}

