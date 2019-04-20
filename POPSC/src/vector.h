/*==============================================================================
vector.h : Vector definitions
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>

/*___________________________________________________________________________*/
/* constants */
#define PI (4.0 * atan(1.0))
#define AN (180 / PI)

/*___________________________________________________________________________*/
/* vector is array of floats,
 * here including a set of transformed coordinates */
typedef struct
{
   float x, y, z;
   /*float tx, ty, tz;*/ /* disabled since 06.01.09 */
} Vec;

/*___________________________________________________________________________*/
/* prototypes */
float v_len(Vec *v);
float v_dot_pro(Vec *v1, Vec *v2);
void v_cro_pro(Vec *v3, Vec *v1, Vec *v2);
void vector_cro_pro(Vec *v1, Vec *v2, Vec *v3);
void v_sum(Vec *v3, Vec *v1, Vec *v2);
void vector_sum(Vec *v1, Vec *v2, Vec *v3);
void v_dif(Vec *v3, Vec *v1, Vec *v2);
void vector_dif(Vec *v1, Vec *v2, Vec *v3);
float v_ang(Vec *v1, Vec *v2);
void v_div_sca(Vec *v2, Vec *v1, float a);
void vector_div_sca(Vec *v1, float a, Vec *v2);
void v_mul_sca(Vec *v2, Vec *v1, float a);
void vector_mul_sca(Vec *v1, float a, Vec *v2);
void v_norm(Vec *v2, Vec *v1);
void vector_norm(Vec *v1, Vec *v2);
void v_zero(Vec *v);
void rotate_2D(float *dim1, float *dim2, float phi);
Vec v_shift_rotate_xy(Vec *v1, Vec *shift, float phi);
Vec v_shift_rotate_yz(Vec *v1, Vec *shift, float phi);
void v_copy(Vec *v2, Vec *v1);
void vector_copy(Vec *v1, Vec *v2);
float v_rmsd(Vec *v1, Vec *v2);
void v_put(Vec *v1);
void v_put_char(Vec *v1, char *str);

#endif
