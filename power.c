/*
 * power.c
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "power.h"
#include "test_vec_gen_emxutil.h"
#include <stdio.h>

/* Function Definitions */
void b_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  unsigned int unnamed_idx_0;
  int k;
  unnamed_idx_0 = (unsigned int)a->size[0];
  k = y->size[0];
  y->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k < (int)unnamed_idx_0; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  unsigned int uv0[2];
  int i2;
  int k;
  for (i2 = 0; i2 < 2; i2++) {
    uv0[i2] = (unsigned int)a->size[i2];
  }

  i2 = y->size[0] * y->size[1];
  y->size[0] = (int)uv0[0];
  y->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
  i2 = (int)uv0[0] * (int)uv0[1];
  for (k = 0; k < i2; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/* End of code generation (power.c) */
