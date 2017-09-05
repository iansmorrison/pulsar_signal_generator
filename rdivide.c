/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "rdivide.h"
#include "test_vec_gen_emxutil.h"
#include <stdio.h>

/* Function Definitions */
void b_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
               emxArray_real_T *z)
{
  int i3;
  int loop_ub;
  i3 = z->size[0] * z->size[1];
  z->size[0] = x->size[0];
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i3, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    z->data[i3] = x->data[i3] / y->data[i3];
  }
}

void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
               emxArray_real_T *z)
{
  int i9;
  int loop_ub;
  i9 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i9, (int)sizeof(double));
  loop_ub = x->size[0];
  for (i9 = 0; i9 < loop_ub; i9++) {
    z->data[i9] = x->data[i9] / y->data[i9];
  }
}

void rdivide(const emxArray_real_T *x, double y, emxArray_real_T *z)
{
  int i1;
  int loop_ub;
  i1 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i1, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    z->data[i1] = x->data[i1] / y;
  }
}

/* End of code generation (rdivide.c) */
