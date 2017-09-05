/*
 * flipud.c
 *
 * Code generation for function 'flipud'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "flipud.h"
#include <stdio.h>

/* Function Definitions */
void flipud(emxArray_creal32_T *x)
{
  int m;
  int md2;
  int i;
  float xtmp_re;
  float xtmp_im;
  m = x->size[0];
  md2 = x->size[0] >> 1;
  for (i = 1; i <= md2; i++) {
    xtmp_re = x->data[i - 1].re;
    xtmp_im = x->data[i - 1].im;
    x->data[i - 1] = x->data[m - i];
    x->data[m - i].re = xtmp_re;
    x->data[m - i].im = xtmp_im;
  }
}

/* End of code generation (flipud.c) */
