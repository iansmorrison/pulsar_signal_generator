/*
 * fftshift.c
 *
 * Code generation for function 'fftshift'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fftshift.h"
#include <stdio.h>

/* Function Definitions */
void fftshift(emxArray_creal32_T *x)
{
  int dim;
  int i1;
  int vlend2;
  int vstride;
  int k;
  int midoffset;
  int j;
  int ia;
  int ib;
  float xtmp_re;
  float xtmp_im;
  int ic;
  for (dim = 0; dim < 2; dim++) {
    if (dim + 1 <= 1) {
      i1 = x->size[0];
    } else {
      i1 = 1;
    }

    if (i1 <= 1) {
    } else {
      vlend2 = i1 / 2;
      vstride = 1;
      k = 1;
      while (k <= dim) {
        vstride *= x->size[0];
        k = 2;
      }

      midoffset = vlend2 * vstride;
      if (vlend2 << 1 == i1) {
        i1 = -1;
        for (j = 1; j <= vstride; j++) {
          i1++;
          ia = i1;
          ib = i1 + midoffset;
          for (k = 1; k <= vlend2; k++) {
            xtmp_re = x->data[ia].re;
            xtmp_im = x->data[ia].im;
            x->data[ia] = x->data[ib];
            x->data[ib].re = xtmp_re;
            x->data[ib].im = xtmp_im;
            ia += vstride;
            ib += vstride;
          }
        }
      } else {
        i1 = -1;
        for (j = 1; j <= vstride; j++) {
          i1++;
          ia = i1;
          ib = (i1 + midoffset) + 1;
          xtmp_re = x->data[ib - 1].re;
          xtmp_im = x->data[ib - 1].im;
          for (k = 1; k <= vlend2; k++) {
            ic = ib + vstride;
            x->data[ib - 1] = x->data[ia];
            x->data[ia] = x->data[ic - 1];
            ia += vstride;
            ib = ic;
          }

          x->data[ib - 1].re = xtmp_re;
          x->data[ib - 1].im = xtmp_im;
        }
      }
    }
  }
}

/* End of code generation (fftshift.c) */
