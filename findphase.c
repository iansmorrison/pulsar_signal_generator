/*
 * findphase.c
 *
 * Code generation for function 'findphase'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "test_vec_gen_emxutil.h"
#include "rdivide.h"
#include "test_vec_gen_rtwutil.h"
#include <stdio.h>

/* Function Definitions */
void findphase(const emxArray_real_T *t, int nbins, double period,
               emxArray_real_T *phase)
{
  unsigned int uv1[2];
  int k;
  emxArray_real_T *r;
  double b_r;
  emxArray_real_T *x;
  int b_x;

  /*  Determines the phase from the absolute time, 't', the number of bins */
  /*  per pulse period 'nbins', and the pulsar period. */
  /*   */
  /*  Inputs: */
  /*  -------- */
  /*    t       - vector of times (in seconds) */
  /*    nbins   - number of bins per pulsar period */
  /*    period  - pulsar period (s) */
  /*   */
  /*  Outputs: */
  /*  -------- */
  /*  */
  /*     dat -  pulsar "phase", given as integer between 1 and nbins */
  /*  */
  /*  Description: */
  /*  ------------ */
  /*  Calculates the pulsar phase corresponding to the requested times  */
  /*  by using a known pulsar calibration */
  /*   */
  /*  Changes: */
  /*  -------- */
  /*  */
  /*  Author           Date         Comments */
  /*  ---------------  -----------  ---------------------------------------- */
  /*  D. Hicks         21-Apr-2014  Original version */
  /*  ---------------------------------------------------------------------- */
  /* phase = 1 + transpose(round(mod(t, pcal.a)/pcal.a*(nbins-1))); */
  /*  fractional period */
  for (k = 0; k < 2; k++) {
    uv1[k] = (unsigned int)t->size[k];
  }

  emxInit_real_T(&r, 2);
  k = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = (int)uv1[1];
  emxEnsureCapacity((emxArray__common *)r, k, (int)sizeof(double));
  for (k = 0; k < (int)uv1[1]; k++) {
    if (period == 0.0) {
      b_r = t->data[k];
    } else if (period == floor(period)) {
      b_r = t->data[k] - floor(t->data[k] / period) * period;
    } else {
      b_r = t->data[k] / period;
      if (fabs(b_r - rt_roundd_snf(b_r)) <= 2.2204460492503131E-16 * fabs(b_r))
      {
        b_r = 0.0;
      } else {
        b_r = (b_r - floor(b_r)) * period;
      }
    }

    r->data[k] = b_r;
  }

  emxInit_real_T(&x, 2);

  /*  bin assignments are the center of each bin (Matlab indexing). */
  rdivide(r, period, x);
  k = x->size[0] * x->size[1];
  x->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)x, k, (int)sizeof(double));
  k = x->size[0];
  b_x = x->size[1];
  b_x *= k;
  for (k = 0; k < b_x; k++) {
    x->data[k] = x->data[k] * (double)nbins + 0.5;
  }

  k = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)r, k, (int)sizeof(double));
  b_x = x->size[0] * x->size[1];
  for (k = 0; k < b_x; k++) {
    r->data[k] = x->data[k];
  }

  for (k = 0; k < x->size[1]; k++) {
    r->data[k] = rt_roundd_snf(r->data[k]);
  }

  emxFree_real_T(&x);
  k = phase->size[0];
  phase->size[0] = r->size[1];
  emxEnsureCapacity((emxArray__common *)phase, k, (int)sizeof(double));
  b_x = r->size[1];
  for (k = 0; k < b_x; k++) {
    phase->data[k] = r->data[r->size[0] * k];
  }

  emxFree_real_T(&r);
}

/* End of code generation (findphase.c) */
