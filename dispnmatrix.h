/*
 * dispnmatrix.h
 *
 * Code generation for function 'dispnmatrix'
 *
 */

#ifndef __DISPNMATRIX_H__
#define __DISPNMATRIX_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "test_vec_gen_types.h"

/* Function Declarations */
extern void b_dispnmatrix(const double frange[2], double usable_BW, int Nout,
  double DD, double Tout, emxArray_creal_T *H, double *f, double *n_hi, double
  *n_lo);
extern void dispnmatrix(const double frange[2], double usable_BW, int Nout, int
  nfreq, double DD, double Tout, double direc, emxArray_creal_T *H,
  emxArray_real_T *f, double *n_hi, double *n_lo);

#endif

/* End of code generation (dispnmatrix.h) */
