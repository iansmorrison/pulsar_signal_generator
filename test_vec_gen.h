/*
 * test_vec_gen.h
 *
 * Code generation for function 'test_vec_gen'
 *
 */

#ifndef __TEST_VEC_GEN_H__
#define __TEST_VEC_GEN_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "test_vec_gen_types.h"

/* Function Declarations */
extern void test_vec_gen(int nbins, double period, double t0, int Nout, int
  nseries, int format, int shift, double f0, double f_sample_out, double DM, int
  out_type);

#endif

/* End of code generation (test_vec_gen.h) */
