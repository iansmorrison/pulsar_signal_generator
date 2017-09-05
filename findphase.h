/*
 * findphase.h
 *
 * Code generation for function 'findphase'
 *
 */

#ifndef __FINDPHASE_H__
#define __FINDPHASE_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "test_vec_gen_types.h"

/* Function Declarations */
extern void findphase(const emxArray_real_T *t, int nbins, double period,
                      emxArray_real_T *phase);

#endif

/* End of code generation (findphase.h) */
