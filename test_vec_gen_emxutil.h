/*
 * test_vec_gen_emxutil.h
 *
 * Code generation for function 'test_vec_gen_emxutil'
 *
 */

#ifndef __TEST_VEC_GEN_EMXUTIL_H__
#define __TEST_VEC_GEN_EMXUTIL_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "test_vec_gen_types.h"

/* Function Declarations */
extern void b_emxInit_creal32_T(emxArray_creal32_T **pEmxArray, int
  numDimensions);
extern void b_emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);
extern void b_emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
extern void b_emxInit_real32_T(emxArray_real32_T **pEmxArray, int numDimensions);
extern void b_emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
extern void emxFree_creal32_T(emxArray_creal32_T **pEmxArray);
extern void emxFree_creal_T(emxArray_creal_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real32_T(emxArray_real32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
extern void emxInit_creal32_T(emxArray_creal32_T **pEmxArray, int numDimensions);
extern void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
extern void emxInit_real32_T(emxArray_real32_T **pEmxArray, int numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/* End of code generation (test_vec_gen_emxutil.h) */
