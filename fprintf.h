/*
 * fprintf.h
 *
 * Code generation for function 'fprintf'
 *
 */

#ifndef __FPRINTF_H__
#define __FPRINTF_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "test_vec_gen_types.h"

/* Function Declarations */
extern void b_fprintf(void);
extern void d_fprintf(void);
extern void f_fprintf(void);
extern void h_fprintf(double formatSpec);
extern void j_fprintf(double formatSpec);
extern void l_fprintf(int formatSpec, int varargin_1);
extern void n_fprintf(double fileID, const creal32_T varargin_1, const creal32_T
                      varargin_2, const creal32_T varargin_3, const creal32_T
                      varargin_4, const creal32_T varargin_5, const creal32_T
                      varargin_6, const creal32_T varargin_7, const creal32_T
                      varargin_8, const creal32_T varargin_9, const creal32_T
                      varargin_10, const creal32_T varargin_11, const creal32_T
                      varargin_12, const creal32_T varargin_13, const creal32_T
                      varargin_14, const creal32_T varargin_15, const creal32_T
                      varargin_16);

#endif

/* End of code generation (fprintf.h) */
