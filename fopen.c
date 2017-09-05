/*
 * fopen.c
 *
 * Code generation for function 'fopen'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fopen.h"
#include "fileManager.h"
#include <stdio.h>

/* Function Definitions */
double b_fopen(void)
{
  return b_fileManager();
}

double c_fopen(void)
{
  return f_fileManager();
}

/* End of code generation (fopen.c) */
