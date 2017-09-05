/*
 * test_vec_gen_initialize.c
 *
 * Code generation for function 'test_vec_gen_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "test_vec_gen_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "fileManager.h"
#include <stdio.h>

/* Function Definitions */
void test_vec_gen_initialize(void)
{
  rt_InitInfAndNaN(8U);
  filedata_init();
  c_eml_rand_mt19937ar_stateful_i();
}

/* End of code generation (test_vec_gen_initialize.c) */
