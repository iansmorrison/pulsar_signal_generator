/*
 * rng.c
 *
 * Code generation for function 'rng'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "rng.h"
#include "test_vec_gen_data.h"
#include <stdio.h>

/* Function Definitions */
void rng(void)
{
  unsigned int r;
  int mti;
  r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = (r ^ r >> 30U) * 1812433253U + (1 + mti);
    state[mti + 1] = r;
  }

  state[624] = 624U;
}

/* End of code generation (rng.c) */
