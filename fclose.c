/*
 * fclose.c
 *
 * Code generation for function 'fclose'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fclose.h"
#include "fileManager.h"
#include <stdio.h>

/* Function Definitions */
void b_fclose(double fileID)
{
  d_fileManager(fileID);
}

/* End of code generation (fclose.c) */
