/*
 * fread.c
 *
 * Code generation for function 'fread'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fread.h"
#include "fileManager.h"
#include <stdio.h>

/* Type Definitions */
#include <stddef.h>

/* Function Definitions */
void b_fread(double fileID, double A_data[], int A_size[1])
{
  size_t nBytes;
  FILE * filestar;
  boolean_T p;
  double b_A_data[4];
  int bytesOut;
  int numRead;
  size_t numReadSizeT;
  int i11;
  nBytes = sizeof(double);
  filestar = c_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  if (filestar == (FILE *)NULL) {
    A_size[0] = 0;
  } else {
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 4) && (numRead > 0)) {
      numReadSizeT = fread(&b_A_data[0], nBytes, 4 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 5; numRead++) {
      b_A_data[numRead] = 0.0;
    }

    A_size[0] = 4;
    for (i11 = 0; i11 < 4; i11++) {
      A_data[i11] = b_A_data[i11];
    }

    if (bytesOut < 4) {
      if (1 > bytesOut) {
        numRead = 0;
      } else {
        numRead = bytesOut;
      }

      A_size[0] = numRead;
      for (i11 = 0; i11 < numRead; i11++) {
        A_data[i11] = b_A_data[i11];
      }
    }
  }
}

/* End of code generation (fread.c) */
