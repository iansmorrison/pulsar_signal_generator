/*
 * fileManager.c
 *
 * Code generation for function 'fileManager'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fileManager.h"
#include "test_vec_gen_rtwutil.h"
#include <stdio.h>

/* Variable Definitions */
static FILE * eml_openfiles[20];
static boolean_T eml_autoflush[20];

/* Function Declarations */
static FILE * e_fileManager(signed char varargin_1);
static signed char filedata(void);
static void getfilestar(double fid, FILE * *filestar, boolean_T *autoflush);

/* Function Definitions */
static FILE * e_fileManager(signed char varargin_1)
{
  FILE * f;
  signed char fileid;
  fileid = varargin_1;
  if (varargin_1 < 0) {
    fileid = -1;
  }

  if (fileid >= 3) {
    f = eml_openfiles[fileid - 3];
  } else if (fileid == 0) {
    f = stdin;
  } else if (fileid == 1) {
    f = stdout;
  } else if (fileid == 2) {
    f = stderr;
  } else {
    f = NULL;
  }

  return f;
}

static signed char filedata(void)
{
  signed char f;
  signed char k;
  boolean_T exitg1;
  f = 0;
  k = 1;
  exitg1 = false;
  while ((!exitg1) && (k < 21)) {
    if (eml_openfiles[k - 1] == (FILE *)NULL) {
      f = k;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return f;
}

static void getfilestar(double fid, FILE * *filestar, boolean_T *autoflush)
{
  signed char fileid;
  fileid = (signed char)rt_roundd_snf(fid);
  if ((fileid < 0) || (fid != fileid)) {
    fileid = -1;
  }

  if (fileid >= 3) {
    fileid = (signed char)(fileid - 2);
    *filestar = eml_openfiles[fileid - 1];
    *autoflush = eml_autoflush[fileid - 1];
  } else if (fileid == 0) {
    *filestar = stdin;
    *autoflush = true;
  } else if (fileid == 1) {
    *filestar = stdout;
    *autoflush = true;
  } else if (fileid == 2) {
    *filestar = stderr;
    *autoflush = true;
  } else {
    *filestar = NULL;
    *autoflush = true;
  }
}

signed char b_fileManager(void)
{
  signed char f;
  signed char j;
  char cv2[10];
  int i10;
  static const char cv3[10] = { 'p', 'u', 'l', 's', 'e', '.', 'b', 'i', 'n',
    '\x00' };

  char cv4[3];
  static const char cv5[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i10 = 0; i10 < 10; i10++) {
      cv2[i10] = cv3[i10];
    }

    for (i10 = 0; i10 < 3; i10++) {
      cv4[i10] = cv5[i10];
    }

    filestar = fopen(cv2, cv4);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i10 = j + 2;
      if (i10 > 127) {
        i10 = 127;
      }

      f = (signed char)i10;
    }
  }

  return f;
}

FILE * c_fileManager(double varargin_1)
{
  FILE * f;
  boolean_T a;
  getfilestar(varargin_1, &f, &a);
  return f;
}

int d_fileManager(double varargin_1)
{
  int f;
  signed char fileid;
  FILE * filestar;
  int cst;
  f = -1;
  fileid = (signed char)rt_roundd_snf(varargin_1);
  if ((fileid < 0) || (varargin_1 != fileid)) {
    fileid = -1;
  }

  filestar = e_fileManager(fileid);
  if ((filestar == (FILE *)NULL) || (fileid < 3)) {
  } else {
    cst = fclose(filestar);
    if (cst == 0) {
      f = 0;
      fileid = (signed char)(fileid - 2);
      eml_openfiles[fileid - 1] = NULL;
      eml_autoflush[fileid - 1] = true;
    }
  }

  return f;
}

signed char f_fileManager(void)
{
  signed char f;
  signed char j;
  char cv6[16];
  int i12;
  static const char cv7[16] = { 't', 'e', 's', 't', '_', 'v', 'e', 'c', 't', 'o',
    'r', '.', 't', 'x', 't', '\x00' };

  char cv8[3];
  static const char cv9[3] = { 'w', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i12 = 0; i12 < 16; i12++) {
      cv6[i12] = cv7[i12];
    }

    for (i12 = 0; i12 < 3; i12++) {
      cv8[i12] = cv9[i12];
    }

    filestar = fopen(cv6, cv8);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i12 = j + 2;
      if (i12 > 127) {
        i12 = 127;
      }

      f = (signed char)i12;
    }
  }

  return f;
}

void fileManager(FILE * *f, boolean_T *a)
{
  *f = stdout;
  *a = true;
}

void filedata_init(void)
{
  FILE * a;
  int i;
  a = NULL;
  for (i = 0; i < 20; i++) {
    eml_autoflush[i] = false;
    eml_openfiles[i] = a;
  }
}

void g_fileManager(double varargin_1, FILE * *f, boolean_T *a)
{
  getfilestar(varargin_1, f, a);
}

/* End of code generation (fileManager.c) */
