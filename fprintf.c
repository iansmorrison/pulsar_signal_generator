/*
 * fprintf.c
 *
 * Code generation for function 'fprintf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fprintf.h"
#include "fileManager.h"
#include <stdio.h>

/* Function Declarations */
static double c_fprintf(void);
static double e_fprintf(void);
static double g_fprintf(void);
static double i_fprintf(double varargin_1);
static double k_fprintf(double varargin_1);
static double m_fprintf(int varargin_1, int varargin_2);

/* Function Definitions */
static double c_fprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[35] = { 'W', 'a', 'r', 'n', 'i', 'n', 'g', ':', ' ',
    'd', 'i', 'r', 'e', 'c', ' ', 'c', 'a', 'n', ' ', 'o', 'n', 'l', 'y', ' ',
    'b', 'e', ' ', '1', ' ', 'o', 'r', ' ', '-', '1', '\x00' };

  nbytesint = 0;
  b_NULL = NULL;
  fileManager(&filestar, &autoflush);
  if (filestar == b_NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

static double e_fprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[48] = { 'D', 'i', 's', 'c', 'a', 'r', 'd', 'i', 'n',
    'g', ' ', 'r', 'e', 'g', 'i', 'o', 'n', 's', ' ', 'o', 'u', 't', 's', 'i',
    'd', 'e', ' ', 't', 'h', 'e', ' ', 'u', 's', 'a', 'b', 'l', 'e', ' ', 'b',
    'a', 'n', 'd', 'w', 'i', 'd', 't', 'h', '\x00' };

  nbytesint = 0;
  b_NULL = NULL;
  fileManager(&filestar, &autoflush);
  if (filestar == b_NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

static double g_fprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[44] = { 'W', 'a', 'r', 'n', 'i', 'n', 'g', ':', ' ',
    't', 'i', 'm', 'e', ' ', 's', 'e', 'r', 'i', 'e', 's', ' ', 't', 'o', 'o',
    ' ', 's', 'h', 'o', 'r', 't', ' ', 'f', 'o', 'r', ' ', 'g', 'i', 'v', 'e',
    'n', ' ', 'D', 'M', '\x00' };

  nbytesint = 0;
  b_NULL = NULL;
  fileManager(&filestar, &autoflush);
  if (filestar == b_NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

static double i_fprintf(double varargin_1)
{
  int nbytesint;
  FILE * b_NULL;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[35] = { 'L', 'o', 's', 't', ' ', 'f', 'r', 'a', 'c',
    't', 'i', 'o', 'n', ' ', 'o', 'f', ' ', 't', 'i', 'm', 'e', ' ', 's', 'e',
    'r', 'i', 'e', 's', ' ', '=', ' ', '%', 'f', '\x0a', '\x00' };

  nbytesint = 0;
  b_NULL = NULL;
  fileManager(&filestar, &autoflush);
  if (filestar == b_NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt, varargin_1);
    fflush(filestar);
  }

  return nbytesint;
}

static double k_fprintf(double varargin_1)
{
  int nbytesint;
  FILE * b_NULL;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[27] = { 'T', 'i', 'm', 'e', ' ', 's', 'e', 'r', 'i',
    'e', 's', ' ', 'l', 'e', 'n', 'g', 't', 'h', ' ', '=', ' ', '%', 'f', ' ',
    's', '\x0a', '\x00' };

  nbytesint = 0;
  b_NULL = NULL;
  fileManager(&filestar, &autoflush);
  if (filestar == b_NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt, varargin_1);
    fflush(filestar);
  }

  return nbytesint;
}

static double m_fprintf(int varargin_1, int varargin_2)
{
  int nbytesint;
  FILE * b_NULL;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[17] = { 'L', 'o', 'o', 'p', ' ', '#', ' ', '%', 'i',
    ' ', 'o', 'f', ' ', '%', 'i', '\x0a', '\x00' };

  nbytesint = 0;
  b_NULL = NULL;
  fileManager(&filestar, &autoflush);
  if (filestar == b_NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2);
    fflush(filestar);
  }

  return nbytesint;
}

void b_fprintf(void)
{
  c_fprintf();
}

void d_fprintf(void)
{
  e_fprintf();
}

void f_fprintf(void)
{
  g_fprintf();
}

void h_fprintf(double formatSpec)
{
  i_fprintf(formatSpec);
}

void j_fprintf(double formatSpec)
{
  k_fprintf(formatSpec);
}

void l_fprintf(int formatSpec, int varargin_1)
{
  m_fprintf(formatSpec, varargin_1);
}

void n_fprintf(double fileID, const creal32_T varargin_1, const creal32_T
               varargin_2, const creal32_T varargin_3, const creal32_T
               varargin_4, const creal32_T varargin_5, const creal32_T
               varargin_6, const creal32_T varargin_7, const creal32_T
               varargin_8, const creal32_T varargin_9, const creal32_T
               varargin_10, const creal32_T varargin_11, const creal32_T
               varargin_12, const creal32_T varargin_13, const creal32_T
               varargin_14, const creal32_T varargin_15, const creal32_T
               varargin_16)
{
  FILE * b_NULL;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[97] = { '%', '.', '1', '5', 'f', ' ', '%', '.', '1',
    '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ',
    '%', '.', '1', '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ', '%', '.', '1',
    '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ',
    '%', '.', '1', '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ', '%', '.', '1',
    '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ',
    '%', '.', '1', '5', 'f', ' ', '%', '.', '1', '5', 'f', ' ', '\x00' };

  b_NULL = NULL;
  g_fileManager(fileID, &filestar, &autoflush);
  if (filestar == b_NULL) {
  } else {
    fprintf(filestar, cfmt, varargin_1.re, varargin_2.re, varargin_3.re,
            varargin_4.re, varargin_5.re, varargin_6.re, varargin_7.re,
            varargin_8.re, varargin_9.re, varargin_10.re, varargin_11.re,
            varargin_12.re, varargin_13.re, varargin_14.re, varargin_15.re,
            varargin_16.re);
    if (autoflush) {
      fflush(filestar);
    }
  }
}

/* End of code generation (fprintf.c) */
