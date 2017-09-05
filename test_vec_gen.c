/*
 * test_vec_gen.c
 *
 * Code generation for function 'test_vec_gen'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fread.h"
#include "test_vec_gen_emxutil.h"
#include "randn.h"
#include "fprintf.h"
#include "ifft.h"
#include "fftshift.h"
#include "flipud.h"
#include "fft.h"
#include "fclose.h"
#include "fopen.h"
#include "rng.h"
#include "test_vec_gen_rtwutil.h"
#include <stdio.h>

/* Function Declarations */
static boolean_T eml_strcmp(const char a_data[], const int a_size[2]);
static double eml_switch_helper(const char expr_data[], const int expr_size[2]);

/* Function Definitions */
static boolean_T eml_strcmp(const char a_data[], const int a_size[2])
{
  boolean_T b_bool;
  int k;
  int32_T exitg1;
  static const char cv11[13] = { 'c', 'o', 'm', 'p', 'l', 'e', 'x', 't', 'o',
    'r', 'e', 'a', 'l' };

  b_bool = false;
  if (a_size[1] != 13) {
  } else {
    k = 0;
    do {
      exitg1 = 0;
      if (k <= 12) {
        if (a_data[k] != cv11[k]) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return b_bool;
}

static double eml_switch_helper(const char expr_data[], const int expr_size[2])
{
  double b_index;
  boolean_T b_bool;
  int k;
  int32_T exitg1;
  static const char cv10[13] = { 'c', 'o', 'm', 'p', 'l', 'e', 'x', 't', 'o',
    'r', 'e', 'a', 'l' };

  b_bool = false;
  if (expr_size[1] != 13) {
  } else {
    k = 0;
    do {
      exitg1 = 0;
      if (k <= 12) {
        if (expr_data[k] != cv10[k]) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    b_index = 0.0;
  } else {
    b_index = -1.0;
  }

  return b_index;
}

void test_vec_gen(int nbins, double period, double t0, int Nout, int nseries,
                  int format, int shift, double f0, double f_sample_out, double
                  DM, int out_type)
{
  int dformat_size[2];
  int i4;
  static const char cv0[13] = { 'c', 'o', 'm', 'p', 'l', 'e', 'x', 't', 'o', 'r',
    'e', 'a', 'l' };

  char dformat_data[16];
  int Nmul;
  static const char cv1[16] = { 'c', 'o', 'm', 'p', 'l', 'e', 'x', 't', 'o', 'c',
    'o', 'm', 'p', 'l', 'e', 'x' };

  double df;
  double Tin;
  int Nin;
  double Pmul;
  double b_df[2];
  double c_df[2];
  emxArray_real_T *S;
  emxArray_creal_T *H;
  double n_lo;
  double n_hi;
  double x;
  int nclip_out;
  int qEnd;
  int ii;
  int tmp_size[1];
  double tmp_data[4];
  emxArray_real_T *S0;
  emxArray_real_T *S1;
  emxArray_real_T *S2;
  emxArray_real_T *S3;
  emxArray_real_T *Jxx;
  emxArray_creal_T *Jxy;
  emxArray_creal_T *Jyx;
  emxArray_creal_T *b_Jxx;
  emxArray_creal_T *b_S0;
  int pEnd;
  int q;
  int p;
  emxArray_creal_T *J;
  int yk;
  long i5;
  int n;
  emxArray_int32_T *y;
  int k;
  emxArray_real_T *trel;
  double out_fid;
  emxArray_real_T *ttclip;
  emxArray_real_T *tt;
  emxArray_creal32_T *z;
  emxArray_creal_T *z0;
  emxArray_creal_T *zjj;
  emxArray_creal32_T *f1a;
  emxArray_creal32_T *f2a;
  emxArray_creal32_T *f1;
  emxArray_creal32_T *z1;
  emxArray_creal32_T *z2;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_real32_T *r3;
  emxArray_real32_T *r4;
  emxArray_real32_T *r5;
  emxArray_real32_T *b_y;
  emxArray_real32_T *b_x;
  emxArray_creal32_T *c_x;
  emxArray_boolean_T *d_x;
  emxArray_int32_T *idx;
  emxArray_int32_T *idx0;
  emxArray_creal32_T *b_z;
  emxArray_creal32_T *c_z;
  emxArray_creal_T *b_S3;
  emxArray_creal_T *r6;
  emxArray_creal_T *c_S3;
  int na;
  unsigned int unnamed_idx_0;
  boolean_T b_p;
  int i;
  int j;
  int kEnd;
  int32_T exitg4;
  int exponent;
  int i6;
  int jj;
  creal_T Jcoh[4];
  boolean_T exitg3;
  boolean_T guard1 = false;
  int32_T exitg1;
  boolean_T exitg2;
  double c_re;
  double c_im;
  float f1a_re;
  float f1a_im;
  float H_re;
  float H_im;
  creal32_T b_f1a;
  creal32_T c_f1a;
  creal32_T d_f1a;
  creal32_T e_f1a;
  creal32_T f_f1a;
  creal32_T g_f1a;
  creal32_T h_f1a;
  creal32_T i_f1a;
  creal32_T j_f1a;
  creal32_T k_f1a;
  creal32_T l_f1a;
  creal32_T m_f1a;
  creal32_T n_f1a;
  creal32_T o_f1a;
  creal32_T p_f1a;
  creal32_T q_f1a;

  /*  Generates a file containing dual noise vectors with phase-dependent */
  /*  partial polarization. Takes its pulse profile from the ASCII file */
  /*  "pulse.txt" which contains nbins rows of full Stokes parameters. */
  /*  Output file is 32-bit floating point with polarizations interleaved at */
  /*  each time step.  */
  /*  */
  /*  Derived from the signalgen() function developed by Damien Hicks. */
  /*  DATA SETTINGS */
  /*  */
  /*  fname     - Ouput filename */
  /*  hdrsize   - Header size */
  /*  hdrtype   - Data type for header ('uint8' = byte) */
  /*  ntype     - Data type for each element in a pair ('single' = float) */
  /*  Nout      - Length of each output vector */
  /*  nbins     - Number of bins within a pulse period */
  /*  npol      - Number of polarizations (should always be 2 when calc Stokes) */
  /*  nseries   - Number of forward FFT's to perform */
  /*  format    - Specifies conversion TO real or complex data */
  /*  shift     - Select whether an fftshift is used before the inverse FFT */
  /*              (don't shift if PFB is in the signal chain) */
  /*  out_type  - Type of output file: 0 for ASCII text, 1 for binary */
  /*  */
  /*  INSTRUMENT SETTINGS */
  /*  f0        - Centre frequency (MHz)  */
  /*  f_sample_out - Sampling frequency of output data (MHz) */
  /*   */
  /*  PULSAR SETTINGS */
  /*  Dconst    - Dispersion constant, s.MHz^2/(pc/cm^3) */
  /*  DM        - Dispersion measure, pc/cm^3 */
  /*  period    - Pulsar period (s) */
  /*  t0        - Time offset (s) for first sample */
  /*  */
  /*  OUTPUTS: */
  /*  -------- */
  /*  */
  /*     fname -  file containing two interleaved floating point test vectors */
  /*  */
  /*  Description: */
  /*  ------------ */
  /*  Generates a file containing dual noise vectors with phase-dependent */
  /*  partial polarization. File is 32-bit floating point with polarizations */
  /*  interleaved at each time step.  */
  /*   */
  /*  Changes: */
  /*  -------- */
  /*  */
  /*  Author           Date         Comments */
  /*  ---------------  -----------  ---------------------------------------- */
  /*  I. Morrison      04-Nov-2015  Changed to read pulse profile from file. */
  /*  */
  /*  ---------------------------------------------------------------------- */
  /*  define fwrite() as extrinsic for stand-alone code generation (not supported) */
  /*  start random number generation from a known seed. */
  /* rng('default'); */
  rng();

  /*  Header size */
  /*  Data type for header ('uint8' = byte) */
  /*  Data type for each element in a pair ('single' = float) */
  /* Nout = 2^20; %Length of each output vector */
  /* nbins = 2^10; % Number of bins within a pulse period */
  /*  Number of polarizations (should always be 2 when calc Stokes) */
  /* nseries = int32(30); % Number of FFT's to perform */
  /* noise = 0.0;  % 0.0 for no noise, 1.0 for noise (max(S/N)=1) */
  /* dformat = 'complextoreal'; %specifies conversion TO real or complex data */
  /* dformat = 'complextocomplex'; %specifies conversion TO real or complex data */
  /* shift = 1; % performs an fftshift before the inverse FFT */
  /*  Instrument settings */
  /* f0 = 1405.; % Centre frequency (MHz) */
  /*  Set bandwidth to 8 x 10 MHz, for testing with 8-channel channelizer */
  /* f_sample_out = 10.; % Sampling frequency of output (MHz) */
  /*  Pulsar settings */
  /*  s.MHz^2/(pc/cm^3) */
  /* =============== */
  /* Multiplying factor going from input to output type */
  if (format == 0) {
    dformat_size[0] = 1;
    dformat_size[1] = 13;
    for (i4 = 0; i4 < 13; i4++) {
      dformat_data[i4] = cv0[i4];
    }

    Nmul = 2;
  } else {
    dformat_size[0] = 1;
    dformat_size[1] = 16;
    for (i4 = 0; i4 < 16; i4++) {
      dformat_data[i4] = cv1[i4];
    }

    Nmul = 1;
  }

  /*  Sample spacing of output (seconds) */
  df = f_sample_out / (double)Nmul;

  /*  Bandwidth/Nyquist frequency (MHz) */
  Tin = 1.0 / fabs(f_sample_out) * 1.0E-6 * (double)Nmul;

  /*  Time spacing between input data elements */
  Nin = (int)rt_roundd_snf((double)Nout / (double)Nmul);

  /*  Number of data elements in input time series */
  Pmul = 1.0 / (double)Nmul;

  /*  Power multiplication factor for all but the DC channel */
  /* =============== */
  /*  Create the dispersion kernel and determine the number of elements to be */
  /*  clipped off the beginning and end. */
  /*  get sign of df (sign() function doesn't work stand-alone) */
  if (df <= 0.0) {
  } else {
    /*  Get matrix to perform dispersion on complex array */
    b_df[0] = -df / 2.0;
    b_df[1] = df / 2.0;
    for (i4 = 0; i4 < 2; i4++) {
      c_df[i4] = b_df[i4] + f0;
    }

    emxInit_real_T(&S, 2);
    emxInit_creal_T(&H, 1);
    b_dispnmatrix(c_df, df, Nin, 4148.804 * DM, Tin, H, &x, &n_hi, &n_lo);

    /*  Calculate the number of elements in the clipped input array */
    /*  Calculate number of elements in the clipped output array */
    df = rt_roundd_snf((double)Nout - n_lo * (double)Nmul);
    if (df < 2.147483648E+9) {
      if (df >= -2.147483648E+9) {
        i4 = (int)df;
      } else {
        i4 = MIN_int32_T;
      }
    } else if (df >= 2.147483648E+9) {
      i4 = MAX_int32_T;
    } else {
      i4 = 0;
    }

    df = rt_roundd_snf((double)i4 - n_hi * (double)Nmul);
    if (df < 2.147483648E+9) {
      if (df >= -2.147483648E+9) {
        nclip_out = (int)df;
      } else {
        nclip_out = MIN_int32_T;
      }
    } else if (df >= 2.147483648E+9) {
      nclip_out = MAX_int32_T;
    } else {
      nclip_out = 0;
    }

    /*  fraction of array that's lost */
    h_fprintf((n_lo + n_hi) / (double)Nin);
    df = rt_roundd_snf((double)Nin - n_lo);
    if (df < 2.147483648E+9) {
      if (df >= -2.147483648E+9) {
        i4 = (int)df;
      } else {
        i4 = MIN_int32_T;
      }
    } else if (df >= 2.147483648E+9) {
      i4 = MAX_int32_T;
    } else {
      i4 = 0;
    }

    df = rt_roundd_snf((double)i4 - n_hi);
    if (df < 2.147483648E+9) {
      if (df >= -2.147483648E+9) {
        i4 = (int)df;
      } else {
        i4 = MIN_int32_T;
      }
    } else if (df >= 2.147483648E+9) {
      i4 = MAX_int32_T;
    } else {
      i4 = 0;
    }

    j_fprintf((double)i4 * Tin);

    /* =============== */
    /*  Read in the phase-dependent Stokes parameters (doubles) from file and */
    /*  generate the coherency matrix */
    df = b_fopen();

    /*  Read Stokes parameters */
    i4 = S->size[0] * S->size[1];
    S->size[0] = nbins;
    S->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)S, i4, (int)sizeof(double));
    qEnd = nbins << 2;
    for (i4 = 0; i4 < qEnd; i4++) {
      S->data[i4] = 0.0;
    }

    for (ii = 1; ii <= nbins; ii++) {
      b_fread(df, tmp_data, tmp_size);
      for (i4 = 0; i4 < 4; i4++) {
        S->data[(ii + S->size[0] * i4) - 1] = tmp_data[i4];
      }
    }

    b_emxInit_real_T(&S0, 1);
    b_fclose(df);
    qEnd = S->size[0];
    i4 = S0->size[0];
    S0->size[0] = qEnd;
    emxEnsureCapacity((emxArray__common *)S0, i4, (int)sizeof(double));
    for (i4 = 0; i4 < qEnd; i4++) {
      S0->data[i4] = S->data[i4];
    }

    b_emxInit_real_T(&S1, 1);
    qEnd = S->size[0];
    i4 = S1->size[0];
    S1->size[0] = qEnd;
    emxEnsureCapacity((emxArray__common *)S1, i4, (int)sizeof(double));
    for (i4 = 0; i4 < qEnd; i4++) {
      S1->data[i4] = S->data[i4 + S->size[0]];
    }

    b_emxInit_real_T(&S2, 1);
    qEnd = S->size[0];
    i4 = S2->size[0];
    S2->size[0] = qEnd;
    emxEnsureCapacity((emxArray__common *)S2, i4, (int)sizeof(double));
    for (i4 = 0; i4 < qEnd; i4++) {
      S2->data[i4] = S->data[i4 + (S->size[0] << 1)];
    }

    b_emxInit_real_T(&S3, 1);
    qEnd = S->size[0];
    i4 = S3->size[0];
    S3->size[0] = qEnd;
    emxEnsureCapacity((emxArray__common *)S3, i4, (int)sizeof(double));
    for (i4 = 0; i4 < qEnd; i4++) {
      S3->data[i4] = S->data[i4 + S->size[0] * 3];
    }

    emxFree_real_T(&S);
    b_emxInit_real_T(&Jxx, 1);

    /*  Create Coherency matrix */
    i4 = Jxx->size[0];
    Jxx->size[0] = S0->size[0];
    emxEnsureCapacity((emxArray__common *)Jxx, i4, (int)sizeof(double));
    qEnd = S0->size[0];
    for (i4 = 0; i4 < qEnd; i4++) {
      Jxx->data[i4] = 0.5 * (S0->data[i4] + S1->data[i4]);
    }

    i4 = S0->size[0];
    emxEnsureCapacity((emxArray__common *)S0, i4, (int)sizeof(double));
    qEnd = S0->size[0];
    for (i4 = 0; i4 < qEnd; i4++) {
      S0->data[i4] = 0.5 * (S0->data[i4] - S1->data[i4]);
    }

    emxInit_creal_T(&Jxy, 1);
    i4 = Jxy->size[0];
    Jxy->size[0] = S2->size[0];
    emxEnsureCapacity((emxArray__common *)Jxy, i4, (int)sizeof(creal_T));
    qEnd = S2->size[0];
    for (i4 = 0; i4 < qEnd; i4++) {
      df = S2->data[i4] + 0.0 * S3->data[i4];
      x = S3->data[i4];
      Jxy->data[i4].re = 0.5 * df;
      Jxy->data[i4].im = 0.5 * x;
    }

    emxInit_creal_T(&Jyx, 1);
    i4 = Jyx->size[0];
    Jyx->size[0] = S2->size[0];
    emxEnsureCapacity((emxArray__common *)Jyx, i4, (int)sizeof(creal_T));
    qEnd = S2->size[0];
    for (i4 = 0; i4 < qEnd; i4++) {
      df = S2->data[i4] - 0.0 * S3->data[i4];
      x = 0.0 - S3->data[i4];
      Jyx->data[i4].re = 0.5 * df;
      Jyx->data[i4].im = 0.5 * x;
    }

    emxInit_creal_T(&b_Jxx, 1);
    i4 = b_Jxx->size[0];
    b_Jxx->size[0] = Jxx->size[0];
    emxEnsureCapacity((emxArray__common *)b_Jxx, i4, (int)sizeof(creal_T));
    qEnd = Jxx->size[0];
    for (i4 = 0; i4 < qEnd; i4++) {
      b_Jxx->data[i4].re = Jxx->data[i4];
      b_Jxx->data[i4].im = 0.0;
    }

    emxInit_creal_T(&b_S0, 1);
    pEnd = Jxx->size[0];
    q = Jyx->size[0];
    p = Jxy->size[0];
    i4 = b_S0->size[0];
    b_S0->size[0] = S0->size[0];
    emxEnsureCapacity((emxArray__common *)b_S0, i4, (int)sizeof(creal_T));
    qEnd = S0->size[0];
    for (i4 = 0; i4 < qEnd; i4++) {
      b_S0->data[i4].re = S0->data[i4];
      b_S0->data[i4].im = 0.0;
    }

    b_emxInit_creal_T(&J, 2);
    yk = S0->size[0];
    i4 = J->size[0] * J->size[1];
    J->size[0] = pEnd;
    J->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)J, i4, (int)sizeof(creal_T));
    for (i4 = 0; i4 < pEnd; i4++) {
      J->data[i4] = b_Jxx->data[i4];
    }

    emxFree_creal_T(&b_Jxx);
    for (i4 = 0; i4 < q; i4++) {
      J->data[i4 + J->size[0]] = Jyx->data[i4];
    }

    emxFree_creal_T(&Jyx);
    for (i4 = 0; i4 < p; i4++) {
      J->data[i4 + (J->size[0] << 1)] = Jxy->data[i4];
    }

    emxFree_creal_T(&Jxy);
    for (i4 = 0; i4 < yk; i4++) {
      J->data[i4 + J->size[0] * 3] = b_S0->data[i4];
    }

    emxFree_creal_T(&b_S0);

    /*  Vector of relative times */
    i5 = Nin - 1L;
    if (i5 > 2147483647L) {
      i5 = 2147483647L;
    } else {
      if (i5 < -2147483648L) {
        i5 = -2147483648L;
      }
    }

    yk = (int)i5;
    if (yk < 0) {
      n = 0;
    } else {
      n = yk + 1;
    }

    emxInit_int32_T(&y, 2);
    i4 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, i4, (int)sizeof(int));
    if (n > 0) {
      y->data[0] = 0;
      yk = 0;
      for (k = 2; k <= n; k++) {
        yk++;
        y->data[k - 1] = yk;
      }
    }

    emxInit_real_T(&trel, 2);
    i4 = trel->size[0] * trel->size[1];
    trel->size[0] = 1;
    trel->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)trel, i4, (int)sizeof(double));
    qEnd = y->size[0] * y->size[1];
    for (i4 = 0; i4 < qEnd; i4++) {
      trel->data[i4] = (double)y->data[i4] * Tin;
    }

    emxFree_int32_T(&y);

    /* =============== */
    /*  Open output file for writing */
    out_fid = c_fopen();

    /*  Write header */
    /* hdr = zeros(hdrsize,1); */
    /* fprintf(fid, '%c ', hdr); */
    df = rt_roundd_snf((double)Nin - n_lo);
    if (df < 2.147483648E+9) {
      if (df >= -2.147483648E+9) {
        i4 = (int)df;
      } else {
        i4 = MIN_int32_T;
      }
    } else if (df >= 2.147483648E+9) {
      i4 = MAX_int32_T;
    } else {
      i4 = 0;
    }

    if ((int)(1.0 + n_hi) > i4) {
      n = 0;
      i4 = 0;
    } else {
      n = (int)(1.0 + n_hi) - 1;
    }

    emxInit_real_T(&ttclip, 2);
    p = ttclip->size[0] * ttclip->size[1];
    ttclip->size[0] = 1;
    ttclip->size[1] = i4 - n;
    emxEnsureCapacity((emxArray__common *)ttclip, p, (int)sizeof(double));
    qEnd = i4 - n;
    for (i4 = 0; i4 < qEnd; i4++) {
      ttclip->data[ttclip->size[0] * i4] = trel->data[n + i4];
    }

    /*  just to initialise the type of ttclip */
    ii = 1;
    emxInit_real_T(&tt, 2);
    emxInit_creal32_T(&z, 2);
    b_emxInit_creal_T(&z0, 2);
    b_emxInit_creal_T(&zjj, 2);
    b_emxInit_creal32_T(&f1a, 1);
    b_emxInit_creal32_T(&f2a, 1);
    b_emxInit_creal32_T(&f1, 1);
    b_emxInit_creal32_T(&z1, 1);
    b_emxInit_creal32_T(&z2, 1);
    b_emxInit_real_T(&r1, 1);
    b_emxInit_real_T(&r2, 1);
    emxInit_real32_T(&r3, 1);
    emxInit_real32_T(&r4, 1);
    emxInit_real32_T(&r5, 1);
    emxInit_real32_T(&b_y, 1);
    b_emxInit_real32_T(&b_x, 2);
    emxInit_creal32_T(&c_x, 2);
    emxInit_boolean_T(&d_x, 1);
    b_emxInit_int32_T(&idx, 1);
    b_emxInit_int32_T(&idx0, 1);
    b_emxInit_creal32_T(&b_z, 1);
    b_emxInit_creal32_T(&c_z, 1);
    emxInit_creal_T(&b_S3, 1);
    emxInit_creal_T(&r6, 1);
    b_emxInit_creal_T(&c_S3, 2);
    while (ii <= nseries) {
      /*  Print loop number */
      l_fprintf(ii, nseries);

      /*  Time vector */
      if (ii == 1) {
        i4 = tt->size[0] * tt->size[1];
        tt->size[0] = 1;
        tt->size[1] = trel->size[1];
        emxEnsureCapacity((emxArray__common *)tt, i4, (int)sizeof(double));
        df = t0 - n_hi * Tin;
        qEnd = trel->size[0] * trel->size[1];
        for (i4 = 0; i4 < qEnd; i4++) {
          tt->data[i4] = df + trel->data[i4];
        }
      } else {
        i4 = tt->size[0] * tt->size[1];
        tt->size[0] = 1;
        tt->size[1] = trel->size[1];
        emxEnsureCapacity((emxArray__common *)tt, i4, (int)sizeof(double));
        df = ttclip->data[ttclip->size[1] - 1] - (n_hi - 1.0) * Tin;
        qEnd = trel->size[0] * trel->size[1];
        for (i4 = 0; i4 < qEnd; i4++) {
          tt->data[i4] = df + trel->data[i4];
        }
      }

      findphase(tt, nbins, period, S0);
      na = S0->size[0];
      n = S0->size[0];
      unnamed_idx_0 = (unsigned int)S0->size[0];
      i4 = idx->size[0];
      idx->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)idx, i4, (int)sizeof(int));
      if (S0->size[0] == 0) {
        for (k = 1; k <= n; k++) {
          idx->data[k - 1] = k;
        }
      } else {
        for (k = 1; k <= n; k++) {
          idx->data[k - 1] = k;
        }

        for (k = 1; k <= n - 1; k += 2) {
          if ((S0->data[k - 1] <= S0->data[k]) || rtIsNaN(S0->data[k])) {
            b_p = true;
          } else {
            b_p = false;
          }

          if (b_p) {
          } else {
            idx->data[k - 1] = k + 1;
            idx->data[k] = k;
          }
        }

        i4 = idx0->size[0];
        idx0->size[0] = S0->size[0];
        emxEnsureCapacity((emxArray__common *)idx0, i4, (int)sizeof(int));
        qEnd = S0->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          idx0->data[i4] = 1;
        }

        i = 2;
        while (i < n) {
          yk = i << 1;
          j = 1;
          for (pEnd = 1 + i; pEnd < n + 1; pEnd = qEnd + i) {
            p = j;
            q = pEnd - 1;
            qEnd = j + yk;
            if (qEnd > n + 1) {
              qEnd = n + 1;
            }

            k = 0;
            kEnd = qEnd - j;
            while (k + 1 <= kEnd) {
              if ((S0->data[idx->data[p - 1] - 1] <= S0->data[idx->data[q] - 1])
                  || rtIsNaN(S0->data[idx->data[q] - 1])) {
                b_p = true;
              } else {
                b_p = false;
              }

              if (b_p) {
                idx0->data[k] = idx->data[p - 1];
                p++;
                if (p == pEnd) {
                  while (q + 1 < qEnd) {
                    k++;
                    idx0->data[k] = idx->data[q];
                    q++;
                  }
                }
              } else {
                idx0->data[k] = idx->data[q];
                q++;
                if (q + 1 == qEnd) {
                  while (p < pEnd) {
                    k++;
                    idx0->data[k] = idx->data[p - 1];
                    p++;
                  }
                }
              }

              k++;
            }

            for (k = 0; k + 1 <= kEnd; k++) {
              idx->data[(j + k) - 1] = idx0->data[k];
            }

            j = qEnd;
          }

          i = yk;
        }
      }

      unnamed_idx_0 = (unsigned int)S0->size[0];
      i4 = S1->size[0];
      S1->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)S1, i4, (int)sizeof(double));
      for (k = 0; k + 1 <= na; k++) {
        S1->data[k] = S0->data[idx->data[k] - 1];
      }

      k = 0;
      while ((k + 1 <= na) && rtIsInf(S1->data[k]) && (S1->data[k] < 0.0)) {
        k++;
      }

      p = k;
      k = S0->size[0];
      while ((k >= 1) && rtIsNaN(S1->data[k - 1])) {
        k--;
      }

      q = S0->size[0] - k;
      while ((k >= 1) && rtIsInf(S1->data[k - 1]) && (S1->data[k - 1] > 0.0)) {
        k--;
      }

      pEnd = (S0->size[0] - k) - q;
      qEnd = -1;
      if (p > 0) {
        qEnd = 0;
      }

      yk = (p + k) - p;
      while (p + 1 <= yk) {
        x = S1->data[p];
        do {
          exitg4 = 0;
          p++;
          if (p + 1 > yk) {
            exitg4 = 1;
          } else {
            df = fabs(x / 2.0);
            if ((!rtIsInf(df)) && (!rtIsNaN(df))) {
              if (df <= 2.2250738585072014E-308) {
                df = 4.94065645841247E-324;
              } else {
                frexp(df, &exponent);
                df = ldexp(1.0, exponent - 53);
              }
            } else {
              df = rtNaN;
            }

            if ((fabs(x - S1->data[p]) < df) || (rtIsInf(S1->data[p]) && rtIsInf
                 (x) && ((S1->data[p] > 0.0) == (x > 0.0)))) {
              b_p = true;
            } else {
              b_p = false;
            }

            if (!b_p) {
              exitg4 = 1;
            }
          }
        } while (exitg4 == 0);

        qEnd++;
        S1->data[qEnd] = x;
      }

      if (pEnd > 0) {
        qEnd++;
        S1->data[qEnd] = S1->data[yk];
      }

      p = yk + pEnd;
      for (j = 1; j <= q; j++) {
        qEnd++;
        S1->data[qEnd] = S1->data[(p + j) - 1];
      }

      i4 = S1->size[0];
      if (1 > qEnd + 1) {
        i6 = -1;
      } else {
        i6 = qEnd;
      }

      S1->size[0] = i6 + 1;
      emxEnsureCapacity((emxArray__common *)S1, i4, (int)sizeof(double));

      /*  Initialize vector sizes */
      i4 = z->size[0] * z->size[1];
      z->size[0] = Nin;
      z->size[1] = 2;
      emxEnsureCapacity((emxArray__common *)z, i4, (int)sizeof(creal32_T));
      qEnd = Nin << 1;
      for (i4 = 0; i4 < qEnd; i4++) {
        z->data[i4].re = 0.0F;
        z->data[i4].im = 0.0F;
      }

      /*  need different sized arrays for standalone code to work in */
      /*  'complextoreal' and 'complextocomplex' cases */
      /*  Loop through groups of data that share the same phase. Random data  */
      /*  in each group are generated from the same coherency matrix */
      i4 = S1->size[0];
      for (jj = 0; jj + 1 <= i4; jj++) {
        /* Get coherency matrix for this pulsar phase */
        Jcoh[0] = J->data[(int)S1->data[jj] - 1];
        Jcoh[2] = J->data[((int)S1->data[jj] + (J->size[0] << 1)) - 1];
        Jcoh[1] = J->data[((int)S1->data[jj] + J->size[0]) - 1];
        Jcoh[3] = J->data[((int)S1->data[jj] + J->size[0] * 3) - 1];

        /*  Indices of elements with a given phase */
        df = S1->data[jj];
        n = d_x->size[0];
        d_x->size[0] = S0->size[0];
        emxEnsureCapacity((emxArray__common *)d_x, n, (int)sizeof(boolean_T));
        qEnd = S0->size[0];
        for (n = 0; n < qEnd; n++) {
          d_x->data[n] = (S0->data[n] == df);
        }

        p = d_x->size[0];
        pEnd = 0;
        n = idx->size[0];
        idx->size[0] = d_x->size[0];
        emxEnsureCapacity((emxArray__common *)idx, n, (int)sizeof(int));
        yk = 1;
        exitg3 = false;
        while ((!exitg3) && (yk <= p)) {
          guard1 = false;
          if (d_x->data[yk - 1]) {
            pEnd++;
            idx->data[pEnd - 1] = yk;
            if (pEnd >= p) {
              exitg3 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            yk++;
          }
        }

        if (d_x->size[0] == 1) {
          if (pEnd == 0) {
            n = idx->size[0];
            idx->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)idx, n, (int)sizeof(int));
          }
        } else {
          n = idx->size[0];
          if (1 > pEnd) {
            idx->size[0] = 0;
          } else {
            idx->size[0] = pEnd;
          }

          emxEnsureCapacity((emxArray__common *)idx, n, (int)sizeof(int));
        }

        n = S2->size[0];
        S2->size[0] = idx->size[0];
        emxEnsureCapacity((emxArray__common *)S2, n, (int)sizeof(double));
        qEnd = idx->size[0];
        for (n = 0; n < qEnd; n++) {
          S2->data[n] = idx->data[n];
        }

        /* Generate two randomly-phased, unit-length phasors   */
        /* z0 = exp(complex(0,1)*2*pi()*rand(nL,npol)); */
        randn(S2->size[0], S3);
        randn(S2->size[0], Jxx);
        randn(S2->size[0], r1);
        randn(S2->size[0], r2);
        n = b_S3->size[0];
        b_S3->size[0] = S3->size[0];
        emxEnsureCapacity((emxArray__common *)b_S3, n, (int)sizeof(creal_T));
        qEnd = S3->size[0];
        for (n = 0; n < qEnd; n++) {
          b_S3->data[n].re = S3->data[n];
          b_S3->data[n].im = Jxx->data[n];
        }

        yk = S3->size[0];
        n = r6->size[0];
        r6->size[0] = r1->size[0];
        emxEnsureCapacity((emxArray__common *)r6, n, (int)sizeof(creal_T));
        qEnd = r1->size[0];
        for (n = 0; n < qEnd; n++) {
          r6->data[n].re = r1->data[n];
          r6->data[n].im = r2->data[n];
        }

        qEnd = r1->size[0];
        n = c_S3->size[0] * c_S3->size[1];
        c_S3->size[0] = yk;
        c_S3->size[1] = 2;
        emxEnsureCapacity((emxArray__common *)c_S3, n, (int)sizeof(creal_T));
        for (n = 0; n < yk; n++) {
          c_S3->data[n] = b_S3->data[n];
        }

        for (n = 0; n < qEnd; n++) {
          c_S3->data[n + c_S3->size[0]] = r6->data[n];
        }

        n = z0->size[0] * z0->size[1];
        z0->size[0] = c_S3->size[0];
        z0->size[1] = 2;
        emxEnsureCapacity((emxArray__common *)z0, n, (int)sizeof(creal_T));
        for (n = 0; n < 2; n++) {
          qEnd = c_S3->size[0];
          for (p = 0; p < qEnd; p++) {
            z0->data[p + z0->size[0] * n].re = 0.70710678118654757 * c_S3->
              data[p + c_S3->size[0] * n].re;
            z0->data[p + z0->size[0] * n].im = 0.70710678118654757 * c_S3->
              data[p + c_S3->size[0] * n].im;
          }
        }

        /* Generate covariant vectors via Cholesky decomposition */
        j = 0;
        do {
          exitg1 = 0;
          if (j + 1 < 3) {
            if (Jcoh[j + (j << 1)].im != 0.0) {
              exitg1 = 1;
            } else {
              j++;
            }
          } else {
            yk = 0;
            q = 0;
            j = 1;
            exitg2 = false;
            while ((!exitg2) && (j < 3)) {
              pEnd = (q + j) - 1;
              c_re = 0.0;
              if (j - 1 < 1) {
              } else {
                c_re = Jcoh[q].re * Jcoh[q].re + Jcoh[q].im * Jcoh[q].im;
              }

              df = Jcoh[pEnd].re - c_re;
              if ((Jcoh[pEnd].im == 0.0) && (df > 0.0)) {
                df = sqrt(df);
                Jcoh[pEnd].re = df;
                Jcoh[pEnd].im = 0.0;
                if (j < 2) {
                  c_re = 1.0 / df;
                  for (k = pEnd + 2; k + 1 <= pEnd + 3; k += 2) {
                    df = Jcoh[k].re;
                    x = Jcoh[k].im;
                    Jcoh[k].re = c_re * Jcoh[k].re - 0.0 * Jcoh[k].im;
                    Jcoh[k].im = c_re * x + 0.0 * df;
                  }

                  q += 2;
                }

                j++;
              } else {
                Jcoh[pEnd].re = df;
                Jcoh[pEnd].im = 0.0;
                yk = j;
                exitg2 = true;
              }
            }

            if (yk == 0) {
              yk = 2;
            } else {
              yk--;
            }

            for (j = 0; j + 1 <= yk; j++) {
              i = j + 2;
              while (i <= yk) {
                Jcoh[1 + (j << 1)].re = 0.0;
                Jcoh[1 + (j << 1)].im = 0.0;
                i = 3;
              }
            }

            exitg1 = 1;
          }
        } while (exitg1 == 0);

        unnamed_idx_0 = (unsigned int)z0->size[0];
        na = z0->size[0];
        n = zjj->size[0] * zjj->size[1];
        zjj->size[0] = (int)unnamed_idx_0;
        zjj->size[1] = 2;
        emxEnsureCapacity((emxArray__common *)zjj, n, (int)sizeof(creal_T));
        qEnd = (int)unnamed_idx_0 << 1;
        for (n = 0; n < qEnd; n++) {
          zjj->data[n].re = 0.0;
          zjj->data[n].im = 0.0;
        }

        if (z0->size[0] == 0) {
        } else {
          yk = 0;
          while ((na > 0) && (yk <= na)) {
            n = yk + na;
            for (kEnd = yk; kEnd + 1 <= n; kEnd++) {
              zjj->data[kEnd].re = 0.0;
              zjj->data[kEnd].im = 0.0;
            }

            yk += na;
          }

          pEnd = 0;
          yk = 0;
          while ((na > 0) && (yk <= na)) {
            p = 0;
            for (q = pEnd; q + 1 <= pEnd + 2; q++) {
              if ((Jcoh[q].re != 0.0) || (Jcoh[q].im != 0.0)) {
                c_re = Jcoh[q].re - 0.0 * Jcoh[q].im;
                df = Jcoh[q].im + 0.0 * Jcoh[q].re;
                qEnd = p;
                n = yk + na;
                for (kEnd = yk; kEnd + 1 <= n; kEnd++) {
                  qEnd++;
                  x = c_re * z0->data[qEnd - 1].re - df * z0->data[qEnd - 1].im;
                  c_im = c_re * z0->data[qEnd - 1].im + df * z0->data[qEnd - 1].
                    re;
                  zjj->data[kEnd].re += x;
                  zjj->data[kEnd].im += c_im;
                }
              }

              p += na;
            }

            pEnd += 2;
            yk += na;
          }
        }

        /* z = transpose(chol(Jcoh, 'lower')*transpose(z0)); %alternative */
        /*  Concatenate with data from other phases */
        for (n = 0; n < 2; n++) {
          qEnd = zjj->size[0];
          for (p = 0; p < qEnd; p++) {
            f1a_re = (float)zjj->data[p + zjj->size[0] * n].re;
            f1a_im = (float)zjj->data[p + zjj->size[0] * n].im;
            z->data[((int)S2->data[p] + z->size[0] * n) - 1].re = f1a_re;
            z->data[((int)S2->data[p] + z->size[0] * n) - 1].im = f1a_im;
          }
        }

        /* iL = iL + nL; % increment to next starting index in z */
      }

      /*  Forward FFT */
      qEnd = z->size[0];
      i4 = c_z->size[0];
      c_z->size[0] = qEnd;
      emxEnsureCapacity((emxArray__common *)c_z, i4, (int)sizeof(creal32_T));
      for (i4 = 0; i4 < qEnd; i4++) {
        c_z->data[i4] = z->data[i4];
      }

      fft(c_z, Nin, f1a);
      qEnd = z->size[0];
      i4 = b_z->size[0];
      b_z->size[0] = qEnd;
      emxEnsureCapacity((emxArray__common *)b_z, i4, (int)sizeof(creal32_T));
      for (i4 = 0; i4 < qEnd; i4++) {
        b_z->data[i4] = z->data[i4 + z->size[0]];
      }

      fft(b_z, Nin, f2a);

      /*  Element-wise multiplication by dispersion matrix. */
      i4 = f1a->size[0];
      emxEnsureCapacity((emxArray__common *)f1a, i4, (int)sizeof(creal32_T));
      qEnd = f1a->size[0];
      for (i4 = 0; i4 < qEnd; i4++) {
        H_re = (float)H->data[i4].re;
        H_im = (float)H->data[i4].im;
        f1a_re = f1a->data[i4].re;
        f1a_im = f1a->data[i4].im;
        f1a->data[i4].re = f1a_re * H_re - f1a_im * H_im;
        f1a->data[i4].im = f1a_re * H_im + f1a_im * H_re;
      }

      i4 = f2a->size[0];
      emxEnsureCapacity((emxArray__common *)f2a, i4, (int)sizeof(creal32_T));
      qEnd = f2a->size[0];
      for (i4 = 0; i4 < qEnd; i4++) {
        H_re = (float)H->data[i4].re;
        H_im = (float)H->data[i4].im;
        f1a_re = f2a->data[i4].re;
        f1a_im = f2a->data[i4].im;
        f2a->data[i4].re = f1a_re * H_re - f1a_im * H_im;
        f2a->data[i4].im = f1a_re * H_im + f1a_im * H_re;
      }

      /*  If complextoreal, then create a Hermitian array */
      switch ((int)eml_switch_helper(dformat_data, dformat_size)) {
       case 0:
        /* Create Hermitian vector */
        if (2 > Nin) {
          i4 = 1;
          n = 3;
        } else {
          i4 = 2;
          n = Nin + 3;
        }

        if (2 > Nin) {
          p = 0;
          pEnd = 0;
        } else {
          p = 1;
          pEnd = Nin;
        }

        yk = z2->size[0];
        z2->size[0] = pEnd - p;
        emxEnsureCapacity((emxArray__common *)z2, yk, (int)sizeof(creal32_T));
        qEnd = pEnd - p;
        for (pEnd = 0; pEnd < qEnd; pEnd++) {
          z2->data[pEnd].re = f1a->data[p + pEnd].re;
          z2->data[pEnd].im = -f1a->data[p + pEnd].im;
        }

        flipud(z2);
        p = f1->size[0];
        f1->size[0] = (n - i4) + z2->size[0];
        emxEnsureCapacity((emxArray__common *)f1, p, (int)sizeof(creal32_T));
        f1a_re = f1a->data[0].re;
        f1->data[0].re = f1a_re;
        f1->data[0].im = 0.0F;
        qEnd = n - i4;
        for (p = 0; p <= qEnd - 3; p++) {
          f1->data[p + 1].re = (float)Pmul * f1a->data[(i4 + p) - 1].re;
          f1->data[p + 1].im = (float)Pmul * f1a->data[(i4 + p) - 1].im;
        }

        f1a_im = f1a->data[0].im;
        f1->data[(n - i4) - 1].re = f1a_im;
        f1->data[(n - i4) - 1].im = 0.0F;
        qEnd = z2->size[0];
        for (p = 0; p < qEnd; p++) {
          f1->data[(p + n) - i4].re = (float)Pmul * z2->data[p].re;
          f1->data[(p + n) - i4].im = (float)Pmul * z2->data[p].im;
        }

        if (2 > Nin) {
          i4 = 1;
          n = 3;
        } else {
          i4 = 2;
          n = Nin + 3;
        }

        if (2 > Nin) {
          p = 0;
          pEnd = 0;
        } else {
          p = 1;
          pEnd = Nin;
        }

        yk = z2->size[0];
        z2->size[0] = pEnd - p;
        emxEnsureCapacity((emxArray__common *)z2, yk, (int)sizeof(creal32_T));
        qEnd = pEnd - p;
        for (pEnd = 0; pEnd < qEnd; pEnd++) {
          z2->data[pEnd].re = f2a->data[p + pEnd].re;
          z2->data[pEnd].im = -f2a->data[p + pEnd].im;
        }

        flipud(z2);
        p = f1a->size[0];
        f1a->size[0] = (n - i4) + z2->size[0];
        emxEnsureCapacity((emxArray__common *)f1a, p, (int)sizeof(creal32_T));
        f1a_re = f2a->data[0].re;
        f1a->data[0].re = f1a_re;
        f1a->data[0].im = 0.0F;
        qEnd = n - i4;
        for (p = 0; p <= qEnd - 3; p++) {
          f1a->data[p + 1].re = (float)Pmul * f2a->data[(i4 + p) - 1].re;
          f1a->data[p + 1].im = (float)Pmul * f2a->data[(i4 + p) - 1].im;
        }

        f1a_im = f2a->data[0].im;
        f1a->data[(n - i4) - 1].re = f1a_im;
        f1a->data[(n - i4) - 1].im = 0.0F;
        qEnd = z2->size[0];
        for (p = 0; p < qEnd; p++) {
          f1a->data[(p + n) - i4].re = (float)Pmul * z2->data[p].re;
          f1a->data[(p + n) - i4].im = (float)Pmul * z2->data[p].im;
        }
        break;

       default:
        i4 = f1->size[0];
        f1->size[0] = f1a->size[0];
        emxEnsureCapacity((emxArray__common *)f1, i4, (int)sizeof(creal32_T));
        qEnd = f1a->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          f1->data[i4] = f1a->data[i4];
        }

        i4 = f1a->size[0];
        f1a->size[0] = f2a->size[0];
        emxEnsureCapacity((emxArray__common *)f1a, i4, (int)sizeof(creal32_T));
        qEnd = f2a->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          f1a->data[i4] = f2a->data[i4];
        }
        break;
      }

      /*  Inverse FFT */
      /*  Optionally include an fftshift before the inverse FFT, as needed */
      if (shift == 1) {
        fftshift(f1);
        fftshift(f1a);
      }

      ifft(f1, Nout, z1);
      ifft(f1a, Nout, z2);

      /*  Remove convolution overlap region */
      df = rt_roundd_snf((double)Nin - n_lo);
      if (df < 2.147483648E+9) {
        if (df >= -2.147483648E+9) {
          i4 = (int)df;
        } else {
          i4 = MIN_int32_T;
        }
      } else if (df >= 2.147483648E+9) {
        i4 = MAX_int32_T;
      } else {
        i4 = 0;
      }

      if ((int)(1.0 + n_hi) > i4) {
        n = 0;
        i4 = 0;
      } else {
        n = (int)(1.0 + n_hi) - 1;
      }

      p = ttclip->size[0] * ttclip->size[1];
      ttclip->size[0] = 1;
      ttclip->size[1] = i4 - n;
      emxEnsureCapacity((emxArray__common *)ttclip, p, (int)sizeof(double));
      qEnd = i4 - n;
      for (i4 = 0; i4 < qEnd; i4++) {
        ttclip->data[ttclip->size[0] * i4] = tt->data[n + i4];
      }

      i4 = (int)(1.0 + n_hi * (double)Nmul);
      df = rt_roundd_snf((double)Nout - n_lo * (double)Nmul);
      if (df < 2.147483648E+9) {
        if (df >= -2.147483648E+9) {
          n = (int)df;
        } else {
          n = MIN_int32_T;
        }
      } else if (df >= 2.147483648E+9) {
        n = MAX_int32_T;
      } else {
        n = 0;
      }

      if (i4 > n) {
        i4 = 1;
        n = 1;
      } else {
        n++;
      }

      p = f1a->size[0];
      f1a->size[0] = n - i4;
      emxEnsureCapacity((emxArray__common *)f1a, p, (int)sizeof(creal32_T));
      qEnd = n - i4;
      for (p = 0; p < qEnd; p++) {
        f1a->data[p] = z1->data[(i4 + p) - 1];
      }

      p = (int)(1.0 + n_hi * (double)Nmul);
      df = rt_roundd_snf((double)Nout - n_lo * (double)Nmul);
      if (df < 2.147483648E+9) {
        if (df >= -2.147483648E+9) {
          pEnd = (int)df;
        } else {
          pEnd = MIN_int32_T;
        }
      } else if (df >= 2.147483648E+9) {
        pEnd = MAX_int32_T;
      } else {
        pEnd = 0;
      }

      if (p > pEnd) {
        p = 1;
        pEnd = 1;
      } else {
        pEnd++;
      }

      yk = f2a->size[0];
      f2a->size[0] = pEnd - p;
      emxEnsureCapacity((emxArray__common *)f2a, yk, (int)sizeof(creal32_T));
      qEnd = pEnd - p;
      for (yk = 0; yk < qEnd; yk++) {
        f2a->data[yk] = z2->data[(p + yk) - 1];
      }

      /*  Interleave polarizations into a single vector */
      if (eml_strcmp(dformat_data, dformat_size)) {
        if (nclip_out > 1073741823) {
          kEnd = MAX_int32_T;
        } else if (nclip_out <= -1073741824) {
          kEnd = MIN_int32_T;
        } else {
          kEnd = nclip_out << 1;
        }

        qEnd = n - i4;
        i4 = c_x->size[0] * c_x->size[1];
        c_x->size[0] = 2;
        c_x->size[1] = qEnd;
        emxEnsureCapacity((emxArray__common *)c_x, i4, (int)sizeof(creal32_T));
        for (i4 = 0; i4 < qEnd; i4++) {
          c_x->data[c_x->size[0] * i4] = f1a->data[i4];
        }

        qEnd = pEnd - p;
        for (i4 = 0; i4 < qEnd; i4++) {
          c_x->data[1 + c_x->size[0] * i4] = f2a->data[i4];
        }

        p = c_x->size[1] << 1;
        i4 = f1a->size[0];
        f1a->size[0] = kEnd;
        emxEnsureCapacity((emxArray__common *)f1a, i4, (int)sizeof(creal32_T));
        for (k = 0; k + 1 <= p; k++) {
          f1a->data[k] = c_x->data[k];
        }
      } else {
        /*  complextocomplex */
        i4 = b_y->size[0];
        b_y->size[0] = f1a->size[0];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int)sizeof(float));
        qEnd = f1a->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          b_y->data[i4] = f1a->data[i4].re;
        }

        i4 = r3->size[0];
        r3->size[0] = f1a->size[0];
        emxEnsureCapacity((emxArray__common *)r3, i4, (int)sizeof(float));
        qEnd = f1a->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          r3->data[i4] = f1a->data[i4].im;
        }

        i4 = r4->size[0];
        r4->size[0] = f2a->size[0];
        emxEnsureCapacity((emxArray__common *)r4, i4, (int)sizeof(float));
        qEnd = f2a->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          r4->data[i4] = f2a->data[i4].re;
        }

        i4 = r5->size[0];
        r5->size[0] = f2a->size[0];
        emxEnsureCapacity((emxArray__common *)r5, i4, (int)sizeof(float));
        qEnd = f2a->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          r5->data[i4] = f2a->data[i4].im;
        }

        if (nclip_out > 536870911) {
          kEnd = MAX_int32_T;
        } else if (nclip_out <= -536870912) {
          kEnd = MIN_int32_T;
        } else {
          kEnd = nclip_out << 2;
        }

        yk = b_y->size[0];
        qEnd = r3->size[0];
        pEnd = r4->size[0];
        q = r5->size[0];
        i4 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 4;
        b_x->size[1] = yk;
        emxEnsureCapacity((emxArray__common *)b_x, i4, (int)sizeof(float));
        for (i4 = 0; i4 < yk; i4++) {
          b_x->data[b_x->size[0] * i4] = b_y->data[i4];
        }

        for (i4 = 0; i4 < qEnd; i4++) {
          b_x->data[1 + b_x->size[0] * i4] = r3->data[i4];
        }

        for (i4 = 0; i4 < pEnd; i4++) {
          b_x->data[2 + b_x->size[0] * i4] = r4->data[i4];
        }

        for (i4 = 0; i4 < q; i4++) {
          b_x->data[3 + b_x->size[0] * i4] = r5->data[i4];
        }

        p = b_x->size[1] << 2;
        i4 = b_y->size[0];
        b_y->size[0] = kEnd;
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int)sizeof(float));
        for (k = 0; k + 1 <= p; k++) {
          b_y->data[k] = b_x->data[k];
        }

        i4 = f1a->size[0];
        f1a->size[0] = b_y->size[0];
        emxEnsureCapacity((emxArray__common *)f1a, i4, (int)sizeof(creal32_T));
        qEnd = b_y->size[0];
        for (i4 = 0; i4 < qEnd; i4++) {
          f1a->data[i4].re = b_y->data[i4];
          f1a->data[i4].im = 0.0F;
        }
      }

      /* Write vector to file */
      /*    - normally use binary file, unless generating stand-alone code */
      /*      generation, where fwite() is unsupported */
      if (out_type == 1) {
      } else {
        /*    - use ASCII text format for stand-alone code generation */
        /*      (convert to binary file later) */
        /*    - do 16 values each call to make it faster */
        i4 = (int)rt_roundd_snf((double)kEnd / 16.0);
        for (i = 1; i <= i4; i++) {
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          b_f1a = f1a->data[n - 16];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          c_f1a = f1a->data[n - 15];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          d_f1a = f1a->data[n - 14];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          e_f1a = f1a->data[n - 13];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          f_f1a = f1a->data[n - 12];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          g_f1a = f1a->data[n - 11];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          h_f1a = f1a->data[n - 10];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          i_f1a = f1a->data[n - 9];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          j_f1a = f1a->data[n - 8];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          k_f1a = f1a->data[n - 7];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          l_f1a = f1a->data[n - 6];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          m_f1a = f1a->data[n - 5];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          n_f1a = f1a->data[n - 4];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          o_f1a = f1a->data[n - 3];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          p_f1a = f1a->data[n - 2];
          if (i > 134217727) {
            n = MAX_int32_T;
          } else {
            n = i << 4;
          }

          q_f1a = f1a->data[n - 1];
          n_fprintf(out_fid, b_f1a, c_f1a, d_f1a, e_f1a, f_f1a, g_f1a, h_f1a,
                    i_f1a, j_f1a, k_f1a, l_f1a, m_f1a, n_f1a, o_f1a, p_f1a,
                    q_f1a);
        }
      }

      ii++;
    }

    emxFree_creal_T(&c_S3);
    emxFree_creal_T(&r6);
    emxFree_creal_T(&b_S3);
    emxFree_creal32_T(&c_z);
    emxFree_creal32_T(&b_z);
    emxFree_int32_T(&idx0);
    emxFree_int32_T(&idx);
    emxFree_boolean_T(&d_x);
    emxFree_creal32_T(&c_x);
    emxFree_real32_T(&b_x);
    emxFree_real32_T(&b_y);
    emxFree_real32_T(&r5);
    emxFree_real32_T(&r4);
    emxFree_real32_T(&r3);
    emxFree_real_T(&r2);
    emxFree_real_T(&r1);
    emxFree_creal_T(&H);
    emxFree_creal32_T(&z2);
    emxFree_creal32_T(&z1);
    emxFree_creal32_T(&f1);
    emxFree_creal32_T(&f2a);
    emxFree_creal32_T(&f1a);
    emxFree_creal_T(&zjj);
    emxFree_creal_T(&z0);
    emxFree_creal32_T(&z);
    emxFree_real_T(&ttclip);
    emxFree_real_T(&tt);
    emxFree_real_T(&trel);
    emxFree_creal_T(&J);
    emxFree_real_T(&Jxx);
    emxFree_real_T(&S3);
    emxFree_real_T(&S2);
    emxFree_real_T(&S1);
    emxFree_real_T(&S0);
    b_fclose(out_fid);
  }
}

/* End of code generation (test_vec_gen.c) */
