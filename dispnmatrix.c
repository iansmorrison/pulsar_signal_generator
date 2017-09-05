/*
 * dispnmatrix.c
 *
 * Code generation for function 'dispnmatrix'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fprintf.h"
#include "test_vec_gen_emxutil.h"
#include "rdivide.h"
#include "power.h"
#include "test_vec_gen_rtwutil.h"
#include <stdio.h>

/* Function Definitions */
void b_dispnmatrix(const double frange[2], double usable_BW, int Nout, double DD,
                   double Tout, emxArray_creal_T *H, double *f, double *n_hi,
                   double *n_lo)
{
  double bandwidth;
  double deltaf;
  long i7;
  int yk;
  int n;
  emxArray_int32_T *y;
  int i8;
  int k;
  emxArray_real_T *b_fabs;
  emxArray_real_T *fc;
  emxArray_real_T *b_fc;
  emxArray_real_T *b;
  emxArray_real_T *b_f;
  emxArray_real_T *r7;
  emxArray_real_T *b_DD;
  emxArray_creal_T *b_y;
  emxArray_real_T *r8;
  double fcmin_idx_1;
  double mtmp;
  boolean_T exitg2;
  boolean_T exitg1;

  /*  Calculates the de-dispersion matrix and the number of leading */
  /*  (n_hi) and trailing (n_lo) elements that need to be removed during */
  /*  overlap-save. */
  /*  */
  /*  Inputs: */
  /*  -------- */
  /*    frange    - 2-element vector containing highest and lowest freqs */
  /*    usable_BW - fraction of frange to be used (e.g. oversampled pass-band) */
  /*    Nout      - length of FFT to be analysed */
  /*    nfreq     - number of frequency channels */
  /*    DD        - DM*Dconst  */
  /*    Tout      - sampling period of data */
  /*    direc     - Dispersion = 1; De-dispersion = -1 (default) */
  /*  */
  /*  Outputs: */
  /*  -------- */
  /*  */
  /*    H       - de-dispersion matrix */
  /*    f       - mean frequency of each frequency channel */
  /*    n_hi    - number of leading elements to be removed */
  /*    n_lo    - number of trailing elements to be removed */
  /*  */
  /*  Description: */
  /*  ------------ */
  /*  Calculates the de-dispersion matrix used to multiply the Fourier  */
  /*  transformed, filterbanked data. Also returns the number of leading and */
  /*  trailing elements that need to be removed. */
  /*   */
  /*  Changes: */
  /*  -------- */
  /*  */
  /*  Author           Date         Comments */
  /*  ---------------  -----------  ---------------------------------------- */
  /*  D. Hicks         21-Apr-2014  Original version */
  /*  I. Morrison      31-Jul-2015  Added "usable_BW" for oversampled case */
  /*  ---------------------------------------------------------------------- */
  /*  if ~exist('direc', 'var'), */
  /*      direc = -1; */
  /*  end; */
  /*  Vector of frequency bin assignments. Note that the highest frequency */
  /*  bin is assigned a value of frange(2)-df where frange(2) is the Nyquist */
  /*  frequency */
  bandwidth = frange[1] - frange[0];
  deltaf = bandwidth / (double)Nout;

  /*  frequency spacing */
  i7 = Nout - 1L;
  if (i7 > 2147483647L) {
    i7 = 2147483647L;
  } else {
    if (i7 < -2147483648L) {
      i7 = -2147483648L;
    }
  }

  yk = (int)i7;
  if (yk < 0) {
    n = 0;
  } else {
    n = yk + 1;
  }

  emxInit_int32_T(&y, 2);
  i8 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity((emxArray__common *)y, i8, (int)sizeof(int));
  if (n > 0) {
    y->data[0] = 0;
    yk = 0;
    for (k = 2; k <= n; k++) {
      yk++;
      y->data[k - 1] = yk;
    }
  }

  emxInit_real_T(&b_fabs, 2);
  i8 = b_fabs->size[0] * b_fabs->size[1];
  b_fabs->size[0] = 1;
  b_fabs->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b_fabs, i8, (int)sizeof(double));
  yk = y->size[0] * y->size[1];
  for (i8 = 0; i8 < yk; i8++) {
    b_fabs->data[i8] = frange[0] + (double)y->data[i8] * deltaf;
  }

  b_emxInit_real_T(&fc, 1);

  /* fabs = linspace(frange(1), frange(2), Nout*nfreq ); %Not quite correct */
  /*  Absolute freqs in each channel */
  i8 = fc->size[0];
  fc->size[0] = Nout;
  emxEnsureCapacity((emxArray__common *)fc, i8, (int)sizeof(double));
  for (k = 0; k + 1 <= b_fabs->size[1]; k++) {
    fc->data[k] = b_fabs->data[k];
  }

  emxFree_real_T(&b_fabs);

  /*  Mean of each freq channel */
  if (fc->size[0] == 0) {
    deltaf = 0.0;
  } else {
    deltaf = fc->data[0];
    for (k = 2; k <= fc->size[0]; k++) {
      deltaf += fc->data[k - 1];
    }
  }

  b_emxInit_real_T(&b_fc, 1);
  *f = deltaf / (double)fc->size[0];

  /*  Expand to create matrix  */
  /*  Phase dispersion */
  i8 = b_fc->size[0];
  b_fc->size[0] = fc->size[0];
  emxEnsureCapacity((emxArray__common *)b_fc, i8, (int)sizeof(double));
  yk = fc->size[0];
  for (i8 = 0; i8 < yk; i8++) {
    b_fc->data[i8] = fc->data[i8] - *f;
  }

  b_emxInit_real_T(&b, 1);
  b_emxInit_real_T(&b_f, 1);
  b_power(b_fc, b);
  i8 = b_f->size[0];
  b_f->size[0] = Nout;
  emxEnsureCapacity((emxArray__common *)b_f, i8, (int)sizeof(double));
  emxFree_real_T(&b_fc);
  for (i8 = 0; i8 < Nout; i8++) {
    b_f->data[i8] = *f;
  }

  b_emxInit_real_T(&r7, 1);
  b_emxInit_real_T(&b_DD, 1);
  b_power(b_f, r7);
  i8 = b_DD->size[0];
  b_DD->size[0] = b->size[0];
  emxEnsureCapacity((emxArray__common *)b_DD, i8, (int)sizeof(double));
  yk = b->size[0];
  emxFree_real_T(&b_f);
  for (i8 = 0; i8 < yk; i8++) {
    b_DD->data[i8] = DD * b->data[i8];
  }

  emxInit_creal_T(&b_y, 1);
  b_emxInit_real_T(&r8, 1);
  c_rdivide(b_DD, r7, r8);
  c_rdivide(r8, fc, b);

  /* Dispersion matrix: De-dispersion has direc = -1; Dispersion has direc = 1 */
  i8 = b_y->size[0];
  b_y->size[0] = b->size[0];
  emxEnsureCapacity((emxArray__common *)b_y, i8, (int)sizeof(creal_T));
  yk = b->size[0];
  emxFree_real_T(&b_DD);
  emxFree_real_T(&r8);
  emxFree_real_T(&r7);
  for (i8 = 0; i8 < yk; i8++) {
    b_y->data[i8].re = b->data[i8] * 1.0E+6 * 0.0;
    b_y->data[i8].im = b->data[i8] * 1.0E+6 * 6.2831853071795862;
  }

  emxFree_real_T(&b);
  i8 = H->size[0];
  H->size[0] = b_y->size[0];
  emxEnsureCapacity((emxArray__common *)H, i8, (int)sizeof(creal_T));
  yk = b_y->size[0];
  for (i8 = 0; i8 < yk; i8++) {
    H->data[i8] = b_y->data[i8];
  }

  for (k = 0; k < b_y->size[0]; k++) {
    if (rtIsInf(H->data[k].im) && rtIsInf(H->data[k].re) && (H->data[k].re < 0.0))
    {
      fcmin_idx_1 = 0.0;
      deltaf = 0.0;
    } else {
      deltaf = exp(H->data[k].re / 2.0);
      fcmin_idx_1 = deltaf * (deltaf * cos(H->data[k].im));
      deltaf *= deltaf * sin(H->data[k].im);
    }

    H->data[k].re = fcmin_idx_1;
    H->data[k].im = deltaf;
  }

  emxFree_creal_T(&b_y);

  /*  Mask to zero those frequency components outside the usable bandwidth. */
  /*  When the de-dispersion kernel is applied to oversampled data, it will */
  /*  cause the transition band (and any spectral leakage it contains) to be */
  /*  discarded. */
  if (fabs(usable_BW) < fabs(bandwidth)) {
    d_fprintf();
    deltaf = (fabs(bandwidth) - fabs(usable_BW)) / fabs(bandwidth) / 2.0;
    fcmin_idx_1 = rt_roundd_snf(deltaf * (double)Nout);
    if (fcmin_idx_1 < 2.147483648E+9) {
      if (fcmin_idx_1 >= -2.147483648E+9) {
        yk = (int)fcmin_idx_1;
      } else {
        yk = MIN_int32_T;
      }
    } else if (fcmin_idx_1 >= 2.147483648E+9) {
      yk = MAX_int32_T;
    } else {
      yk = 0;
    }

    if (1 > yk) {
      yk = 0;
    }

    i8 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = yk;
    emxEnsureCapacity((emxArray__common *)y, i8, (int)sizeof(int));
    for (i8 = 0; i8 < yk; i8++) {
      y->data[y->size[0] * i8] = i8;
    }

    yk = y->size[0] * y->size[1];
    for (i8 = 0; i8 < yk; i8++) {
      H->data[y->data[i8]].re = 0.0;
      H->data[y->data[i8]].im = 0.0;
    }

    fcmin_idx_1 = rt_roundd_snf((1.0 - deltaf) * (double)Nout);
    if (fcmin_idx_1 < 2.147483648E+9) {
      if (fcmin_idx_1 >= -2.147483648E+9) {
        yk = (int)fcmin_idx_1;
      } else {
        yk = MIN_int32_T;
      }
    } else if (fcmin_idx_1 >= 2.147483648E+9) {
      yk = MAX_int32_T;
    } else {
      yk = 0;
    }

    if (yk > Nout) {
      i8 = 0;
      k = 0;
    } else {
      i8 = yk - 1;
      k = Nout;
    }

    yk = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = k - i8;
    emxEnsureCapacity((emxArray__common *)y, yk, (int)sizeof(int));
    yk = k - i8;
    for (k = 0; k < yk; k++) {
      y->data[y->size[0] * k] = i8 + k;
    }

    yk = y->size[0] * y->size[1];
    for (i8 = 0; i8 < yk; i8++) {
      H->data[y->data[i8]].re = 0.0;
      H->data[y->data[i8]].im = 0.0;
    }
  }

  emxFree_int32_T(&y);

  /*  Calculate convolution overlap region */
  /* fcmin = [fabs(nnmin(1)) fabs(nnmax(1))]; */
  deltaf = fc->data[0];
  fcmin_idx_1 = fc->data[Nout - 1];

  /* CHECK THIS !!!!!! */
  bandwidth = (deltaf + fcmin_idx_1) / 2.0;

  /*  CHECK WHETHER n_hi and n_lo ARE FLIPPED DEPENDING ON THE DIRECTION */
  /*  OF THE DISPERSION (i.e. DISPERSING OR DE-DISPERSING) */
  yk = 1;
  mtmp = deltaf;
  emxFree_real_T(&fc);
  if (rtIsNaN(deltaf)) {
    k = 2;
    exitg2 = false;
    while ((!exitg2) && (k < 3)) {
      yk = 2;
      if (!rtIsNaN(fcmin_idx_1)) {
        mtmp = fcmin_idx_1;
        exitg2 = true;
      } else {
        k = 3;
      }
    }
  }

  if ((yk < 2) && (fcmin_idx_1 > mtmp)) {
    mtmp = fcmin_idx_1;
  }

  *n_hi = ceil(DD * (1.0 / (bandwidth * bandwidth) - 1.0 / (mtmp * mtmp)) / Tout);
  yk = 1;
  if (rtIsNaN(deltaf)) {
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 3)) {
      yk = 2;
      if (!rtIsNaN(fcmin_idx_1)) {
        deltaf = fcmin_idx_1;
        exitg1 = true;
      } else {
        k = 3;
      }
    }
  }

  if ((yk < 2) && (fcmin_idx_1 < deltaf)) {
    deltaf = fcmin_idx_1;
  }

  *n_lo = ceil(DD * (1.0 / (deltaf * deltaf) - 1.0 / (bandwidth * bandwidth)) /
               Tout);

  /*  If overlap exceeds length of vector, force there to be no overlap */
  /*  (useful for debugging with short time series). */
  if (*n_hi + *n_lo >= Nout) {
    *n_hi = 0.0;
    *n_lo = 0.0;
    f_fprintf();
  }
}

void dispnmatrix(const double frange[2], double usable_BW, int Nout, int nfreq,
                 double DD, double Tout, double direc, emxArray_creal_T *H,
                 emxArray_real_T *f, double *n_hi, double *n_lo)
{
  double bandwidth;
  long i0;
  double deltaf;
  int yk;
  int n;
  emxArray_int32_T *y;
  int i;
  int k;
  emxArray_real_T *b_fabs;
  emxArray_real_T *fc;
  double sz[2];
  int ix;
  emxArray_real_T *f0c;
  emxArray_real_T *b_fc;
  emxArray_real_T *b;
  emxArray_real_T *r0;
  emxArray_real_T *b_DD;
  double mtmp;
  boolean_T exitg2;
  boolean_T exitg1;

  /*  Calculates the de-dispersion matrix and the number of leading */
  /*  (n_hi) and trailing (n_lo) elements that need to be removed during */
  /*  overlap-save. */
  /*  */
  /*  Inputs: */
  /*  -------- */
  /*    frange    - 2-element vector containing highest and lowest freqs */
  /*    usable_BW - fraction of frange to be used (e.g. oversampled pass-band) */
  /*    Nout      - length of FFT to be analysed */
  /*    nfreq     - number of frequency channels */
  /*    DD        - DM*Dconst  */
  /*    Tout      - sampling period of data */
  /*    direc     - Dispersion = 1; De-dispersion = -1 (default) */
  /*  */
  /*  Outputs: */
  /*  -------- */
  /*  */
  /*    H       - de-dispersion matrix */
  /*    f       - mean frequency of each frequency channel */
  /*    n_hi    - number of leading elements to be removed */
  /*    n_lo    - number of trailing elements to be removed */
  /*  */
  /*  Description: */
  /*  ------------ */
  /*  Calculates the de-dispersion matrix used to multiply the Fourier  */
  /*  transformed, filterbanked data. Also returns the number of leading and */
  /*  trailing elements that need to be removed. */
  /*   */
  /*  Changes: */
  /*  -------- */
  /*  */
  /*  Author           Date         Comments */
  /*  ---------------  -----------  ---------------------------------------- */
  /*  D. Hicks         21-Apr-2014  Original version */
  /*  I. Morrison      31-Jul-2015  Added "usable_BW" for oversampled case */
  /*  ---------------------------------------------------------------------- */
  /*  if ~exist('direc', 'var'), */
  /*      direc = -1; */
  /*  end; */
  if ((direc != -1.0) && (direc != 1.0)) {
    direc = -1.0;
    b_fprintf();
  }

  /*  Vector of frequency bin assignments. Note that the highest frequency */
  /*  bin is assigned a value of frange(2)-df where frange(2) is the Nyquist */
  /*  frequency */
  bandwidth = frange[1] - frange[0];
  i0 = (long)Nout * nfreq;
  if (i0 > 2147483647L) {
    i0 = 2147483647L;
  } else {
    if (i0 < -2147483648L) {
      i0 = -2147483648L;
    }
  }

  deltaf = bandwidth / (double)(int)i0;

  /*  frequency spacing */
  i0 = (long)Nout * nfreq;
  if (i0 > 2147483647L) {
    i0 = 2147483647L;
  } else {
    if (i0 < -2147483648L) {
      i0 = -2147483648L;
    }
  }

  i0 = (int)i0 - 1L;
  if (i0 > 2147483647L) {
    i0 = 2147483647L;
  } else {
    if (i0 < -2147483648L) {
      i0 = -2147483648L;
    }
  }

  yk = (int)i0;
  if (yk < 0) {
    n = 0;
  } else {
    n = yk + 1;
  }

  emxInit_int32_T(&y, 2);
  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(int));
  if (n > 0) {
    y->data[0] = 0;
    yk = 0;
    for (k = 2; k <= n; k++) {
      yk++;
      y->data[k - 1] = yk;
    }
  }

  emxInit_real_T(&b_fabs, 2);
  i = b_fabs->size[0] * b_fabs->size[1];
  b_fabs->size[0] = 1;
  b_fabs->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b_fabs, i, (int)sizeof(double));
  yk = y->size[0] * y->size[1];
  for (i = 0; i < yk; i++) {
    b_fabs->data[i] = frange[0] + (double)y->data[i] * deltaf;
  }

  emxInit_real_T(&fc, 2);

  /* fabs = linspace(frange(1), frange(2), Nout*nfreq ); %Not quite correct */
  /*  Absolute freqs in each channel */
  i = fc->size[0] * fc->size[1];
  fc->size[0] = Nout;
  fc->size[1] = nfreq;
  emxEnsureCapacity((emxArray__common *)fc, i, (int)sizeof(double));
  for (k = 0; k + 1 <= b_fabs->size[1]; k++) {
    fc->data[k] = b_fabs->data[k];
  }

  /*  Mean of each freq channel */
  for (i = 0; i < 2; i++) {
    sz[i] = fc->size[i];
  }

  i = b_fabs->size[0] * b_fabs->size[1];
  b_fabs->size[0] = 1;
  b_fabs->size[1] = (int)sz[1];
  emxEnsureCapacity((emxArray__common *)b_fabs, i, (int)sizeof(double));
  if ((fc->size[0] == 0) || (fc->size[1] == 0)) {
    i = b_fabs->size[0] * b_fabs->size[1];
    b_fabs->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)b_fabs, i, (int)sizeof(double));
    i = b_fabs->size[0] * b_fabs->size[1];
    b_fabs->size[1] = (int)sz[1];
    emxEnsureCapacity((emxArray__common *)b_fabs, i, (int)sizeof(double));
    yk = (int)sz[1];
    for (i = 0; i < yk; i++) {
      b_fabs->data[i] = 0.0;
    }
  } else {
    ix = 0;
    n = -1;
    for (i = 1; i <= fc->size[1]; i++) {
      yk = ix;
      ix++;
      deltaf = fc->data[yk];
      for (k = 2; k <= fc->size[0]; k++) {
        ix++;
        deltaf += fc->data[ix - 1];
      }

      n++;
      b_fabs->data[n] = deltaf;
    }
  }

  emxInit_real_T(&f0c, 2);
  rdivide(b_fabs, fc->size[0], f);

  /*  Expand to create matrix  */
  i = f0c->size[0] * f0c->size[1];
  f0c->size[0] = Nout;
  f0c->size[1] = f->size[1];
  emxEnsureCapacity((emxArray__common *)f0c, i, (int)sizeof(double));
  yk = f->size[1];
  emxFree_real_T(&b_fabs);
  if ((Nout == 0) || (yk == 0)) {
  } else {
    for (yk = 0; yk + 1 <= f->size[1]; yk++) {
      n = yk * Nout;
      for (ix = 1; ix <= Nout; ix++) {
        f0c->data[(n + ix) - 1] = f->data[yk];
      }
    }
  }

  emxInit_real_T(&b_fc, 2);

  /*  Phase dispersion */
  i = b_fc->size[0] * b_fc->size[1];
  b_fc->size[0] = fc->size[0];
  b_fc->size[1] = fc->size[1];
  emxEnsureCapacity((emxArray__common *)b_fc, i, (int)sizeof(double));
  yk = fc->size[0] * fc->size[1];
  for (i = 0; i < yk; i++) {
    b_fc->data[i] = fc->data[i] - f0c->data[i];
  }

  emxInit_real_T(&b, 2);
  emxInit_real_T(&r0, 2);
  emxInit_real_T(&b_DD, 2);
  power(b_fc, b);
  power(f0c, r0);
  i = b_DD->size[0] * b_DD->size[1];
  b_DD->size[0] = b->size[0];
  b_DD->size[1] = b->size[1];
  emxEnsureCapacity((emxArray__common *)b_DD, i, (int)sizeof(double));
  yk = b->size[0] * b->size[1];
  emxFree_real_T(&b_fc);
  for (i = 0; i < yk; i++) {
    b_DD->data[i] = DD * b->data[i];
  }

  b_rdivide(b_DD, r0, b);
  b_rdivide(b, fc, f0c);

  /* Dispersion matrix: De-dispersion has direc = -1; Dispersion has direc = 1 */
  deltaf = 3.1415926535897931 * (2.0 * direc);
  i = H->size[0] * H->size[1];
  H->size[0] = f0c->size[0];
  H->size[1] = f0c->size[1];
  emxEnsureCapacity((emxArray__common *)H, i, (int)sizeof(creal_T));
  yk = f0c->size[0] * f0c->size[1];
  emxFree_real_T(&b_DD);
  emxFree_real_T(&r0);
  emxFree_real_T(&b);
  for (i = 0; i < yk; i++) {
    H->data[i].re = f0c->data[i] * 1.0E+6 * 0.0;
    H->data[i].im = f0c->data[i] * 1.0E+6 * deltaf;
  }

  emxFree_real_T(&f0c);
  i = H->size[0] * H->size[1];
  for (k = 0; k < i; k++) {
    if (rtIsInf(H->data[k].im) && rtIsInf(H->data[k].re) && (H->data[k].re < 0.0))
    {
      mtmp = 0.0;
      deltaf = 0.0;
    } else {
      deltaf = exp(H->data[k].re / 2.0);
      mtmp = deltaf * (deltaf * cos(H->data[k].im));
      deltaf *= deltaf * sin(H->data[k].im);
    }

    H->data[k].re = mtmp;
    H->data[k].im = deltaf;
  }

  /*  Mask to zero those frequency components outside the usable bandwidth. */
  /*  When the de-dispersion kernel is applied to oversampled data, it will */
  /*  cause the transition band (and any spectral leakage it contains) to be */
  /*  discarded. */
  if (fabs(usable_BW) < fabs(bandwidth)) {
    d_fprintf();
    deltaf = (fabs(bandwidth) - fabs(usable_BW)) / fabs(bandwidth) / 2.0;
    mtmp = rt_roundd_snf(deltaf * (double)Nout);
    if (mtmp < 2.147483648E+9) {
      if (mtmp >= -2.147483648E+9) {
        yk = (int)mtmp;
      } else {
        yk = MIN_int32_T;
      }
    } else if (mtmp >= 2.147483648E+9) {
      yk = MAX_int32_T;
    } else {
      yk = 0;
    }

    if (1 > yk) {
      yk = 0;
    }

    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = yk;
    emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(int));
    for (i = 0; i < yk; i++) {
      y->data[y->size[0] * i] = i;
    }

    yk = y->size[0] * y->size[1];
    for (i = 0; i < yk; i++) {
      H->data[y->data[i]].re = 0.0;
      H->data[y->data[i]].im = 0.0;
    }

    mtmp = rt_roundd_snf((1.0 - deltaf) * (double)Nout);
    if (mtmp < 2.147483648E+9) {
      if (mtmp >= -2.147483648E+9) {
        yk = (int)mtmp;
      } else {
        yk = MIN_int32_T;
      }
    } else if (mtmp >= 2.147483648E+9) {
      yk = MAX_int32_T;
    } else {
      yk = 0;
    }

    if (yk > Nout) {
      i = 0;
      ix = 0;
    } else {
      i = yk - 1;
      ix = Nout;
    }

    yk = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = ix - i;
    emxEnsureCapacity((emxArray__common *)y, yk, (int)sizeof(int));
    yk = ix - i;
    for (ix = 0; ix < yk; ix++) {
      y->data[y->size[0] * ix] = i + ix;
    }

    yk = y->size[0] * y->size[1];
    for (i = 0; i < yk; i++) {
      H->data[y->data[i]].re = 0.0;
      H->data[y->data[i]].im = 0.0;
    }
  }

  emxFree_int32_T(&y);

  /*  Calculate convolution overlap region */
  /* fcmin = [fabs(nnmin(1)) fabs(nnmax(1))]; */
  sz[0] = fc->data[0];
  sz[1] = fc->data[Nout - 1];

  /* CHECK THIS !!!!!! */
  deltaf = (sz[0] + sz[1]) / 2.0;

  /*  CHECK WHETHER n_hi and n_lo ARE FLIPPED DEPENDING ON THE DIRECTION */
  /*  OF THE DISPERSION (i.e. DISPERSING OR DE-DISPERSING) */
  yk = 1;
  mtmp = sz[0];
  emxFree_real_T(&fc);
  if (rtIsNaN(sz[0])) {
    ix = 2;
    exitg2 = false;
    while ((!exitg2) && (ix < 3)) {
      yk = 2;
      if (!rtIsNaN(sz[1])) {
        mtmp = sz[1];
        exitg2 = true;
      } else {
        ix = 3;
      }
    }
  }

  if ((yk < 2) && (sz[1] > mtmp)) {
    mtmp = sz[1];
  }

  *n_hi = ceil(DD * (1.0 / (deltaf * deltaf) - 1.0 / (mtmp * mtmp)) / Tout);
  yk = 1;
  mtmp = sz[0];
  if (rtIsNaN(sz[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 3)) {
      yk = 2;
      if (!rtIsNaN(sz[1])) {
        mtmp = sz[1];
        exitg1 = true;
      } else {
        ix = 3;
      }
    }
  }

  if ((yk < 2) && (sz[1] < mtmp)) {
    mtmp = sz[1];
  }

  *n_lo = ceil(DD * (1.0 / (mtmp * mtmp) - 1.0 / (deltaf * deltaf)) / Tout);

  /*  If overlap exceeds length of vector, force there to be no overlap */
  /*  (useful for debugging with short time series). */
  if (*n_hi + *n_lo >= Nout) {
    *n_hi = 0.0;
    *n_lo = 0.0;
    f_fprintf();
  }
}

/* End of code generation (dispnmatrix.c) */
