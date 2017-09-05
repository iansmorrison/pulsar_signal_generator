/*
 * fft.c
 *
 * Code generation for function 'fft'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "dispnmatrix.h"
#include "findphase.h"
#include "test_vec_gen.h"
#include "fft.h"
#include "test_vec_gen_emxutil.h"
#include "ifft.h"
#include "test_vec_gen_rtwutil.h"
#include <stdio.h>

/* Function Definitions */
void fft(const emxArray_creal32_T *x, int varargin_1, emxArray_creal32_T *y)
{
  int nd2;
  int u1;
  int ixDelta;
  emxArray_real32_T *costab1q;
  int nRowsD2;
  int nRowsD4;
  int lastChan;
  float e;
  int k;
  emxArray_real32_T *costab;
  emxArray_real32_T *sintab;
  int n;
  int n2;
  int ix;
  int chanStart;
  int i;
  boolean_T tst;
  float temp_re;
  float temp_im;
  int iDelta2;
  int iheight;
  int ihi;
  float twid_im;
  nd2 = y->size[0];
  y->size[0] = varargin_1;
  emxEnsureCapacity((emxArray__common *)y, nd2, (int)sizeof(creal32_T));
  if (varargin_1 > x->size[0]) {
    nd2 = y->size[0];
    y->size[0] = varargin_1;
    emxEnsureCapacity((emxArray__common *)y, nd2, (int)sizeof(creal32_T));
    for (nd2 = 0; nd2 < varargin_1; nd2++) {
      y->data[nd2].re = 0.0F;
      y->data[nd2].im = 0.0F;
    }
  }

  if (x->size[0] == 0) {
  } else {
    u1 = x->size[0];
    if (varargin_1 <= u1) {
      u1 = varargin_1;
    }

    nd2 = (x->size[0] - u1) + 1;
    if (1 >= nd2) {
      ixDelta = 1;
    } else {
      ixDelta = nd2;
    }

    b_emxInit_real32_T(&costab1q, 2);
    nRowsD2 = varargin_1 / 2;
    nRowsD4 = nRowsD2 / 2;
    lastChan = varargin_1 * (div_s32(x->size[0], x->size[0]) - 1);
    e = 6.28318548F / (float)varargin_1;
    nd2 = costab1q->size[0] * costab1q->size[1];
    costab1q->size[0] = 1;
    costab1q->size[1] = nRowsD4 + 1;
    emxEnsureCapacity((emxArray__common *)costab1q, nd2, (int)sizeof(float));
    costab1q->data[0] = 1.0F;
    nd2 = nRowsD4 / 2;
    for (k = 1; k <= nd2; k++) {
      costab1q->data[k] = (real32_T)cos(e * (float)k);
    }

    for (k = nd2 + 1; k < nRowsD4; k++) {
      costab1q->data[k] = (real32_T)sin(e * (float)(nRowsD4 - k));
    }

    b_emxInit_real32_T(&costab, 2);
    b_emxInit_real32_T(&sintab, 2);
    costab1q->data[nRowsD4] = 0.0F;
    n = costab1q->size[1] - 1;
    n2 = (costab1q->size[1] - 1) << 1;
    nd2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = n2 + 1;
    emxEnsureCapacity((emxArray__common *)costab, nd2, (int)sizeof(float));
    nd2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = n2 + 1;
    emxEnsureCapacity((emxArray__common *)sintab, nd2, (int)sizeof(float));
    costab->data[0] = 1.0F;
    sintab->data[0] = 0.0F;
    for (k = 1; k <= n; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[n - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - n];
    }

    emxFree_real32_T(&costab1q);
    ix = 0;
    chanStart = 0;
    while ((varargin_1 > 0) && (chanStart <= lastChan)) {
      n2 = 0;
      nd2 = chanStart;
      for (i = 1; i < u1; i++) {
        y->data[nd2] = x->data[ix];
        n = varargin_1;
        tst = true;
        while (tst) {
          n >>= 1;
          n2 ^= n;
          tst = ((n2 & n) == 0);
        }

        nd2 = chanStart + n2;
        ix++;
      }

      y->data[nd2] = x->data[ix];
      ix += ixDelta;
      nd2 = (chanStart + varargin_1) - 2;
      if (varargin_1 > 1) {
        for (i = chanStart; i <= nd2; i += 2) {
          temp_re = y->data[i + 1].re;
          temp_im = y->data[i + 1].im;
          y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
          y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
        }
      }

      n = 2;
      iDelta2 = 4;
      k = nRowsD4;
      iheight = 1 + ((nRowsD4 - 1) << 2);
      while (k > 0) {
        i = chanStart;
        ihi = chanStart + iheight;
        while (i < ihi) {
          nd2 = i + n;
          temp_re = y->data[nd2].re;
          temp_im = y->data[nd2].im;
          y->data[i + n].re = y->data[i].re - y->data[nd2].re;
          y->data[i + n].im = y->data[i].im - y->data[nd2].im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
          i += iDelta2;
        }

        nd2 = chanStart + 1;
        for (n2 = k; n2 < nRowsD2; n2 += k) {
          e = costab->data[n2];
          twid_im = sintab->data[n2];
          i = nd2;
          ihi = nd2 + iheight;
          while (i < ihi) {
            temp_re = e * y->data[i + n].re - twid_im * y->data[i + n].im;
            temp_im = e * y->data[i + n].im + twid_im * y->data[i + n].re;
            y->data[i + n].re = y->data[i].re - temp_re;
            y->data[i + n].im = y->data[i].im - temp_im;
            y->data[i].re += temp_re;
            y->data[i].im += temp_im;
            i += iDelta2;
          }

          nd2++;
        }

        k /= 2;
        n = iDelta2;
        iDelta2 <<= 1;
        iheight -= n;
      }

      chanStart += varargin_1;
    }

    emxFree_real32_T(&sintab);
    emxFree_real32_T(&costab);
  }
}

/* End of code generation (fft.c) */
