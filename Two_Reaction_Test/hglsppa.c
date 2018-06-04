#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <glsp.h>

/* Test program for the function glsppa. */

int main (void)
{
  int  i, n, error;
  int  marke_t = 0, marke_w = 2, rand;
  REAL xn[13], fn[13], wx[13], wf[13], t[13];
  REAL ax[13], bx[13], cx[13], dx[13];
  REAL ay[13], by[13], cy[13], dy[13];
  REAL alpha[3][4], beta[3][4], help[159];
  REAL dummy1[3];
  REAL dummy2[3];

  for (rand = 1; rand <= 4; rand++)
  {
    if (rand == 4)                                        /* periodic */
    {
      n = 12;
      xn[ 0] = 32.0;  xn[1] =  25.0;  xn[ 2] =  16.0;  xn[ 3] =   5.0;
      xn[ 4] = -6.0;  xn[5] = -18.0;  xn[ 6] = -16.0;  xn[ 7] = -13.0;
      xn[ 8] = -6.0;  xn[9] =   4.0;  xn[10] =  21.0;  xn[11] =  33.0;
      xn[12] = xn[0];

      fn[ 0] =  1.5;  fn[1] =  2.5;  fn[ 2] =  3.2;  fn[ 3] =  3.0;
      fn[ 4] =  3.1;  fn[5] =  1.8;  fn[ 6] = -0.7;  fn[ 7] = -2.2;
      fn[ 8] = -3.4;  fn[9] = -3.1;  fn[10] = -2.5;  fn[11] = -1.1;
      fn[12] = fn[0];

      wx[ 0] = 2.00;  wx[1] = 0.20;  wx[ 2] = 0.15;  wx[ 3] = 0.15;
      wx[ 4] = 0.40;  wx[5] = 0.10;  wx[ 6] = 0.50;  wx[ 7] = 0.50;
      wx[ 8] = 0.50;  wx[9] = 0.50;  wx[10] = 0.50;  wx[11] = 0.50;
      wx[12] = wx[0];

      wf[ 0] = 0.60; wf[1] = 0.08; wf[ 2] = 0.01; wf[ 3] = 0.01;
      wf[ 4] = 0.10; wf[5] = 0.20; wf[ 6] = 0.10; wf[ 7] = 0.10;
      wf[ 8] = 0.05; wf[9] = 0.10; wf[10] = 0.10; wf[11] = 0.10;
      wf[12] = wf[0];

      error =  glsppa (n, xn, fn, wx, wf, t, marke_t, rand,
                       &dummy1[0], &dummy2[0], marke_w,
                       ax, bx, cx, dx, ay, by, cy, dy, help);
    }
    else
    {
      n = 8;
      alpha[1][1] =  0.1;  alpha[1][2] =  1.0;
      alpha[2][1] =  0.0;  alpha[2][2] =  0.0;
      alpha[3][1] =  1.5;

      beta [1][1] =  1.0;  beta [1][2] = -1.0;
      beta [2][1] =  0.0;  beta [2][2] =  0.0;
      beta [3][1] = -0.6;

      xn[0] = 1.0;  xn[1] = 1.5;  xn[2] = 2.0;  xn[3] = 2.5;
      xn[4] = 2.5;  xn[5] = 2.0;  xn[6] = 2.5;  xn[7] = 3.0;
      xn[8] = 4.0;

      fn[0] = 1.0;  fn[1] = 2.0;  fn[2] = 2.5;  fn[3] = 2.0;
      fn[4] = 1.0;  fn[5] = 1.5;  fn[6] = 2.0;  fn[7] = 3.0;
      fn[8] = 3.0;

      wx[0] =  10.0; wx[1] =  50.0; wx[2] = 100.0;  wx[3] =  80.0;
      wx[4] = 100.0; wx[5] = 300.0; wx[6] = 200.0;  wx[7] = 100.0;
      wx[8] =  10.0;

      wf[0] =  10.0; wf[1] =  50.0; wf[2] = 200.0; wf[3] =  80.0;
      wf[4] = 300.0; wf[5] = 100.0; wf[6] = 200.0; wf[7] = 100.0;
      wf[8] =  10.0;

      error =  glsppa (n, xn, fn, wx, wf, t, marke_t, rand,
                       &alpha[rand][0], &beta[rand][0],
                       marke_w, ax, bx, cx, dx, ay, by, cy,
                       dy, help);
    }

    printf ("\n%seriodic, ", rand==4 ? "P" : "Not p");
    printf ("end point condition %d, ", rand);
    printf ("error: %d\n\n", error);
    if (!error)
    {
      for (i=0; i<=n; i++)
        printf ("t[%d]=%"LZP"f\n", i, t[i]);
      printf ("\n");
      for (i=0; i<n; i++)
        printf ("ax[%2d]=%10.6"LZP"f   bx[%2d]=%10.6"LZP"f   "
                "cx[%2d]=%10.6"LZP"f   dx[%2d]=%10.6"LZP"f\n",
                 i, ax[i], i, bx[i], i, cx[i], i, dx[i]);
      printf ("\n");
      for (i=0; i<n; i++)
        printf ("ay[%2d]=%10.6"LZP"f   by[%2d]=%10.6"LZP"f   "
                "cy[%2d]=%10.6"LZP"f   dy[%2d]=%10.6"LZP"f\n",
                 i, ay[i], i, by[i], i, cy[i], i, dy[i]);
    }
  }

  return 0;
}
