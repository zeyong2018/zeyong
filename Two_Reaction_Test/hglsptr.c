#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <glsp.h>

/* Test program for the function glsptr. */

int main (void)
{
  int  i, n, marke_v, error;
  REAL xn[13], fn[13], w[13], a[13], b[13], c[13], d[13];
  REAL phin[13], r[13], help[161];
  REAL px = ZERO, py = ZERO, phid = ZERO;

  n = 12;
  marke_v = 2;

  xn[0] = 32.0;  xn[1] = 16.0;  xn[2] = -6.0;  xn[3] = -18.0;
  xn[4] = -20.0;  xn[5] = -20.0;  xn[6] = -13.0;  xn[7] = -4.0;
  xn[8] = 4.0;  xn[9] = 12.0;  xn[10] = 21.0;  xn[11] = 33.0;
  xn[12] = 32.0;

  fn[0] = 15.0;  fn[1] = 32.0;  fn[2] = 31.0;  fn[3] = 18.0;
  fn[4] = 10.0;  fn[5] = -7.0;  fn[6] = -22.0;  fn[7] = -26.0;
  fn[8] = -31.0; fn[9] = -28.0; fn[10] = -25.0; fn[11] = -11.0;
  fn[12] = 15.0;

  w[0] = 10.0;  w[1] = 5.0;  w[2] = 7.0;  w[3] = 10.0;
  w[4] = 1.0;  w[5] = 1.0;  w[6] = 2.0;  w[7] = 5.0;
  w[8] = 1.0;  w[9] = 5.0;  w[10] = 5.0;  w[11] = 10.0;  w[12] = 10.0;

  error =  glsptr (n, xn, fn, w, marke_v, &px, &py, a, b, c, d, phin, r,
                   &phid, help);

  printf ("error: %d\n", error);
  printf ("\n");
  for (i=0; i<n; i++)
    printf ("a[%2d]=%10.6"LZP"f   b[%2d]=%10.6"LZP"f   "
            "c[%2d]=%10.6"LZP"f   d[%2d]=%10.6"LZP"f\n",
             i, a[i], i, b[i], i, c[i], i, d[i]);

  return 0;
}
