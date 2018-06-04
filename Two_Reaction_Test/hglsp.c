#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <glsp.h>

/* Test program for the functions glspnp and glsppe. */

int main (void)
{
  int i, error, n=10;
  REAL   x[11] = {0.0,1.5,2.5,5.5,8.0,10.0,12.0,14.5,17.5,18.5,20.0};
  REAL   y[11] = {0.0,2.0,0.0,-4.0,-2.0,0.0,2.0,4.0,0.0,-2.0,0.0};
  REAL   w[11] = {10.,7.5,5.,10.,5.,1.,5.,10.,5.,7.5,10.};
  REAL   a[11], b[11], c[11], d[11];
  REAL   hup[86], h[11], h1[11], h2[11], h3[11], rs[11];

  error = glspnp (n, x, y, w, 4, 0.0, 0.0, a, b, c, d);

  if (error != 0)
    printf ("error %d has occured in  glspnp \n",error);
  else
  {
    printf (" i       x          y            w\n");
    printf ("---------------------------------------\n");
    for (i=0; i<=n; ++i)
      printf ("%2d  %10.6"LZP"f  %10.6"LZP"f  %10.6"LZP"f\n",
              i, x[i], y[i], w[i]);
    printf ("\n\nSpline coefficients from glspnp:\n\n");
    printf (" i        a               b               c               "
            "d\n");
    printf ("----------------------------------------------------------"
            "---------\n");
    for (i=0; i<=n-1; ++i)
      printf ("%2d  %14.6"LZP"e  %14.6"LZP"e  %14.6"LZP"e  "
              "%14.6"LZP"e\n",
              i,   a[i],    b[i],    c[i],    d[i]      );
  }
 /*
     test glsppe with  rep = 1
 */
  error = glsppe (n, x, y, w, 0, a, b, c, d, h, h1, h2, h3, rs, hup);
  error = glsppe (n, x, y, w, 1, a, b, c, d, h, h1, h2, h3, rs, hup);

  if (error != 0)
    printf ("error %d has occured\n",error);
  else
  {
    printf ("\n\nSpline coefficients from glsppe:\n\n");
    printf (" i        a               b               c               "
            "d\n");
    printf ("----------------------------------------------------------"
            "---------\n");
    for (i=0; i<=n-1; ++i)
      printf ("%2d  %14.6"LZP"e  %14.6"LZP"e  %14.6"LZP"e  "
              "%14.6"LZP"e\n",
              i,   a[i],    b[i],    c[i],    d[i]      );
  }

  return 0;
}
