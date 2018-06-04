#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <glsp.h>

/* Test program for the function glspnp. */

int main (void)
{
  int  i, mode, error, n=8;
  REAL x[11] = {1.0,2.0,3.0,3.5,4.8,5.7,7.0,8.5,9.2};
  REAL y[11] = {4.0,5.0,4.5,3.2,2.7,2.0,1.0,1.5,2.5};
  REAL w[11] = {1000.0,500.0,100.0,10.0,1.0,100.0,10.0,500.0,1000.0};
  REAL alpha[3] = {1., 0.,-1.2},
       beta [3] = {1.5,0.,-1.7};
  REAL a[11], b[11], c[11], d[11];

  printf (" i       x          y            w\n");
  printf ("---------------------------------------\n");
  for (i=0; i<=n; ++i)
    printf ("%2d  %10.6"LZP"f  %10.6"LZP"f  %11.6"LZP"f\n",
            i, x[i], y[i], w[i]);

  for (mode = 0; mode <= 2; mode ++)
  {
    error = glspnp (n, x,y,w, mode+1, alpha[mode],beta[mode], a,b,c,d);
    if (error != 0)
      printf ("error %d has occured in glspnp \n", error);
    else
    {
      printf("\nSpline coefficients for glspnp mode = %d\n\n", mode);
      printf(" i       a           b            c            d\n");
      printf("-----------------------------------------------------\n");
      for (i=0; i<=n-1; ++i)
        printf ("%2d  %10.6"LZP"f  %10.6"LZP"f  %10.6"LZP"f  "
                "%10.6"LZP"f\n",
                i,   a[i],    b[i],    c[i],    d[i]      );
    }
  }

  return 0;
}
