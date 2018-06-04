#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <vmblock.h>
#include <bikub.h>

/*
   Test program:  Bicubic splines 2   and
                  compute the value of a bicubic
                  splines at a point.

   (testing subroutines  bikub2 and bsval here.)
*/

int main (int argc, char *argv[])
{
  #define IMAX 20
  int  m, n, i, j, k, error, num;
  REAL x [IMAX], y [IMAX], value, tx, ty;
  mat4x4 **b;
  int    fehler;
  void   *vmblock;

  if ((fehler = umleiten(argc, argv))    /* assign input/output files */
      != 0)                              /* to standard ones          */
    return fehler;  /* 1 or 2 */

  printf ("Bicubic spline 2\n");

  do {
    printf ("number of x and y intervals : \n");
    scanf  ("%d", &m);
    scanf  ("%d", &n);
    if (m < 0 || n < 0 || m > IMAX || n > IMAX)
      printf ("No values in excess of %d please !\n", IMAX);
  } while (m < 0 || n < 0 || m > IMAX || n > IMAX);

  printf ("end points for x intervals:\n");
  for (i=0; i<=m; i++)
    scanf ("%"LZS"f",&x[i]);

  printf ("end points for y intervals:\n");
  for (i=0; i<=n; i++)
    scanf ("%"LZS"f",&y[i]);

  vmblock = vminit();
  b = (mat4x4 **)vmalloc(vmblock, MMATRIX, m + 1, n + 1);
  if (! vmcomplete(vmblock))   /* allocations partially unsuccessful? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  printf ("Function values b[i][j][0][0] (i=0,...,m, j=0,...,n):\n");
  for (i=0; i<=m; i++)
    for (j=0; j<=n; j++)
      scanf ("%"LZS"f", &b[i][j][0][0]);

  printf ("Input b[i][j][0][0]:\n");
  for (i=0; i<=m; i++)
  {
    for (j=0; j<=n; j++) printf ("%10.5"LZP"f ", b [i][j][0][0]);
    printf ("\n");
  }

  error = bikub2 (m,n,b,x,y);

  if (! error)
  {
    printf ("\nOutput:\n");
    for (i=0; i<=m-2; i++)
    {
      printf ("\n\n");
      for (j=0; j<=n-2; j++)
      {
        printf ("\n");
        for (k=0; k<=3; k++)
           printf ("%10.7"LZP"f  %10.7"LZP"f  %10.7"LZP"f  "
                   "%10.7"LZP"f\n",
                   b[i][j][k][0], b[i][j][k][1],
                   b[i][j][k][2], b[i][j][k][3]);
      }
    }

    printf ("Some test values:\n");

    for (num = 0; num <= 10; num++)
    {
      error = bsval (m, n, b, x, y,
                     tx = x[0] + (x[m]-x[0])/TEN * num,
                     ty = y[0] + (y[n]-y[0])/TEN * num,
                     &value);
      printf ("\n (%"LZP"g,%"LZP"g) => bsval = %"LZP"g (error = %d)\n",
              tx,ty, value,error);
    }
  }
  else
    printf ("error code  = %d\n", error);

  vmfree(vmblock);

  return 0;
}
