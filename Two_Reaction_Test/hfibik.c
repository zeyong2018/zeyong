#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <vmblock.h>       /*  for  mat4x4, vminit, vmalloc, MMATRIX, */
                           /*       PMATRIX, vmcomplete               */
#include <bikub.h>
#include <kubatur.h>

/*
   Test program: Compute the surface integral over the domain
                 of a spline and compute a sub-rectangle
*/

int main (int argc, char *argv[])
{
  int  n, m, i, j, error;
  REAL x [20], y [20], hilf, xlow, ylow, xup, yup;
  mat4x4 **b;
  int    fehler;
  void   *vmblock;

  if ((fehler = umleiten(argc, argv))    /* assign input/output files */
      != 0)                              /* to standard ones          */
    return fehler;  /* 1 or 2 */

  printf ("Bicubic splines: surface integral:\n\n");

  printf ("number of  x- and y-intervals\n");
  scanf  ("%d %d", &n,&m);
  printf ("%d %d\n",  n, m);

  printf ("\nInput of %d x-nodes:\n", n+1);
  for (i=0; i<=n; i++)
    scanf  ("%"LZS"f",&x[i]);
  for (i=0; i<=n; i++)
    printf ("%"LZP"f\n", x[i]);

  printf ("\nInput of %d y-nodes:\n", m+1);
  for (i=0; i<=m; i++)
    scanf ("%"LZS"f",&y[i]);
  for (i=0; i<=m; i++)
    printf ("%"LZP"f\n", y[i]);

  vmblock = vminit();
  b = (mat4x4 **)vmalloc(vmblock, MMATRIX, n + 1, m + 1);
  if (! vmcomplete(vmblock))   /* allocations partially unsuccessful? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }
  printf ("\nInput of functional values (b[i][j][0][0], i=0..n, "
          "j=0..m):\n");
  for (i=0; i<=n; i++)
    for (j=0; j<=m; j++)
    {
      printf ("for index [%d][%d]: ", i, j);
      scanf  ("%"LZS"f ", &b[i][j][0][0]);
      printf ("%"LZP"f \n",  b[i][j][0][0]);
    }

  printf ("\nderivatives (b[i][0][1][0] und b[i][m][1][0]):\n");
  for (i=0; i<=n; i++)
   scanf ("%"LZS"f %"LZS"f", &b[i][0][0][1],&b[i][m][0][1]);
  for (i=0; i<=n; i++)
   printf ("%"LZP"f %"LZP"f\n",  b[i][0][0][1], b[i][m][0][1]);

  printf ("\nderivatives (b[0][j][0][1] und b[n][j][0][1]):\n");
  for (j=0; j<=m; j++)
    scanf ("%"LZS"f %"LZS"f", &b[0][j][1][0],&b[n][j][1][0]);
  for (j=0; j<=m; j++)
    printf ("%"LZP"f %"LZP"f\n", b[0][j][1][0], b[n][j][1][0]);

  printf ("\nderivatives (b[0][j][1][1] und b[n][j][1][1]):\n");
  for (j=0; j<=m; j++)
    scanf ("%"LZS"f %"LZS"f", &b[0][j][1][1], &b[n][j][1][1]);
  for (j=0; j<=m; j++)
    printf ("%"LZP"f %"LZP"f\n", b[0][j][1][1], b[n][j][1][1]);

  printf ("\nderivatives (b[i][0][1][1] und b[i][m][1][1]):\n");
  for (i=0; i<=n; i++)
    scanf ("%"LZS"f %"LZS"f", &b[i][0][1][1], &b[i][m][1][1]);
  for (i=0; i<=n; i++)
    printf ("%"LZP"f %"LZP"f\n", b[i][0][1][1], b[i][m][1][1]);

  printf ("\nInput b[i][j]:\n");
  for (i=0; i<=n; i++)
  {
    for (j=0; j<=m; j++)
      printf ("%10.5"LZP"f ", b [i][j][0][0]);
    printf ("\n");
  }

  error = bikub1 (n,m,b,x,y);

  if (error == 0)
  {
    printf ("\nOutput:\n");

    hilf = fibiku (n, m, b, x, y);

    printf ("\nvalue fibiku = %"LZP"f\n",hilf);

    printf ("\nenter lower boundary (x,y); upper boundary (x,y)\n");
    scanf  ("%"LZS"f %"LZS"f %"LZS"f %"LZS"f", &xlow, &ylow, &xup, &yup);
    printf ("%"LZP"f %"LZP"f %"LZP"f %"LZP"f\n", xlow, ylow, xup, yup);

    error = fibik2 (n, m, b, x, y, xlow, ylow, xup, yup, &hilf);

    if (error != 0)
      printf ("error in fibik2 = %d\n",error);
    else
      printf ("\nvalue fibik2 = %"LZP"f\n\n",hilf);
  }
  else
    printf ("error number = %d\n\n", error);

  return 0;
}
