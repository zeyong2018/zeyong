#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>      /*  for  REAL, freopen, stdin, NULL, fprintf, */
                        /*       stderr, printf, scanf, LZS, LZP,     */
                        /*       fehler_melden, umleiten              */
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vminit, MATRIX  */
#include <bikub.h>      /*  for  mokube, valbez                       */
#include <zeigkrv2.h>   /*  for  zeigkrv2                             */



/* ------------------------------------------------------------------ */

int main
        (
         int  argc,
         char *argv[]
        )

/*
   Test program: one-dimensional Bezier splines, modified method
                 (the subroutine kubbez is used in mokube and
                  thus tested indirectly here.)
*/

{
  int    n,              /* number of spline pieces                   */
         tablen,         /* size of table of curve points             */
         i, j;
  REAL   **b,            /* [0..3n] vector with weight points         */
         **d,            /* [0..n] vector with Bezier points          */
         **kurvtab,      /* [0..tablen] vector with curve points      */
         eps;            /* accuracy of interpolation                 */
  int    dim;            /* space dimension (2 or 3)                  */
  int    fehler;         /* error code of umleiten(), mokube() resp.  */
                         /* zeigkrv2()                                */
  REAL   t,
         x, y, z;
  void   *vmblock;       /* liste of dynamically allocated vectors    */
                         /* and matrices                              */


  if ((fehler = umleiten(argc, argv))    /* assign input/output files */
      != 0)                              /* to standard ones          */
    return fehler;  /* 1 or 2 */

  printf("Test program: one-dimensional Bezier splines, "
         "modified method\n");

  dim    = 3;                             /* compute curve in space   */
  tablen = 100;                           /* 100 points on the curve  */

  printf("Input n\n");
  scanf("%d", &n);

  vmblock = vminit();                      /* initialize memory block */
  b       = (REAL **)vmalloc(vmblock, MATRIX, 3 * n + 1, dim);
  d       = (REAL **)vmalloc(vmblock, MATRIX, n + 1,     dim);
  kurvtab = (REAL **)vmalloc(vmblock, MATRIX, tablen,    dim);
  if (! vmcomplete(vmblock))   /* allocations partially unsuccessful? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  printf("Input d [i][k] (i=0,...,n, k=0,1,2) \n");
  for (i = 0; i <= n; i++)
    for (j = 0; j < dim; j++)
      scanf("%"LZS"f", &d[i][j]);
  printf("Input eps\n");
  scanf("%"LZS"f", &eps);

  printf("\nn  = %d\n", n);

  printf("Input:\n");
  for (i = 0; i <= n; i++)
  {
    for (j = 0; j < dim; j++)
      printf("%21.13"LZP"f ", d[i][j]);
    printf("\n");
  }
  printf("eps = %"LZP"f\n",eps);

  fehler = mokube(b, d, n, dim, eps);               /* modified cubic */
                                                    /* Bezier spline  */

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("mokube(): m < 2  or  dim < 2  or  dim > 3",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("mokube(): eps too small",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("mokube(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;

  printf("\nOutput:\n");
  for (i = 0; i <= 3 * n; i++)
  {
    for (j = 0; j < dim; j++)
      printf("%21.13"LZP"f ", b[i][j]);
    printf("\n");
  }

  printf("\nEvaluation at the parameter values of "
         "given nodes:\n"
         "   t               x                     "
         "y                     z\n"
         "-----------------------------------------"
         "-------------------------------\n");
  for (i = 0; i <= n; i++)
  {
    t = (REAL)i / (REAL)n;
    valbez(t, n, dim, b, &x, &y, &z);
    printf("%6.3"LZP"f%22.13"LZP"f%22.13"LZP"f", t, x, y);
    if (dim > 2)
      printf("%22.13"LZP"f", z);
    printf("\n");
  }


  for (i = 0; i < tablen; i++)                 /* compute `tablen'    */
  {                                            /* points on the curve */
    t = (REAL)i / (REAL)(tablen - 1);
    valbez(t, n, dim, b, kurvtab[i], kurvtab[i] + 1, kurvtab[i] + 2);
  }

  /* ----------- plot spline curve if desired              ---------- */
  /* ----------- (attention: the z component of spatial    ---------- */
  /* ----------- curves is ignored.)                       ---------- */

  if (argc <= 3 || *argv[3] != 'n')     /* plot not suppressed?       */
  {
    fehler = zeigkrv2(tablen, kurvtab, n + 1, d);
    switch (fehler)
    {
      case 0:
        break;
      case 3:
        fehler_melden("zeigkrv2(): lack of memory",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 4:
        fehler_melden("zeigkrv2(): graphics error",
                      30 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("zeigkrv2(): other error",
                      30 + fehler, __FILE__, __LINE__);
    }
    if (fehler != 0)
      return 30 + fehler;
  }


  return 0;
}



/*
                   T E S T  E X A M P L E


Input:

5
1.   2.   3.
1.12 4.67 5.12
0.4  5.2  8.9
4.5  3.2  3.2
2.6  7.3  0.2
4.7  5.2 -3.8
1.e-5

Output:

n  = 5
Input:
      1.0000000000000       2.0000000000000       3.0000000000000
      1.1200000000000       4.6700000000000       5.1200000000000
      0.4000000000000       5.2000000000000       8.9000000000000
      4.5000000000000       3.2000000000000       3.2000000000000
      2.6000000000000       7.3000000000000       0.2000000000000
      4.7000000000000       5.2000000000000      -3.8000000000000
eps = 0.000010

Output:
      1.0000000000000       2.0000000000000       3.0000000000000
      1.2749900460119       2.9718045249147       3.3127922202957
      1.5499800920238       3.9436090498293       3.6255844405914
      1.1199966887912       4.6700038340300       5.1200019265261
      0.6900132855587       5.3963986182306       6.6144194124608
     -0.4449435669184       5.8773836617172       9.2904621640344
      0.4000013394208       5.1999984490893       8.8999992207017
      1.2449462457599       4.5226132364614       8.5095362773689
      4.0697929109152       2.6868577677190       5.0525676391297
      4.4999996651448       3.2000003877277       3.2000001948246
      4.9302064193744       3.7131430077364       1.3474327505195
      2.9657732626784       6.5751837164961       1.0992665001484
      2.6000000000000       7.3000000000000       0.2000000000000
      2.2342267373216       8.0248162835039      -0.6992665001484
      3.4671133686608       6.6124081417519      -2.2496332500742
      4.7000000000000       5.2000000000000      -3.8000000000000

Evaluation at the parameter values of given nodes:
   t               x                     y                     z
------------------------------------------------------------------------
 0.000       1.0000000000000       2.0000000000000       3.0000000000000
 0.200       1.1199966887912       4.6700038340300       5.1200019265261
 0.400       0.4000013394208       5.1999984490893       8.8999992207017
 0.600       4.4999996651448       3.2000003877277       3.2000001948246
 0.800       2.6000000000000       7.3000000000000       0.2000000000000
 1.000       4.7000000000000       5.2000000000000      -3.8000000000000
*/
