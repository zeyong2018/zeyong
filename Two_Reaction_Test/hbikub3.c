#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>      /*  for  REAL, umleiten, printf, scanf,       */
                        /*       fehler_melden, LZS, LZP              */
#include <vmblock.h>    /*  for  mat4x4, vminit, vmalloc, VEKTOR,     */
                        /*       MMATRIX, PMATRIX, vmcomplete, vmfree */
#include <bikub.h>      /*  for  bikub3                               */
#include <zeigkrv2.h>   /*  for  zeigflaeche                          */

/*
   Test program: Bicubic splines 3
*/

int main (int argc, char *argv[])
{
  int  m, n, i, j, k, error;
  REAL *x, *y, ***fn;
  mat4x4 **b;
  int    fehler;
  REAL   ***d,
         ***flaeche,
         x0,
         y0,
         xmin,
         xmax,
         ymin,
         ymax,
         wert,
         hx,
         hy;
  int    xtablen,
         ytablen;
  void   *vmblock;

  if ((fehler = umleiten(argc, argv))    /* assign input/output files */
      != 0)                              /* to standard ones          */
    return fehler;  /* 1 or 2 */

  printf ("Bicubic splines via normal vectors:\n\n");

  printf ("Input of number of x- and y-intervals\n");
  scanf  ("%d %d", &m,&n);

  vmblock = vminit();
  x  = (REAL *)   vmalloc(vmblock, VEKTOR,  m + 1, 0);
  y  = (REAL *)   vmalloc(vmblock, VEKTOR,  n + 1, 0);
  b  = (mat4x4 **)vmalloc(vmblock, MMATRIX, m + 1, n + 1);
  fn = (REAL ***) vmalloc(vmblock, PMATRIX, m + 1, n + 1);
  if (! vmcomplete(vmblock))   /* allocations partially unsuccessful? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  printf ("Input of x-interval end points\n");
  for (i=0; i<=m; i++)
    scanf ("%"LZS"f",&x[i]);

  printf ("Input of end points for y-intervals\n");
  for (i=0; i<=n; i++)
    scanf ("%"LZS"f ",&y[i]);

  printf ("Function values:\n");
  for (i=0; i<=m; i++)
    for (j=0; j<=n; j++)
      scanf ("%"LZS"f ", &b[i][j][0][0]);

  printf ("Normal vectors \n");
  for (i=0; i<=m; i++)
    for (j=0; j<=n; j++)
      scanf ("%"LZS"f %"LZS"f %"LZS"f",
             &fn[i][j][0],&fn[i][j][1],&fn[i][j][2]);


  printf ("Input of b[i][j][0][0]:\n");
  for (i=0; i<=m; i++)
  {
    for (j=0; j<=n; j++)
      printf ("%10.5"LZP"f ", b [i][j][0][0]);
    printf ("\n");
  }

  error = bikub3 (m,n,b,x,y,fn);

  if (error == 0)
  {
    printf ("\nOutput:\n");
    printf ("mat = :\n");
    for (i=0; i<=m-1; i++)
    {
      printf ("\n\n");
      for (j=0; j<=n-1; j++)
      {
        printf ("\n");
        for (k=0; k<=3; k++)
          printf ("%10.7"LZP"f  %10.7"LZP"f  %10.7"LZP"f  "
                  "%10.7"LZP"f\n",
                  b[i][j][k][0], b[i][j][k][1],
                  b[i][j][k][2], b[i][j][k][3]);
      }
    }
  }
  else
    printf ("error code = %d\n", error);


  /* ---------------------- Tabulation and plot --------------------- */

  if (argc <= 3 || *argv[3] != 'n')     /* plot not suppressed?       */
  {

    xtablen = 40;
    ytablen = 25;

    d       = (REAL ***)vmalloc(vmblock, PMATRIX, m + 1,   n + 1);
    flaeche = (REAL ***)vmalloc(vmblock, PMATRIX, xtablen, ytablen);
    if (! vmcomplete(vmblock))
    {
      fehler_melden("lack of memory", 0, __FILE__, __LINE__);
      return 3;
    }

    for (i = 0; i <= m; i++)
      for (j = 0; j <= n; j++)
        d[i][j][0] = x[i],
        d[i][j][1] = y[j],
        d[i][j][2] = b[i][j][0][0];

    xmin = x[0];
    xmax = x[m];
    ymin = y[0];
    ymax = y[n];

    hx = (xmax - xmin) / (REAL)xtablen;
    hy = (ymax - ymin) / (REAL)ytablen;
    for (x0 = xmin, i = 0; i < xtablen; i++, x0 += hx)
      for (y0 = ymin, j = 0; j < ytablen; j++, y0 += hy)
      {
        flaeche[i][j][0] = x0,
        flaeche[i][j][1] = y0,
        bsval(m, n, b, x, y, x0, y0, &wert);
        flaeche[i][j][2] = wert;
      }

    fehler = zeigflaeche(flaeche, xtablen, ytablen, d, m + 1, n + 1,
                         1, 1, 1);

    switch (fehler)
    {
      case 0:
        break;
      case 1:
        fehler_melden("zeigflaeche(): nv too small",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 2:
        fehler_melden("zeigflaeche(): nw too small",
                     30 + fehler, __FILE__, __LINE__);
        break;
      case 3:
        fehler_melden("zeigflaeche(): m too small",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 4:
        fehler_melden("zeigflaeche(): n too small",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 5:
        fehler_melden("zeigflaeche(): nw or n too large",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 6:
        fehler_melden("zeigflaeche(): graphics error",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 7:
        fehler_melden("zeigflaeche(): lack of memory",
                      30 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("zeigflaeche(): other error",
                      30 + fehler, __FILE__, __LINE__);
    }
  }


  vmfree(vmblock);

  return 0;
}
