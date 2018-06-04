#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>      /*  for  REAL, umleiten, printf, scanf,       */
                        /*       fehler_melden, LZS, LZP              */
#include <vmblock.h>    /*  for  mat4x4, vminit, vmalloc, VEKTOR,     */
                        /*       MMATRIX, PMATRIX, vmcomplete, vmfree */
#include <bikub.h>      /*  for  bikub1                               */
#include <zeigkrv2.h>   /*  for  zeigflaeche                          */

/*
   Test program: Bicubic splines 1
*/

int main (int argc, char *argv[])
{
  int  m, n, i, j, k, error;
  REAL *x, *y;
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

#if 0
    FILE *inp = fopen ("hbikub1.inp", "r");    /* test */
    #define scanf(f,v) fscanf (inp,f,v)        /* test */
#endif

  printf ("Bicubic splines via derivatives:\n\n");

  printf ("number of x- and y-intervals\n");
  scanf  ("%d", &m);
  scanf  ("%d", &n);

  vmblock = vminit();
  x = (REAL *)   vmalloc(vmblock, VEKTOR,  m + 1, 0);
  y = (REAL *)   vmalloc(vmblock, VEKTOR,  n + 1, 0);
  b = (mat4x4 **)vmalloc(vmblock, MMATRIX, m + 1, n + 1);
  if (! vmcomplete(vmblock))   /* allocations partially unsuccessful? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  printf ("%d x interval end points:\n", m+1);
  for (i=0; i<=m; i++)
    scanf ("%"LZS"f",&x[i]);

  printf ("%d y interval end points :\n", n+1);
  for (i=0; i<=n; i++)
    scanf ("%"LZS"f",&y[i]);

  printf ("Function values  b[i][j][0][0] (i=0(1)%d, j=0(1)%d):\n",
          m, n);
  for (i=0; i<=m; i++)
    for (j=0; j<=n; j++)
      scanf ("%"LZS"f", &b[i][j][0][0]);

  printf ("derivatives (b[0][j][1][0] and b[%d][j][1][0], "
          "j=0(1)%d)\n", m, n);
  for (j=0; j<=n; j++)
  {
    scanf ("%"LZS"f", &b[0][j][1][0]);
    scanf ("%"LZS"f", &b[m][j][1][0]);
  }
  printf ("derivatives (b[i][0][0][1] and b[i][%d][0][1], "
          "i=0(1)%d)\n", n, m);
  for (i=0; i<=m; i++)
  {
   scanf ("%"LZS"f", &b[i][0][0][1]);
   scanf ("%"LZS"f", &b[i][n][0][1]);
  }
  printf ("derivatives (b[0][j][1][1] and b[%d][j][1][1], "
          "j=0(1)%d)\n", m, n);
  for (j=0; j<=n; j++)
  {
    scanf ("%"LZS"f", &b[0][j][1][1]);
    scanf ("%"LZS"f", &b[m][j][1][1]);
  }
  printf ("derivatives (b[i][0][1][1] and b[i][%d][1][1], "
          "i=0(1)%d)\n", n, m);
  for (i=0; i<=m; i++)
  {
    scanf ("%"LZS"f", &b[i][0][1][1]);
    scanf ("%"LZS"f", &b[i][n][1][1]);
  }
#if 0
  {
    REAL hilf;
    x[0] = y[0] = b[0][0][0][0];
    hilf = (b[m][n][0][0] - x[0]) / m;
    for (i=1; i<=m; ++i)
      x[i] = x[i-1] + hilf;
    hilf = (b[m][n][0][0] - y[0]) / n;
    for (j=1; j<=n; ++j)
      y[j] = y[j-1] + hilf;
  }
#endif
  printf ("Input  b[i][j][0][0] (i=0(1)%d, j=0(1)%d):\n", m,n);
  for (i=0; i<=m; i++)
  {
    for (j=0; j<=n; j++)
      printf ("%10.5"LZP"g ", b [i][j][0][0]);
    printf ("\n");
  }

  error = bikub1 (m,n,b,x,y);

  if (!error)
  {
    printf ("\nOutput :\n");
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
    printf ("error code for (bikub1) = %d\n", error);


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
