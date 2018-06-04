#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE thnpltst.c ------------------------ */

#include <basis.h>       /*  for  REAL, umleiten, scanf, LZS, ONE,    */
                         /*       ZERO, printf, LZP, stderr           */
#include <vmblock.h>     /*  for  vminit, vmalloc, VEKTOR, vmcomplete */
#include <thinplat.h>    /*  for  prob2, apprx2                       */
#include <zeigkrv2.h>    /*  for  zeigflaeche                         */



/* ------------------------------------------------------------------ */

int main
        (
         int  argc,
         char *argv[]
        )

/***********************************************************************
* Test program for the function prob2() which computes two-dimensional *
* surface splines, or thin plate splines.                              *
* The program takes the input from the first entry of the command line *
* and computes the spline function.                                    *
* Then the function  apprx2() makes a table of values for one parameter*
* line. The output goes into the second entry on the command line.     *
* If these are missing, standard in/output files are used.             *
*                                                                      *
* Form of an input file:                                               *
* =======================                                              *
* NX                 number of nodes                                   *
* M                  given derivative order                            *
* rho                smoothing parameter                               *
*                        = 0: Interpolation                            *
*                        > 0: smoothing spline                         *
* x[1] y[1] z[1]     nodes that are to be interpolated or smoothed     *
* x[2] y[2] z[2]                                                       *
*  ......                                                              *
* x[NX],y[NX],z[NX]                                                    *
*                                                                      *
* only if   rho > 0:                                                   *
* w[1]               smoothing value for first node                    *
* w[2]               ditto for second node                             *
* ...                                                                  *
* w[NX]              ... and for last node                             *
***********************************************************************/

{
  REAL *x,
       *y,
       *z,
       *w,
       *c,
       rho,
       x0,
       y0,
       y0min,
       y0max;
  REAL ***d,
       ***flaeche,
       xmin,
       xmax;
  int  xtablen,
       ytablen,
       j;
  int  fehler;
  int  NX;
  int  M;
  int  i;
  void *vmblock;   /* List of dynamically allocated vectors/matrices  */


  if ((fehler = umleiten(argc, argv))   /* asign input/output files   */
      != 0)                             /* to standard in/outputs     */
    return fehler;  /* 1 or 2 */


  scanf("%d%d%"LZS"f", &NX, &M, &rho);

  vmblock = vminit();                           /* initialize buffers */
  x = (REAL *)vmalloc(vmblock, VEKTOR, NX,  0);
  y = (REAL *)vmalloc(vmblock, VEKTOR, NX,  0);
  z = (REAL *)vmalloc(vmblock, VEKTOR, NX,  0);
  w = (REAL *)vmalloc(vmblock, VEKTOR, NX,  0);
  c = (REAL *)vmalloc(vmblock, VEKTOR, NX + (M * (M + 1)) / 2, 0);
  x--; y--; z--; w--; c--;                       /* Index shift !!!!! */
  if (! vmcomplete(vmblock))                      /* lack of memory ? */
  {
    fehler_melden("lack of available memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 1; i <= NX; i++)
    scanf("%"LZS"f%"LZS"f%"LZS"f", x + i, y + i, z + i),
    w[i] = ONE;
  if (rho > ZERO)
    for (i = 1; i <= NX; i++)
      scanf("%"LZS"f", w + i);


  /* --------------------- Print input data ------------------------- */

  printf("\n\n"
         "Surface spline\n"
         "==============\n\n\n"
         "given triples:\n\n"
         "  i      X(i)          Y(i)          Z(i)\n"
         "---------------------------------------------\n"
        );
  for (i = 1; i <= NX; i++)
    printf("%3d  %12.4"LZP"g  %12.4"LZP"g  %12.4"LZP"g\n",
           i, x[i], y[i], z[i]);
  printf("\n");
  if (rho > ZERO)
  {
    printf("smoothing with  rho = %12.4"LZP"g\n"
           "smoothing weights:\n\n"
           "  i      W(i)\n"
           "-----------------\n",
           rho
          );
    for (i = 1; i <= NX; i++)
      printf("%3d  %12.4"LZP"g\n", i, w[i]);
  }
  else
    printf("Computed interpolating spline function\n");
  printf("\ngiven derivative order = %d\n\n", M);


  fprintf(stderr, "Start prob2()...\n");
  fehler = prob2(NX, x, y, z, M, rho, w, c);
  fprintf(stderr, "prob2(): error = %d\n", fehler);

  if (fehler)
  {
    printf("*** error: error = %d\n", fehler);
    return 10 + fehler;
  }

  printf("\n"
         "coefficients of spline:\n"
         "=======================\n\n"
         "coefficients of kernel function part:\n\n"
         "  i       C(i)\n"
         "--------------------\n"
        );

  for (i = 1; i <= NX; i++)
    printf("%3d  %15.7"LZP"g\n", i, c[i]);
  printf("\n"
         "coefficients of polynomial part:\n\n"
         "  i       C(i)\n"
         "--------------------\n"
        );
  for (i = NX + 1; i <= NX + M * (M + 1) / 2; i++)
    printf("%3d  %15.7"LZP"g\n", i, c[i]);
  printf("\n");

  /* ----- tabulate along parameter line x = const = x(NX/2) -------- */
  x0    = x[NX / 2];
  y0min = y[1];
  y0max = y[1];
  for (i = 2; i <= NX; i++)
  {
    if (y[i] < y0min)
      y0min = y[i];
    if (y[i] > y0max)
      y0max = y[i];
  }

  printf("Table of values along x = const = %12.4"LZP"g:\n\n"
         "    X             Y              S(X,Y)\n"
         "----------------------------------------\n",
         x0
        );
  for (y0 = y0min; y0 <= y0max; y0 += (y0max - y0min) / (REAL)20.0)
    printf("%12.4"LZP"g  %12.4"LZP"g  %12.4"LZP"g\n",
           x0, y0, apprx2(x0, y0, NX, M, x, y, c));


  /* ---------------------- Tabulation and plot --------------------- */

  if (argc <= 3 || *argv[3] != 'n')     /* plot not suppressed?       */
  {

    xtablen = 40;
    ytablen = 20;

    d       = (REAL ***)vmalloc(vmblock, PMATRIX, 1,       NX);
    flaeche = (REAL ***)vmalloc(vmblock, PMATRIX, xtablen, ytablen);
    if (! vmcomplete(vmblock))
    {
      fehler_melden("lack of memory", 0, __FILE__, __LINE__);
      return 3;
    }

    xmin = x[1];
    xmax = x[1];
    for (i = 2; i <= NX; i++)
      xmin = (x[i] < xmin) ? x[i] : xmin,
      xmax = (x[i] > xmax) ? x[i] : xmax;

    xmin  -= (xmax  - xmin)  / TEN;
    xmax  += (xmax  - xmin)  / TEN;
    y0min -= (y0max - y0min) / TEN;
    y0max += (y0max - y0min) / TEN;

    for (i = 0; i < NX; i++)
      d[0][i][0] = x[i + 1],
      d[0][i][1] = y[i + 1],
      d[0][i][2] = z[i + 1];

    for (x0 = xmin, i = 0; i < xtablen;
         i++, x0 += (xmax - xmin) / (REAL)xtablen)
      for (y0 = y0min, j = 0; j < ytablen;
           j++, y0 += (y0max - y0min) / (REAL)ytablen)
        flaeche[i][j][0] = x0,
        flaeche[i][j][1] = y0,
        flaeche[i][j][2] = apprx2(x0, y0, NX, M, x, y, c);

    fehler = zeigflaeche(flaeche, xtablen, ytablen, d, 1, NX, 1, 1, 0);

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


  return 0;
}

/* -------------------------- END thnpltst.c ------------------------ */
