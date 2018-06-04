#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /* for   umleiten, printf, scanf, NULL,      */
                         /*       fprintf, stderr, REAL, LZS, LZP,    */
                         /*       fehler_melden                       */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <kubsplin.h>    /*  for  spline                              */
#include <spliwert.h>    /*  for  spwert                              */
#include <splintab.h>    /*  for  sptab                               */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define TABELLENLAENGE  90



/***********************************************************************
*                                                                      *
* Test program: non periodic interpolating polynomial splines of       *
*               degree 3                                               *
*               Compute coefficients, tabulate and plot for Turbo-C    *
*               compilers                                              *
*                                                                      *
* Input is read from the first file in the command line.               *
* If none is given input comes from the keyboard.                      *
* Output other than an error code appears in the second file of the    *
* command line.                                                        *
* For calls with less than 2 files in the command line, all output     *
* appears on the screen.                                               *
* For Turbo C compiler the spline and its nodes are plotted.           *
* This can be voided by adding a third entry of "n" on the command line*
*                                                                      *
* Construction of input files:                                         *
* n                 (number of nodes)                                  *
* x1  y1            (1st node: x-value, y-value)                       *
* x2  y2            (2nd node: ditto)                                  *
* ...               ...                                                *
* xn  yn            (nth node: ditto)                                  *
* mass              (type of end point condition)                      *
* alfa beta         (left and right end point conditions)              *
*                                                                      *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL *x, *y,
       *b, *c, *d,
       *xtab, *ytab,
       alfa, beta,
       xmit,
       hilf,
       ausg[3],
       laenge;           /* Length of interpolating interval          */
  int  n,
       i,
       fehler,
       mass,
       lentab;
  void *vmblock;         /* List of dynamically allocated vectors     */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "number of nodes:                     ");
#endif
  scanf("%d", &n);

  vmblock = vminit();                             /* allocate storage */
  x    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  b    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  c    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  d    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  xtab = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  ytab = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "nodes (x[%2d],y[%2d]):               ", i, i);
#endif
    scanf("%"LZS"f%"LZS"f", x + i, y + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "type of end point conditions (0...4): ");
#endif
  scanf("%d", &mass);

  if (mass != 4)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "end point conditions alfa, beta:   ");
#endif
    scanf("%"LZS"f%"LZS"f", &alfa, &beta);
  }


  /* ------------ print input for control purposes ------------------ */

  printf("\n"
         "Compute and tabulate an interpolating cubic spline\n"
         "==================================================\n"
        );

  printf("\n\n"
         "number of nodes:              %d (= n+1)\n"
         "type of end point conditions: %d\n",
         n, mass
        );
  if (mass != 4)
    printf("end point conditions:         alfa = %20.13"LZP"f\n"
           "                              beta = %20.13"LZP"f\n",
           alfa, beta);

  printf("\n\n"
         "nodes:\n\n"
         " i          x[i]                y[i]\n"
         "--------------------------------------------\n"
        );
  for (i = 0; i < n; i++)
    printf("%2d %20.13"LZP"f %20.13"LZP"f\n", i, x[i], y[i]);


  /* -------- call  spline() to compute the coefficients  ----------- */

  fehler = spline(n, x, y, mass, alfa, beta, 0, b, c, d);

  if (fehler != 0)
  {
    fehler_melden("spline()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }


  /* ---------- print output ---------------------------------------- */

  printf("\n\n"
         "Spline coefficients:\n\n"
          " i        a[i]               b[i]               "
          "c[i]               d[i]\n"
          "------------------------------------------------"
          "-----------------------------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f"
           "  %17.13"LZP"f\n", i, y[i], b[i], c[i], d[i]);


  /* --- print function values and derivatives at 5 places --------- */

  printf("\n\n"
         "some functional and derivative values of the spline:\n\n"
         "      x0             y              "
         "y'             y''            y'''\n"
         "------------------------------------"
         "--------------------------------------\n");
  xmit = (x[n - 1] - x[0]) / (REAL)7.0;
  for (i = 1; i <= 5; i++)
  {
    hilf = spwert(n - 1, x[0] + i * xmit, y, b, c, d, x, ausg);
    printf("%14.10"LZP"f %14.10"LZP"f %14.10"LZP"f %14.10"LZP"f"
           " %14.10"LZP"f\n",
           x[0] + i * xmit, hilf, ausg[0], ausg[1], ausg[2]
          );
  }


  /* ------- as control, compute function values and values for ----- */
  /* ------- the derivatives of the spline at the nodes       ------- */

  printf("\n\n"
         "Function and derivative values at the nodes:\n\n"
         " i    x[i]          y[i]          "
         "y'[i]        y''[i]         y'''[i]\n"
         "----------------------------------"
         "--------------------------------------\n");
  for (i = 0; i < n; i++)
  {
    hilf = spwert(n - 1, x[i], y, b, c, d, x, ausg);
    printf("%2d %13.8"LZP"f %13.8"LZP"f %13.8"LZP"f %13.8"LZP"f"
           " %13.8"LZP"f\n",
           i, x[i], hilf, ausg[0], ausg[1], ausg[2]
          );
  }


  /* tabulate the spline from  x[0] to  x[n-1] with minimal step size */

  laenge = x[n - 1] - x[0];
  alfa = x[0]     - (REAL)0.015 * laenge;
  beta = x[n - 1] + (REAL)0.015 * laenge;

  fehler = sptab(n - 1, alfa, beta,
                 (beta - alfa) / (TABELLENLAENGE - 1 - n),
                 TABELLENLAENGE - 1, x, y, b, c, d, xtab, ytab,
                 &lentab);

  if (fehler != 0)
  {
    fehler_melden("sptab()", 20 + fehler, __FILE__, __LINE__);
    return 20 + fehler;
  }

  printf("\n\n"
         "Table of values for the spline from about x[0] to x[n] "
         "with minimal step size:\n\n"
         " i       xtab[i]              ytab[i]\n"
         "--------------------------------------------\n"
        );
  for (i = 0; i <= lentab; i++)
    printf("%2d %20.14"LZP"f %20.14"LZP"f\n", i, xtab[i], ytab[i]);


  /* ------------ plot spline if desired ---------------------------- */

  if (argc <= 3 || *argv[3] != 'n')
  {
    fehler = zeigkurv(lentab + 1, n, xtab, ytab, x, y);
    switch (fehler)
    {
      case 0:
        break;
      case 3:
        fehler_melden("zeigkurv(): lack of memory",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 4:
        fehler_melden("zeigkurv(): BGI graphics error",
                      30 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("zeigkurv(): other error",
                      30 + fehler, __FILE__, __LINE__);
    }
    if (fehler != 0)
      return 30 + fehler;
  }


  return 0;
}
