#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, fprintf,   */
                         /*       stderr, NULL, REAL, LZS, ZERO, LZP, */
                         /*       fehler_melden                       */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <subsplin.h>    /*  for  akima                               */
#include <spliwert.h>    /*  for  spwert                              */
#include <splintab.h>    /*  for  sptab                               */
#include <zeigkurv.h>    /*  for  zeigkurv                            */


#define TABELLENLAENGE  160

/***********************************************************************
* Test program for the function akima() from the module subsplin to    *
* compute interpolating Akima subsplines.                              *
*                                                                      *
* Scope of the program :                                               *
* ======================                                               *
* The program reads the input data from the first entry file in the    *
* command line and computes the corresponding subspline.               *
* Then using the function spwert() from the module splintab we compute *
* several function and derivative values and make a table of values for*
* the spline using sptab() from splintab.                              *
* The output is recorded onto the second file in the command line.     *
* If there is only one entry, output is dumped onto the screen.        *
* If there is not even one entry on the command line, input is from    *
* the keyboard.                                                        *
* Several test data sets are provided in akitst.ei*. The output should *
* be identical to the contents of akitst.au* .                         *
* For Turbo C compilers we plot the subspline and the nodes, unless    *
* the third command entry is "n".                                      *
*                                                                      *
* Form of input file :                                                 *
* ====================                                                 *
* n               number of nodes                                      *
* x[0],y[0]       1st node                                             *
* x[1],y[1]       2nd node                                             *
*   ...           ...                                                  *
* x[n-1],y[n-1]   last node                                            *
* rund            Parameter for rounding :                             *
*                 0 < rund < 1: round corners,                         *
*                 for other values of rund : no rounding               *
* perio           1: periodic Akima interpolation                      *
*                 0: non periodic Akima interpolation                  *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL *x,           /* the vectors x, y, b, c and d must have room   */
       *y,           /* for maximally  n+[n/2] entries comprised of   */
       *b,           /* n original nodes and [n/2] potential new nodes*/
       *c,           /* from rounding.                                */
       *d,
       *xtab,
       *ytab,
       ausg[3],
       rund,
       x0,
       xanf,
       xend,
       xstep,
       hilf;
  int  n,
       n0,
       i,
       perio,
       fehler,
       lentab;
  void *vmblock;                 /* List of dynamical allocations     */


  if ((fehler = umleiten(argc, argv))   /* assign in/output files to  */
      != 0)                             /* standard ones              */
    return fehler;  /* 1 or 2 */


  /* ---------------------- read input ------------------------------ */

#ifdef INTERAKTIV
  fprintf(stderr, "\nNumber of nodes :                     ");
#endif
  scanf("%d", &n);

  n0 = n + n / 2;
  vmblock = vminit();                         /* initialize storage  */
  x    = (REAL *)vmalloc(vmblock, VEKTOR, n0,             0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, n0,             0);
  b    = (REAL *)vmalloc(vmblock, VEKTOR, n0,             0);
  c    = (REAL *)vmalloc(vmblock, VEKTOR, n0,             0);
  d    = (REAL *)vmalloc(vmblock, VEKTOR, n0,             0);
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
    fprintf(stderr, "node (x[%2d],y[%2d]):                    ",
                    i, i);
#endif
    scanf("%"LZS"f%"LZS"f", x + i, y + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "rounding corners (0 < rund < 1):                  ");
#endif
  scanf("%"LZS"f", &rund);

#ifdef INTERAKTIV
  fprintf(stderr, "periodic interpolation (1) or not (0): ");
#endif
  scanf("%d", &perio);


  /* ----------------- print out input data ------------------------- */

  printf("interpolating Akima subspline\n"
         "=============================\n\n\n"
         "nodes:\n\n"
         "  i     x[i]        y[i]\n"
         "---------------------------\n"
        );
  for (i = 0; i < n; i++)
    printf("%3d%12.5"LZP"f%12.5"LZP"f\n", i, x[i], y[i]);
  printf("\n");
  if (rund == ZERO)
    printf("Corners are not rounded.\n");
  else
    printf("rounding parameter: %7.3"LZP"f\n", rund);
  if (perio)
    printf("periodic Akima interpolation\n\n");


  /* ------------------ compute subspline -------------------------- */

  n--;
  fehler = akima(&n, n0 - 1, x, y, perio, rund, b, c, d);
  n++;

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("akima(): less than 5 nodes",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("akima(): nodes not monotone "
                    "steigend", 10 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("akima(): y[0] != y[n-1]  for a periodic spline",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 4:
      fehler_melden("akima(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 5:
      fehler_melden("akima(): lack of space for all roundings",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("akima(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ---------------- print subspline coefficients ------------------ */

  if (rund > ZERO)
  {
    printf("\n\nnodes with extra nodes at corners:\n\n"
            "  i     x[i]        y[i]\n"
            "---------------------------\n"
          );
    for (i = 0; i < n; i++)
      printf("%3d%12.5"LZP"f%12.5"LZP"f\n", i, x[i], y[i]);
  }

  printf("\n\nSubspline coefficients:\n\n"
         "  i       a[i]           b[i]           c[i]"
         "           d[i]\n"
         "-----------------------------------------------"
         "----------------\n"
        );

  for (i = 0; i < n - 1; i++)
    printf("%3d%15.6"LZP"f%15.6"LZP"f%15.6"LZP"f%15.6"LZP"f\n",
           i, y[i], b[i], c[i], d[i]);


  /* ------- compute several function and derivative values --------- */

  printf("\n\nExample of using spwert():\n"
         "Function and derivative values at several points\n\n"
         "   x            S(x)         S'(x)       "
         " S''(x)       S'''(x)\n"
         "----------------------------------------"
         "----------------------\n"
        );
  for (x0 = x[0]; x0 <= x[n - 1]; x0 += (x[n - 1] - x[0]) / (REAL)15.0)
  {
    hilf = spwert(n - 1, x0, y, b, c, d, x, ausg);
    printf("%10.5"LZP"f%13.5"LZP"f%13.5"LZP"f%13.5"LZP"f"
           "%13.5"LZP"f\n", x0, hilf, ausg[0], ausg[1], ausg[2]);
  }


  /* -------------- make table of values with  sptab() -------------- */

  xstep  = x[n - 1] - x[0];
#if 0
  xanf   = x[0]     - (REAL)0.015 * xstep;
  xend   = x[n - 1] + (REAL)0.015 * xstep;
#else
  xanf   = x[0]     - (REAL)0.000 * xstep;
  xend   = x[n - 1] + (REAL)0.000 * xstep;
#endif
  xstep /= TABELLENLAENGE - 4 - n;

  fehler = sptab(n - 1, xanf, xend, xstep, TABELLENLAENGE - 1,
                 x, y, b, c, d, xtab, ytab, &lentab);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("sptab(): xanf > xend",
                    20 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("sptab(): step size <= 0",
                    20 + fehler, __FILE__, __LINE__);
    default:
      fehler_melden("sptab(): unknown error",
                    20 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 20 + fehler;

  printf("\n\nExample for using sptab():\n"
         "Table of values from x =%10.5"LZP"f to %11.5"LZP"f,"
         " Step size = %10.5"LZP"f\n\n", xanf, xend, xstep);

  printf("  i        xtab[i]        ytab[i]\n"
         "---------------------------------\n");
  for (i = 0; i <= lentab; i++)
    printf("%3d%15.6"LZP"f%15.6"LZP"f\n", i, xtab[i], ytab[i]);


  /* ------------ plot spline if desired                  ----------- */

  if (argc <= 3 || *argv[3] != 'n')
  {
    fehler = zeigkurv(lentab + 1, n, xtab, ytab, x, y);
    switch (fehler)
    {
      case 0:
        break;
      case 3:
        fehler_melden("zeigkurv(): lack of memory ",
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
