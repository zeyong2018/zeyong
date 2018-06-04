#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, NULL,      */
                         /*       fprintf, stderr, REAL, LZS, LZP,    */
                         /*       ZERO, fehler_melden                 */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <approx.h>      /*  for  gfq, horner                         */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define MAXTLEN  21                        /* max length of table     */



/***********************************************************************
* A test program for gfq()  to compute the coefficients of a root-mean-*
* square approximation using the discrete Gaussian least square method.*
*                                                                      *
* Scope of the test program :                                          *
* ===========================                                          *
* The program reads from the first mentioned file and computes the     *
* root-mean-square approximating polynomial.                           *
* Subsequently  horner() creates a table of values for this polynomial.*
* The output is stored in the second file name in the command line.    *
* If there is only one file, the output is dumped to the screen.       *
* If the first file is missing also, the input must be typed in.       *
* Several test data sets sre contained in the files gfqtst.ei*.        *
* The output must have the values in  gfqtst.au* .                     *
* The examples  gfqtst.ei0 and gfqtst.ei1 come from :                  *
*           Numerische Mathematik fuer Ingenieure;                     *
*           G. Engeln-Muellges, F. Reutter, 6th edition 1988;          *
*           Example 8.1.3.2 (P8), p. 548 - 553. [ENGE88]               *
* The examples  gfqtst.ei2, ..., gfqtst.ei9 come from                  *
*           J. Stoer, R Bulirsch:                                      *
*           Introduction to Numerical Analysis, 3rd edition, Springer, *
*           1991., Example 2, ch. 4.8.3, p. 207 ff; [STOE91]           *
* In  each of these the theoretical solution is given by the           *
* coefficients 0, 1, 1, 1.                                             *
* If a Turbo C-Compiler is used, the approximating function and its    *
* nodes are displayed graphically. This can be avoided by setting the  *
* third command line parameter equal to "n".                           *
*                                                                      *
* Construction of input files:                                         *
* ============================                                         *
* n                         max degree of approximating polynomial     *
* m                         number of nodes                            *
* x[0],y[0],w[0]            node, function value, weight               *
* x[1],y[1],w[1]            ditto                                      *
* ...                       ...                                        *
* x[m-1],y[m-1],w[m-1]      ditto                                      *
***********************************************************************/

int main(int argc, char *argv[])

{
  int  i,              /* Loop variable                               */
       fehler,         /* error code for gfq()                        */
       n,              /* max degrre of approximating polynomial      */
       m;              /* number of nodes                             */
  REAL *x,             /* vector of nodes                             */
       *y,             /* vector of function values at nodes          */
       *w,             /* vector of weights                           */
       *c,             /* coefficients for the polynomial             */
       tabanf,         /* left end point for tabulating               */
       tabend,         /* right end point for tabulating              */
       xx,             /* Loop variable for tabulating                */
       schritt,        /* step size for tabulating                    */
       *xtab,          /* table of values; x-values and               */
       *ytab;          /*                  y-values                   */
  void *vmblock;       /* List of dynamically allocated vectors       */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nPolynomial degree n:          ");
#endif
  scanf("%d", &n);                /* read polynomial degree          */
#ifdef INTERAKTIV
  fprintf(stderr, "number m of nodes:            ");
#endif
  scanf("%d", &m);                            /* read number of nodes */

  vmblock = vminit();                             /* allocate storage */
  x    = (REAL *)vmalloc(vmblock, VEKTOR, m, 0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, m, 0);
  w    = (REAL *)vmalloc(vmblock, VEKTOR, m, 0);
  c    = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  xtab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN, 0);
  ytab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN, 0);
  if (! vmcomplete(vmblock))                     /* lack of storage ? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i < m; i++)               /* reads nodes, function      */
  {                                     /* values and weights         */
#ifdef INTERAKTIV
    fprintf(stderr, "node, function value, weight: ");
#endif
    scanf("%"LZS"f%"LZS"f%"LZS"f", x + i, y + i, w + i);
  }


  /* ------------------- put out input ----------------------------- */

  printf("\n"
         "Compute least square polynomial using the "
         "normal equations\n"
         "=========================================="
         "================\n\n\n"
         "max degree of polynomial: %5d\n"
         "number of nodes:          %5d\n\n\n"
         "Table of values: nodes, function values, weights:\n\n"
         " i            x[i]                   f[i]                   "
         "w[i]\n"
         "------------------------------------------------------------"
         "-----------\n",
         n, m
        );

  for (i = 0; i < m; i++)
    printf("%2d   %20.12"LZP"f   %20.12"LZP"f   %20.12"LZP"f\n",
           i, x[i], y[i], w[i]);


  /* ------ compute polynomial coeficients -------------------------- */

  fehler = gfq(n, m - 1, x, y, w, c);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("gfq(): n < 1 or m <= n",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("gfq(): Normal equations could not be solved",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("gfq(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("gfq(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /*  print coefficients of the least square approximating polynomial */

  printf("\n\n"
         "coefficients of the approximating polynomial:\n\n"
         " i           c[i]\n"
         "-------------------------\n");
  for (i = 0; i <= n; i++)
    printf("%2d   %20.12"LZP"f\n", i, c[i]);


  /* ------ form a table of values for the polynomial --------------- */

  schritt = x[m - 1] - x[0];
#if 0
  tabanf = x[0]     - (REAL)0.05 * schritt;       /* general case     */
  tabend = x[m - 1] + (REAL)0.05 * schritt;
#else                       /* these values are chosen specifically   */
  tabanf = ZERO;            /* for  gfqtst.ei2 ... gfqtst.ei9,        */
  tabend = (REAL)0.11;      /* in order to see the increasing de-     */
#endif                      /* viations from the theoretical solution.*/
  schritt = (tabend - tabanf) / (MAXTLEN - 1);

  printf("\n\n"
         "table of values:\n\n"
         "          x                   f(x)\n"
         "-----------------------------------------\n");

  for (xx = tabanf, i = 0; i < MAXTLEN; i++, xx += schritt)
  {
    xtab[i] = xx;
    ytab[i] = horner(n, c, xx);
    printf("%20.12"LZP"f %20.12"LZP"f\n", xtab[i], ytab[i]);
  }


  /* ---------------- if desired, plot polynomial ------------------- */

  if (argc <= 3 || *argv[3] != 'n')
  {
    fehler = zeigkurv(MAXTLEN, m, xtab, ytab, x, y);
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
