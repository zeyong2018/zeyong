#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, NULL, ONE, */
                         /*       fprintf, stderr, REAL, LZS, LZP,    */
                         /*       TWO, POW, FIVE, TEN, fehler_melden, */
                         /*       horner                              */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <approx.h>      /*  for  pol_appr, opolkoeff, opolwert       */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define MAXTLEN  501



int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
* Test program for  pol_appr()  to compute the coefficients of an nth  *
* degree approximate polynomial based on discrete orthogonal           *
* polynomials                                                          *
*                                                                      *
* Scope of this program:                                               *
* ======================                                               *
* The input is read from the first named file in the command line and  *
* dimensions the neede vector dynamically from the input data.         *
* The nodes are equidistant in [0,1]; the y-values are computed for the*
* test function                                                        *
*                f(x) = 1 / (1 + (10x - 5) * (10x - 5))                *
* All weights are set to 1.                                            *
* After computing the least square polynomial, it is evaluated at      *
* MAXTLEN equistant points in the interval [0;1].                      *
* The output is directed to the second named fiel in the command line. *
* If this does not exist, all output is dumped onto the screen. If     *
* no file names occur there, input must be executed from the keyboard. *
* Several test data sets are available in the file  opoltst.ei* .      *
* The output should be nearly identical to the data in  opoltst.au* .  *
* For Turbo C compilers the resulting function is plotted. This can be *
* avoided by entering "n" as the third parameter on the command line.  *
*                                                                      *
* Construction of input file:                                          *
* ===========================                                          *
* stuetz    number of nodes                                            *
* grad      degree of least squares polynomial                         *
***********************************************************************/

{
  int  stuetz;          /* number of nodes                            */
  int  grad;            /* degree of least squares polynomial         */
  int  fehler;
  int  i;
  REAL *c;
  REAL *b;
  REAL *d;
  REAL *x;
  REAL *y;
  REAL *w;
  REAL *a;              /* [0..grad]-vector with coefficients of      */
                        /* least squares polynomial in representation */
                        /* a[0] + a[1]*x^1 + ... + a[grad]*x^grad     */
  REAL *xtab;           /* Table of values for solution: x-values     */
  REAL *ytab;           /*                               y-values     */
  void *vmblock;        /* List of dynamically allocated vectors      */


  if ((fehler = umleiten(argc, argv))   /* asign the input/output     */
      != 0)                             /* files to the standard files*/
    return fehler;  /* 1 or 2 */


  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "number of nodes:      ");
#endif
  scanf("%d", &stuetz);

#ifdef INTERAKTIV
  fprintf(stderr, "degree of polynomial: ");
#endif
  scanf("%d", &grad);

  vmblock = vminit();                 /* initialize storage buffers   */
  c    = (REAL *)vmalloc(vmblock, VEKTOR, grad + 1, 0);
  b    = (REAL *)vmalloc(vmblock, VEKTOR, grad + 1, 0);
  d    = (REAL *)vmalloc(vmblock, VEKTOR, grad + 1, 0);
  x    = (REAL *)vmalloc(vmblock, VEKTOR, stuetz,   0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, stuetz,   0);
  w    = (REAL *)vmalloc(vmblock, VEKTOR, stuetz,   0);
  a    = (REAL *)vmalloc(vmblock, VEKTOR, grad + 1, 0);
  xtab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN,  0);
  ytab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN,  0);
  if (! vmcomplete(vmblock))                       /* lack of memory? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }


  /* -------------- compute nodes and weights ----------------------- */

  for (i = 0; i < stuetz; i++)   /* use function 1/(1+(10x-5)(10x-5)) */
  {                              /* at equidistant nodes in [0,1]     */
    if (stuetz != 1)
      x[i] = (REAL)i / (REAL)(stuetz - 1);
    else
      x[i] = (REAL)i;
    y[i] = ONE / (ONE + sqr(TEN * x[i] - FIVE));
  }

  for (i = 0; i < (stuetz + 1) / 2; i++)
#if 1
    w[(stuetz - 1) / 2 - i] = ONE,               /* set all weights   */
    w[stuetz / 2       + i] = ONE;               /* equal to one      */
#endif
#if 0
    w[(stuetz - 1) / 2 - i] = sqr(i + 1),    /* increase weights near */
    w[stuetz / 2       + i] = sqr(i + 1);    /* ends                  */
#endif
#if 0
    w[(stuetz - 1) / 2 - i] = POW(TWO, (REAL)(i + 1)),  /* increase   */
    w[stuetz / 2       + i] = POW(TWO, (REAL)(i + 1));  /* weights    */
                                                        /* even more  */
                                                        /* at ends    */

#endif


  /* -------------------- print test data --------------------------- */

  printf("\n"
         "finding a least square approximating polynomial"
         " using orthogonal polynomials\n"
         "==============================================="
         "=============================\n\n\n"
         "maximal degree of polynomial: %5d\n"
         "number of nodes:              %5d\n\n\n"
         "Table of x-values, y-values and weights:\n\n"
         " i            x[i]                   y[i]                   "
         "w[i]\n"
         "------------------------------------------------------------"
         "-----------\n",
         grad, stuetz
        );
  for (i = 0; i < stuetz; i++)
    printf("%2d   %20.12"LZP"f   %20.12"LZP"f   %20.12"LZP"f\n",
           i, x[i], y[i], w[i]);


  /* ------ compute coefficients of least square ------------------- */
  /* ------ approximating polynomial ------------------------------- */

  fehler = pol_appr(grad, stuetz - 1, x, y, w, c, b, d);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("pol_appr(): Polynomial degree exceeds number of "
                    "nodes or\n                          "
                    "less than two nodes or negative pol. degree\n",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("pol_appr(): When interpolating two x[i] "
                    "are identical", 10 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("pol_appr(): not all weights are positive",
                     10 + fehler, __FILE__, __LINE__);
    default:
      fehler_melden("pol_appr(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* -------- print coefficients of solution polynomial ------------- */

  printf("\n\n"
         "coefficients of least square polynomial:\n\n"
         " i            c[i]                   b[i]                   "
         "d[i]\n"
         "------------------------------------------------------------"
         "-----------\n");
  for (i = 0; i <= grad; i++)
    if (i == 0)
      printf("%2d  %21.12"LZP"f"
             "             -                     -\n", i, c[i]);
    else if (i == 1)
      printf("%2d  %21.12"LZP"f  %21.12"LZP"f            -\n",
             i, c[i], b[i]);
    else
      printf("%2d  %21.12"LZP"f  %21.12"LZP"f  %21.12"LZP"f\n",
             i, c[i], b[i], d[i]);


  /* --- print standard coefficients a of approximating polynomial -- */

  fehler = opolkoeff(grad, b, d, c, a);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("opolkoeff(): degree of polynomial negative",
                    20 + fehler, __FILE__, __LINE__);
    case 3:
      fehler_melden("opolkoeff(): lack of memory",
                    20 + fehler, __FILE__, __LINE__);
    default:
      fehler_melden("opolkoeff(): unknown error",
                    20 + fehler, __FILE__, __LINE__);
  }


  /* -------- print a-coefficients of approximating polynomial ------ */

#ifdef DEBUG
  printf("\n\n"
         "standard coefficients a of approximating polynomial in "
         "representation\n"
         "a[0] + a[1]*x + a[2]*x^2 + ... + a[n]*x^n:\n\n"
         " i            a[i]\n"
         "-------------------------\n"
         );
  for (i = 0; i <= grad; i++)
    printf("%2d  %21.12"LZP"f\n", i, a[i]);
#endif


  /* ------ compute table of values for polynomial ------------------ */

  printf("\n\n"
         "table of values for polynomial:\n\n"
         "  i       x             f(x)\n"
         "--------------------------------\n");

  for (i = 0; i < MAXTLEN; i++)
  {
    xtab[i] = (REAL)i / (MAXTLEN - 1);
    ytab[i] = opolwert(grad, xtab[i], b, d, c);
    printf("%3d   %12.6"LZP"e  %12.6"LZP"e\n", i, xtab[i], ytab[i]);
  }


  /* ---------------- plot polynomial if desired -------------------- */

  if (argc <= 3 || *argv[3] != 'n')       /* graphics not suppressed? */
  {
    fehler = zeigkurv(MAXTLEN, stuetz, xtab, ytab, x, y);
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


  /* ------ compute table of values for polynomial,             ----- */
  /* ------ but now based on the a-cofficients                  ----- */

#ifdef DEBUG
  printf("\n\n"
         "table of values for polynomial (based on "
         "a-coefficients):\n\n"
         "  i       x             f(x)\n"
         "--------------------------------\n");
#endif

  for (i = 0; i < MAXTLEN; i++)
  {
    xtab[i] = (REAL)i / (MAXTLEN - 1);
    ytab[i] = horner(grad, a, xtab[i]);
#ifdef DEBUG
    printf("%3d   %12.6"LZP"e  %12.6"LZP"e\n", i, xtab[i], ytab[i]);
#endif
  }


  return 0;
}
