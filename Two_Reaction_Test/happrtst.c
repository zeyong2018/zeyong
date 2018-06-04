#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, NULL,      */
                         /*       fprintf, stderr, REAL, LZS, LZP,    */
                         /*       fehler_melden                       */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <approx.h>      /*  for  lin_happr, lin_hwert                */
#include <happrans.h>    /*  for  bsptyp2, linansf_waehlen            */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define MAXTLEN  80         /* maximal size of table of values        */



/***********************************************************************
* This is a test program for the program  lin_happr()  to compute the  *
* optimal coefficients of a linear sum of n user given functions using *
* Householder transformations.                                         *
*                                                                      *
* Scope of program:                                                    *
* =================                                                    *
* The program reads the input fromthe first file in the command line.  *
* Then it computes the linear least square approximation.              *
* Using  lin_hwert() it constructs a table of values for the polynomial*
* Output is assigned to the second named file. If there is none, output*
* appears on the screen. If the first file is also missing, the input  *
* must be entered from the keyboard.                                   *
* Several model problems are in the files happrtst.ei*. The output     *
* should be the same as that in happrtst.au* .                         *
* The model functions must be defined inside happrans.c .              *
* For Turbo C compiler the polynomial will be plotted at its nodes,    *
* except when the third command line input consists of "n".            *
*                                                                      *
* Construction of an input file:                                       *
* ==============================                                       *
* bspnummer             Number of model function in happrans.c         *
* n                     number of functions used (only if this can be  *
*                       user specified)                                *
* m                     number of nodes                                *
* x[0],y[0],w[0]        nodes, function values, weights                *
* x[1],y[1],w[1]        ditto                                          *
* ...                   ...                                            *
* x[m-1],y[m-1],w[m-1]  ditto                                          *
***********************************************************************/

int main(int argc, char *argv[])

{
  int     m,           /* number of nodes                             */
          n,           /* number of test functions in system          */
          i,           /* Loop variable                               */
          fehler;      /* error code for lin_happr()                  */
  unsigned int
          bspnummer;   /* Number of example                           */
  REAL    *x,          /* vector of nodes                             */
          *y,          /* vector of functional values                 */
          *w,          /* vector of weights                           */
          *c,          /* coefficient vector for approx. polynomial   */
          tabanf,      /* left end point of tabulated interval        */
          tabend,      /* right end point of tabulated interval       */
          xx,          /* Loop variable for tabulation                */
          schritt,     /* atep size of tabulation                     */
          *xtab,       /* table of values : x-values                  */
          *ytab,       /*                   y-values                  */
          mqf;         /* mean least square error of approximation    */
  void    *vmblock;    /* List of dynamically allocated vectors       */
  bsptyp2 *beispiel;   /* pointer to structure that describes the     */
                       /* system of model functions                   */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nExample:                                ");
#endif
  scanf("%u", &bspnummer);
  if ((beispiel = linansf_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("Example number not registered",
                  0, __FILE__, __LINE__);
    return 3;
  }

#ifdef INTERAKTIV
  if (beispiel->n == 0)  /* free choice of number of model functions? */
    fprintf(stderr, "Number of model functions (1,2,...):  ");
#endif
  scanf("%d", &n);

  if (beispiel->n > 0 && n != beispiel->n)        /* improper number? */
  {
    fehler_melden("too few or too many model functions",
                  0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "number of nodes:                        ");
#endif
  scanf("%d", &m);

  vmblock = vminit();                             /* allocate storage */
  x    = (REAL *)vmalloc(vmblock, VEKTOR, m,       0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, m,       0);
  w    = (REAL *)vmalloc(vmblock, VEKTOR, m,       0);
  c    = (REAL *)vmalloc(vmblock, VEKTOR, n,       0);
  xtab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN, 0);
  ytab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN, 0);
  if (! vmcomplete(vmblock))               /* allocations all right ? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 5;
  }

  for (i = 0; i < m; i++)              /* read nodes, function values */
  {                                    /* and weights                 */
#ifdef INTERAKTIV
    fprintf(stderr, "node, function value, weight:           ");
#endif
    scanf("%"LZS"f%"LZS"f%"LZS"f", x + i, y + i, w + i);
  }


  /* ------------------- print input -------------------------------- */

  printf("\n"
         "general linear approximation using "
         "Householder transformations\n"
         "==================================="
         "===========================\n\n\n"
         "number of model functions: %5d\n"
         "number of nodes:           %5d\n\n"
         "system of functions:\n"
         "%s\n\n"
         "Table of values for: nodes, function values, weights:\n\n"
         " i            x[i]                   y[i]                   "
         "w[i]\n"
         "------------------------------------------------------------"
         "-----------\n",
         n, m, (*beispiel->phi_text)()
        );

  for (i = 0; i < m; i++)
    printf("%2d   %20.12"LZP"f   %20.12"LZP"f   %20.12"LZP"f\n",
           i, x[i], y[i], w[i]);


  /* ------- compute approximating function ------------------------- */

  fehler = lin_happr(m, n - 1, x, y, w, beispiel->phi, c, &mqf);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("lin_happr(): m <= n or n < 1",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("lin_happr(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("lin_happr(): model functions linearly dependent",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("lin_happr(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* -------- put out results --------------------------------------- */

  printf("\n\n"
         "coefficients of the approximating function PHI:\n\n"
         " i           c[i]\n"
         "-------------------------\n");
  for (i = 0; i < n; i++)
    printf("%2d   %20.12"LZP"f\n", i, c[i]);

  printf("\nmean square error of PHI: "
         "%7.3"LZP"f\n", mqf);


  /* ------ generate a table of values for PHI ---------------------- */

  schritt = x[m - 1] - x[0];
  tabanf = x[0]     - (REAL)0.05 * schritt;
  tabend = x[m - 1] + (REAL)0.05 * schritt;
  schritt = (tabend - tabanf) / (MAXTLEN - 1);

  printf("\n\n"
         "table of values for PHI:\n\n"
         "  i           x                   PHI(x)\n"
         "---------------------------------------------\n");

  for (xx = tabanf, i = 0; i < MAXTLEN; i++, xx += schritt)
  {
    xtab[i] = xx;
    ytab[i] = lin_hwert(xx, n - 1, beispiel->phi, c);
    printf("%3d %20.12"LZP"f %20.12"LZP"f\n",
           i, xtab[i], ytab[i]);
  }


  /* ----------- show the spline function for the optimal  ---------- */
  /* ----------- coefficients, if desired                  ---------- */

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
