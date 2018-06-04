#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, NULL,      */
                         /*       fprintf, stderr, REAL, LZS, LZP,    */
                         /*       ZERO, fehler_melden                 */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <fft.h>         /*  for  nli_happr                           */
#include <nlappphi.h>    /*  for  bsptyp3, nliansf_waehlen            */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define MAXTLEN  80         /* maximal length of the table of values  */



/***********************************************************************
* This is a test program for  nli_happr()  for finding the optimal     *
* coefficients of a nonlinear model function for a least square problem*
*                                                                      *
* Scope of program:                                                    *
* =================                                                    *
* The program reads the input from the first file on the command line. *
* Then it computes the least square approximation and tabulates its    *
* function values. The output is sent to the second named file on the  *
* command line. If there is none, all output appears on the screen.    *
* If there is no input file either, the input is to be given from the  *
* keyboard.                                                            *
* Several test examples are available in the files  nlapptst.ei*.      *
* Their output should coincide with the data in  nlapptst.au* .        *
* The model function and its derivatives wrt. its coefficients must be *
* supplied in nlappphi.c.                                              *
* For Turbo C compilers we plot the result as well as the original     *
* starting function. This can be avoided by adding a third entry of    *
* "n" on the command line.                                             *
*                                                                      *
* Construction of input file:                                          *
* ===========================                                          *
* bspnummer             Number of the model function in  nlappphi.c    *
* m                     number of nodes                                *
* x[0],y[0],w[0]        nodes, function values, weights                *
* x[1],y[1],w[1]        ...                                            *
* ...                   ...                                            *
* x[m-1],y[m-1],w[m-1]  ...                                            *
* maxit                 maximal number of Newton iterations            *
* c[0]   \                                                             *
* ...     >             Starting values for optimal coefficients       *
* c[n-1] /                                                             *
*                                                                      *
* The number n of coefficients in the model function is part of the    *
* input from  nlappphi.c.                                              *
* muss.                                                                *
***********************************************************************/

int main(int argc, char *argv[])

{
  int    n,            /* number of coefficients in model function    */
         m,            /* number of nodes                             */
         i,            /* Loop variable                               */
         fehler,       /* error code of  nli_happr()                  */
         maxit;        /* maximal number of Newton iterations for     */
                       /* nli_happr() / actual number of iterations   */
  unsigned int
          bspnummer;   /* Number of example for model function        */
  REAL    *x,          /* node vector                                 */
          *y,          /* function value vector                       */
          *w,          /* vector of weights                           */
          *c,          /* coefficient vector for least square solution*/
          tabanf,      /* left end point of table of values           */
          tabend,      /* right end point ...                         */
          xx,          /* Loop variable                               */
          schritt,     /* step size for table                         */
          *xtab,       /* Table of values for solution: x-values      */
          *ytab,       /*                               y-values      */
          mqf;         /* mean square error                           */
  void    *vmblock;    /* List of dynamically allocated vectors       */
  bsptyp3 *beispiel;   /* pointer to the structure that describes     */
                       /* the nonlinear model function                */


  if ((fehler = umleiten(argc, argv))   /* assign the input/output    */
      != 0)                             /* files to the standard      */
    return fehler;  /* 1 or 2 */        /* files                      */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nExample:                                      ");
#endif
  scanf("%u", &bspnummer);
  if ((beispiel = nliansf_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("non registered example",
                  0, __FILE__, __LINE__);
    return 3;
  }
  n = beispiel->n;

#ifdef INTERAKTIV
  fprintf(stderr, "number of nodes:                              ");
#endif
  scanf("%d", &m);

  vmblock = vminit();                       /* initialize storage     */
  x    = (REAL *)vmalloc(vmblock, VEKTOR, m,       0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, m,       0);
  w    = (REAL *)vmalloc(vmblock, VEKTOR, m,       0);
  c    = (REAL *)vmalloc(vmblock, VEKTOR, n,       0);
  xtab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN, 0);
  ytab = (REAL *)vmalloc(vmblock, VEKTOR, MAXTLEN, 0);
  if (! vmcomplete(vmblock))                /* lack of memory ?       */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

  for (i = 0; i < m; i++)               /* read nodes and weights     */
  {
#ifdef INTERAKTIV
    fprintf(stderr, "x-value, y-value, weight:                   ");
#endif
    scanf("%"LZS"f%"LZS"f%"LZS"f", x + i, y + i, w + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of iterations:                 ");
#endif
  scanf("%d", &maxit);               /* read maxit                    */

  for (i = 0; i < n; i++)            /* read starting values          */
  {
#ifdef INTERAKTIV
  fprintf(stderr, "Starting values for the optimal coefficients: ");
#endif
    scanf("%"LZS"f", c + i);
  }


  /* ------------------ print out the input data ------------------- */

  printf("\n"
         "nonlinear approximation using "
         "Householder transformations\n"
         "=============================="
         "===========================\n\n\n"
         "number of coefficients in the model function: %5d\n"
         "number of nodes:                              %5d\n"
         "maximal number of Newton steps:               %5d\n\n"
         "model function:\n"
         "%s\n\n"
         "Table of x-values, y-values and weights :\n\n"
         " i            x[i]                   y[i]                   "
         "w[i]\n"
         "------------------------------------------------------------"
         "-----------\n",
         beispiel->n, m, maxit, beispiel->PHI_text()
        );

  for (i = 0; i < m; i++)
    printf("%2d   %20.12"LZP"f   %20.12"LZP"f   %20.12"LZP"f\n",
           i, x[i], y[i], w[i]);

  printf("\n\n"
         "Starting vector for coeffficients:\n\n"
         " i           c[i]\n"
         "-------------------------\n");
  for (i = 0; i < n; i++)
    printf("%2d   %20.12"LZP"f\n", i, c[i]);


  /* -------- set up x-values for table with the solution ----------- */

  schritt = x[m - 1] - x[0];
  tabanf  = x[0]     - ZERO * schritt;
  tabend  = x[m - 1] + ZERO * schritt;
  schritt = (tabend - tabanf) / (MAXTLEN - 1);

  for (xx = tabanf, i = 0; i < MAXTLEN; i++, xx += schritt)
    xtab[i] = xx,
    ytab[i] = (*beispiel->PHI)(c, xx);


  /* ------ plot model function for starting coefficients ----------- */
  /* ------ if desired                                    ----------- */

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


  /* ------ compute optimal coefficients for this example ----------- */

  fehler = nli_happr(m, n - 1, x, y, w, beispiel->PHI, TRUE,
                     beispiel->ABL, &maxit, (REAL)0.0005, c, &mqf);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("nli_happr(): m <= n or n < 1",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("nli_happr(): error in Householder transform",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("nli_happr(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 4:
      fehler_melden("nli_happr(): max. iteration number exceeded",
                    10 + fehler, __FILE__, __LINE__);
    default:
      fehler_melden("nli_happr(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ------ record result of nonlinear approximation --------------- */

  printf("\n\n"
         "coefficients of the optimal function PHI:\n\n"
         " i           c[i]\n"
         "-------------------------\n");
  for (i = 0; i < n; i++)
    printf("%2d   %20.12"LZP"f\n", i, c[i]);

  printf("\n"
         "mean square error of PHI:         %7.3"LZP"f\n"
         "number of Newton steps performed: %7d\n", mqf, maxit);


  /* ----------- compute table of values for optimal solution ------- */

  printf("\n\n"
         "table of values for optimal solution PHI:\n\n"
         "  i           x                   PHI(x)\n"
         "---------------------------------------------\n");

  for (i = 0; i < MAXTLEN; i++)
  {
    ytab[i] = (beispiel->PHI)(c, xtab[i]);
    printf("%3d %20.12"LZP"f %20.12"LZP"f\n", i, xtab[i], ytab[i]);
  }


  /* ----------- plot optimal solution, if desired ------------------ */

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
        fehler_melden("zeigkurv(): BGI graphic error",
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
