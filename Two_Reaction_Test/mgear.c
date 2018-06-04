#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE mgear.c ------------------------- */
#ifndef VOLLTEST                            /* Test with input file ? */

#include <basis.h>         /*  for umleiten, fprintf, stderr, scanf, */
                           /*       printf, NULL, REAL, LZS, LZP,     */
                           /*       fehler_melden, fehler_t           */
#include <vmblock.h>       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include <gear.h>          /*  for  gear4, gear_fehlertext            */
#include <t_dgls.h>        /*  for  bsptyp, dgls_waehlen              */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
* Testprogramm for the function gear4() .                              *
*                                                                      *
* Mode of operation:                                                   *
* ==================                                                   *
* Standard input and output fiels are stdin and stdout. The first enrty*
* of the command line is taken to be the input file, while the second  *
* one serves as output file, if it exists. Calls for input and error   *
* messages are in stderr.                                              *
* Following reading of the inoput, it is printed out with the computed *
* results.                                                             *
*
* To test other systems of DEs, please proceed as explained in t_dgls.c*
*                                                                      *
* Input file :                                                         *
* ============                                                         *
* bspnummer  Number of DE system from t_dgls.c                         *
* epsabs     desired absolute error bound                              *
* epsrel     desired relative error bound                              *
* x0         left edge of integration                                  *
* y0[0]   \  known approximation for the solution at x0                *
* ...      >                                                           *
* y0[n-1] /                                                            *
* h          initial step size                                         *
* xend       right endpoint of integration                             *
* fmax       maximal number of calls of the right hand side            *
*                                                                      *
* The size n of the DE system is passed on from t_dgls.c.              *
***********************************************************************/

{
  REAL     epsabs;       /* absolute error bound                      */
  REAL     epsrel;       /* relative error bound                      */
  REAL     x0;           /* left edge of integration interval         */
  REAL     *y0;          /* [0..n-1]-vector: initial value, approxim. */
  REAL     *yex;         /* [0..n-1]-vector: exact solution           */
  REAL     h;            /* initial, final step size                  */
  REAL     xend;         /* right edge of integration interval        */
  long     fmax;         /* maximal number of calls of right side     */
                         /* in gear4()                                */
  long     aufrufe;      /* actual number of function calls           */
  int      bspnummer;    /* Number of the system of DEs from t_dgls.c */
  int      n;            /* number of DEs in system                   */
  int      fehler;       /* error code from umleiten(), gear4()       */
  int      i;            /* loop counter                              */
  fehler_t fehlerart;    /* error classification for gear4()          */
  bsptyp   *beispiel;    /* pointer to the structure that describes   */
                         /* the actual system of DEs                  */
  void     *vmblock;     /* List of dynamically allocated vectors     */


  if ((fehler = umleiten(argc, argv))/* assign input/output files     */
      != 0)                          /* to standard ones              */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input  -------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "Test example number:                             ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = dgls_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("non registered example",
                  0, __FILE__, __LINE__);
    return 3;
  }

  n = beispiel->n;
  vmblock = vminit();                 /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  yex = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))   /* out of memory? */
  {
    fehler_melden("out of memory", 0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "absolute error bound epsabs:                 ");
#endif
  scanf("%"LZS"f", &epsabs);

#ifdef INTERAKTIV
  fprintf(stderr, "relative error bound epsrel:                 ");
#endif
  scanf("%"LZS"f", &epsrel);

#ifdef INTERAKTIV
  fprintf(stderr, "left edge x0:      ");
#endif
  scanf("%"LZS"f", &x0);

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "Function value y0[%d] at x0:               ",
                    i);
#endif
    scanf("%"LZS"f", y0 + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "initial step size h:                    ");
#endif
  scanf("%"LZS"f", &h);

#ifdef INTERAKTIV
  fprintf(stderr, "right edge  xend:   ");
#endif
  scanf("%"LZS"f", &xend);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of calls of right hand side: ");
#endif
  scanf("%ld", &fmax);


  /* ------------ put out the input data ----------- */

  printf("\n"
         "Solve a first order ordinary system of DEs\n"
         "==========================================\n"
         "using the implicit method of Gear of 4th order\n"
         "==============================================\n\n\n"
         "System of DEs:\n"
         "--------------\n"
         "%s\n\n"
         "Input data:\n"
         "-----------\n"
         "Example  = %24d\n"
         "n        = %24d\n"
         "x0       = %24.15"LZP"e\n"
         "xend     = %24.15"LZP"e\n"
         "epsabs   = %24.15"LZP"e\n"
         "epsrel   = %24.15"LZP"e\n"
         "fmax     = %24ld\n"
         "h        = %24.15"LZP"e\n",
         (*beispiel->dgl_text)(), bspnummer, n, x0, xend, epsabs,
         epsrel, fmax, h);

  for (i = 0; i < n; i++)
    printf("y0[%d]    = %24.15"LZP"e\n", i, y0[i]);


  /* ------------ Solve system of DEs -------------- */

  fehler = gear4(&x0, xend, n, beispiel->rechte_seite, y0, epsabs,
                 epsrel, &h, fmax, &aufrufe);

  if (fehler != 0)
  {
    fehler_melden(gear_fehlertext(fehler, &fehlerart), 10 + fehler,
                  __FILE__, __LINE__);
    if (fehlerart > WARNUNG)
      return 10 + fehler;
  }


  /* -------------------- put out results ------------------- */

  printf("\n\n"
         "Output data:\n"
         "------------\n"
         "error code from gear4():                  %24d\n"
         "final local step size:                    %24.15"LZP"e\n"
         "number of calls of right hand side:       %24ld\n"
         "Integration stopped at x =                %24.15"LZP"e\n\n",
         fehler, h, aufrufe, x0);

  for (i = 0; i < n; i++)
    printf("approximate solution y%d(x) = %24.15"LZP"e\n",
           i + 1, y0[i]);

  if (beispiel->exakte_loesung != NULL)       /* "exact" solution     */
  {                                           /* available?           */
    (*beispiel->exakte_loesung)(x0, yex);
    printf("\n");
    for (i = 0; i < n; i++)
      printf("'exact' solution     y%d(x) = %24.15"LZP"e\n",
             i + 1, yex[i]);
    printf("\nDifference  approximate solution - 'exact' solution:\n");
    for (i = 0; i < n; i++)
      printf("%24.15"LZP"g\n", y0[i] - yex[i]);
  }


  return 0;
}
#else

#include <basis.h>         /*  for  REAL                              */
#include <gear.h>          /*  for  gear4                             */
#include <t_dgls.h>        /*  for  bsptyp, dgls_waehlen              */



/* ------------------------------------------------------------------ */

int main(void)

/***********************************************************************
* Program to test gear4() with the stiff IVPs 5 (default) or 11 (if    *
* macro FNUM is defined before compilation) from `t_dgls.c'.           *
* We use various values for the absolute and relative error bounds.    *
***********************************************************************/

{
  static
    REAL eps[4][2] = {{ (REAL)1.0e-6,  (REAL)1.0e-6  },
                      { (REAL)1.0e-10, (REAL)1.0e-10 },
                      { (REAL)1.0e-11, ZERO          },
                      { ZERO,          (REAL)1.0e-11 }
                     };
  REAL   y[2];
  REAL   yex[2];       /* [0..n-1] vector: exact solution             */
  int    k;
  REAL   x;
  REAL   h;
  bsptyp  *beispiel;   /* pointer to the structure with the DE system */
  int    n;
  long   aufrufe;
  int    i;
  REAL   xx;
  int    fehler;


#ifdef FNUM
  if ((beispiel = dgls_waehlen(11)) == NULL)
#else
  if ((beispiel = dgls_waehlen(5)) == NULL)
#endif
  {
    fehler_melden("non registered example",
                  0, __FILE__, __LINE__);
    return 3;
  }

  n = beispiel->n;

  for (k = 0; k < 4; k++)
  {
    x    = ZERO;
#ifdef FNUM
    y[0] =  ONE;
    y[1] = -ONE;
#else
    y[0] = FOUR;
    y[1] = TWO;
#endif
    h    =  (REAL)0.01;
    printf("epsabs = %8"LZP"g  epsrel = %8"LZP"g\n",
           eps[k][0], eps[k][1]);
    for (i = 5; i <= 50; i += 5)
    {
      xx = (REAL)0.1 * (REAL)i;
      fehler = gear4(&x, xx, n, beispiel->rechte_seite, y, eps[k][0],
                    eps[k][1], &h, 50000l, &aufrufe);
#ifndef FNUM
      if (beispiel->exakte_loesung != NULL)
        (*beispiel->exakte_loesung)(x, yex);
      else
        printf("There is no analytical solution!\n");
#endif
      printf("error code = %d  # calls = %ld\n", fehler, aufrufe);
      printf("x = %4.1"LZP"f  y = %13.9"LZP"g  computational error = "
             "%16.9"LZP"g\n", x, y[0],
#ifdef FNUM
             FABS(y[0] - EXP(-x)));
#else
             FABS(y[0] - yex[0]) + FABS(y[1] - yex[1]));
#endif
    }
    printf("\n");
  }


  return 0;
}
#endif

/* --------------------------- END mgear.c -------------------------- */
