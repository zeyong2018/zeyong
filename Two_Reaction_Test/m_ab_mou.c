#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>         /*  for  umleiten, fprintf, stderr, scanf, */
                           /*       printf, NULL, FABS, REAL, LZS,    */
                           /*       fehler_melden, LZP, TRUE          */
#include <vmblock.h>       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include <ab_mou.h>        /*  for  prae_korr                         */
#include <t_dgls.h>        /*  for  bsptyp, dgls_waehlen              */
#ifdef __TURBOC__
#ifdef __MSDOS__
#include <alloc.h>         /*  for  coreleft                          */
#endif
#endif



/***********************************************************************
*                                                                      *
* Test program for the function prae_korr() from the module ab_mou     *
* to solve a first order system of DEs using the multi-step method of  *
* Adams-Bashforth-Moulton.                                             *
*                                                                      *
* Scope of program:                                                    *
* =================                                                    *
* The program reads the input data from the file stdin and writes      *
* output onto the file stdout.                                         *
* if the first command line entry exists, it is interpreted to be the  *
* input file and this is assigned to stdin.                            *
* The second command line parameter is treated analogously wrt. stdout.*
* Calls for input are directed to stderr which also collects error     *
* messages.                                                            *
*                                                                      *
* After reading the input, it is put out for control purposes; then    *
* the multi-step method is executed and the computed results are put   *
* out.                                                                 *
*                                                                      *
* To solve a differential equation, please proceed as in example       *
* t_dgls.c .                                                           *
*                                                                      *
* Construction of input data files:                                    *
* =================================                                    *
* bspnummer  Number of DE system in  t_dgls.c                          *
* epsabs     desired absolute error bound                              *
* epsrel     desired relative error bound                              *
* x0         initial x-value                                           *
* y0[0]   \  initial y-value at x0                                     *
* ...      >                                                           *
* y0[n-1] /                                                            *
* h          starting step size                                        *
* hmax       maximal step size                                         *
* fmax       maximal number of calls of right hand side                *
* xend       final desired x-value                                     *
*                                                                      *
* The number n of differential equations follows from the number of    *
* the DE system chosen; it is stored in conjunction with the right     *
* hand side function in t_dgls.c .                                     *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL    epsabs,
          epsrel,
          x0,
          *y0,
          h,
          hmax,
          xend;
  long    fmax,
          aufrufe,
          j;
  int     bspnummer,
          n,
          fehler,        /* error code from prae_korr()               */
          i;             /* loop counter                              */
  bsptyp  *beispiel;     /* pointer to the stucture describing the DE */
                         /* system                                    */
  void    *vmblock;      /* List of dynamic allocations               */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard ones           */
    return fehler;                                          /* 1 or 2 */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "Example :                             ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = dgls_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("improper example number",
                  0, __FILE__, __LINE__);
    return 3;
  }

  n = beispiel->n;
  vmblock = vminit();                          /* initialize storage */
  y0 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "absolute error bound epsabs:          ");
#endif
  scanf("%"LZS"f", &epsabs);

#ifdef INTERAKTIV
  fprintf(stderr, "relative error bound epsrel:          ");
#endif
  scanf("%"LZS"f", &epsrel);

#ifdef INTERAKTIV
  fprintf(stderr, "initial x-value x0:                   ");
#endif
  scanf("%"LZS"f", &x0);

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "initial y-value y0[%d] at x0:          ",
            i);
#endif
    scanf("%"LZS"f", y0 + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "starting step size h:                 ");
#endif
  scanf("%"LZS"f", &h);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal step size hmax:               ");
#endif
  scanf("%"LZS"f", &hmax);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of calls of r.h. side: ");
#endif
  scanf("%ld", &fmax);

#ifdef INTERAKTIV
  fprintf(stderr, "desired final x-value :               ");
#endif
  scanf("%"LZS"f", &xend);


  /* ------------ print input data --------------------------------- */

  printf("\n"
         "Solution of a first order ordinary system of DEs\n"
         "================================================\n"
         "using the multi-step method of Adams-Bashforth-Moulton\n"
         "======================================================\n"
         "\n\n"
         "DE system:\n"
         "----------\n"
         "%s\n\n"
         "Input data:\n"
         "-----------\n"
         "Example  = %24d\n"
         "n        = %24d\n"
         "epsabs   = %24.15"LZP"e\n"
         "epsrel   = %24.15"LZP"e\n"
         "x0       = %24.15"LZP"e\n"
         "xend     = %24.15"LZP"e\n"
         "h        = %24.15"LZP"e\n",
         (*beispiel->dgl_text)(), bspnummer, n, epsabs, epsrel,
         x0, xend, h);

  for (i = 0; i < n; i++)
    printf("y0[%d]    = %24.15"LZP"e\n", i, y0[i]);


  /* ------------ Solve DE system ---------------------------------- */
#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "\nfree before:  %u\n", coreleft());
#endif
#endif
#endif

  fehler = prae_korr(&x0, y0, n, beispiel->rechte_seite, xend, &h,
                     epsabs, epsrel, fmax, &aufrufe, hmax, TRUE);

  if (fehler != 0)
  {
    fehler_melden("prae_korr()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }

#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "free after: %u\n", coreleft());
#endif
#endif
#endif
  /*  The following call of "prae_korr()" is done to free allocations */
  /*  We achieve this by setting n = 0 illegally here.                */
  prae_korr(&x0, y0, 0, beispiel->rechte_seite, xend, &h, epsabs,
            epsrel, fmax, &j, hmax, TRUE);
#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "free after: %u\n", coreleft());
#endif
#endif
#endif


  /* -------------------- put out the results ---------------------- */

  printf("\n\n"
         "Output data:\n"
         "------------\n"
         "error code from prae_korr():       %24d\n"
         "final local step size:             %24.15"LZP"e\n"
         "final x-value for integration x =  %24.15"LZP"e\n\n",
         fehler, h, x0);

  for (i = 0; i < n; i++)
    printf("approximate solution y%d(x) = %24.15"LZP"e\n",
           i + 1, y0[i]);

  if (beispiel->exakte_loesung != NULL)  /* exact solution available? */
  {
    (*beispiel->exakte_loesung)(x0, y0);
    printf("\n");
    for (i = 0; i < n; i++)
      printf("'exact' solution  y%d(x)    = %24.15"LZP"e\n",
             i + 1, y0[i]);
  }


  return 0;
}
