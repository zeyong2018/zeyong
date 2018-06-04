#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE m_awp.c ------------------------- */

#include <basis.h>         /*  for  umleiten, fprintf, stderr, scanf, */
                           /*       printf, NULL, REAL, LZS, LZP,     */
                           /*       fehler_melden                     */
#include <vmblock.h>       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include <awp.h>           /*  for  awp, fehler_t, awp_fehlertext     */
#include <t_dgls.h>        /*  for  bsptyp, dgls_waehlen              */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
* Test program for the function awp() from the module awp to solve     *
* a first order ordinary system of differential equations with         *
* automatic step size control.                                         *
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
* the adaptive initial value solver is executed and the computed       *
* results are put out.                                                 *
*                                                                      *
* To solve a differential equation, please proceed as in example       *
* t_dgls.c.                                                            *
*                                                                      *
* Construction of input data files:                                    *
* =================================                                    *
* bspnummer  Number of DE system in  t_dgls.c                          *
* epsabs     desired absolute error bound                              *
* epsrel     desired relative error bound                              *
* methode    desired method for awp():                                 *
*            = 3 for Runge-Kutta embedding formula of order  2 and 3,  *
*            = 6 for the embedding formula of England of order 4 and 5,*
*            = 7 for the embedding formula of Prince and Dormand  of   *
*                order 4 and 5, combined with test for stiffness.      *
* x0         initial x-value                                           *
* y0[0]   \  initial y-value at x0                                     *
* ...      >                                                           *
* y0[n-1] /                                                            *
* h          starting step size                                        *
* fmax       maximal number of calls of right hand side                *
* xend       final desired x-value                                     *
*                                                                      *
* The number n of differential equations follows from the number of    *
* the DE system chosen; it is stored in conjunction with the right     *
* hand side function in t_dgls.c.                                      *
***********************************************************************/

{
  REAL   epsabs,
         epsrel,
         x0,
         *y0,
         *yex,           /* [0..n-1] vector: exact solution           */
         h,
         xend;
  long   fmax,         /* maximal number of calls of right hand side  */
                       /* of DE in awp()                              */
         aufrufe;      /* actual number of calls of r h side in awp() */
  int    bspnummer,
         methode,
         n,
         fehler,         /* error code from awp()                     */
         i;              /* loop counter                              */
  fehler_t fehlerart;    /* error classification   for awp()          */
  bsptyp *beispiel;      /* pointer to the structure describing the   */
                         /* problem at hand                           */
  void   *vmblock;       /* List of dynamic allocations               */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard files          */
    return fehler;                                          /* 1 or 2 */


  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nExample:                               ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = dgls_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("improper example number",
                  0, __FILE__, __LINE__);
    return 3;
  }

  n = beispiel->n;
  vmblock = vminit();                           /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  yex = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "absolute error bound epsabs:           ");
#endif
  scanf("%"LZS"f", &epsabs);

#ifdef INTERAKTIV
  fprintf(stderr, "relative error bound epsrel:           ");
#endif
  scanf("%"LZS"f", &epsrel);

#ifdef INTERAKTIV
  fprintf(stderr, "Method (3 or 6):                       ");
#endif
  scanf("%d", &methode);

#ifdef INTERAKTIV
  fprintf(stderr, "initial x-value x0:                    ");
#endif
  scanf("%"LZS"f", &x0);

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "Function value y0[%d] at x0:            ",
                    i);
#endif
    scanf("%"LZS"f", y0 + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "starting step size h:                  ");
#endif
  scanf("%"LZS"f", &h);

#ifdef INTERAKTIV
  fprintf(stderr, "desired final x-value:                 ");
#endif
  scanf("%"LZS"f", &xend);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of calls for r.h. side: ");
#endif
  scanf("%ld", &fmax);


  /* ------------ print out input ----------------------------------- */

  printf("\n"
         "Solution of a first order ordinary differential equations "
         "system\n"
         "=========================================================="
         "======\n"
         "with automatic step size control\n"
         "================================\n\n\n"
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
         "methode  = %24d\n"
         "h        = %24.15"LZP"e\n",
         (*beispiel->dgl_text)(), bspnummer, n, x0, xend, epsabs,
         epsrel, fmax, methode, h);

  for (i = 0; i < n; i++)
    printf("y0[%d]    = %24.15"LZP"e\n", i, y0[i]);


  /* ------------ solve DE system ----------------------------------- */

  fehler = awp(&x0, xend, n, beispiel->rechte_seite, y0, epsabs,
               epsrel, &h, methode, fmax, &aufrufe);

  if (fehler != 0)
  {
    fehler_melden(awp_fehlertext(fehler, &fehlerart), 10 + fehler,
                  __FILE__, __LINE__);
    if (fehlerart > WARNUNG)
      return 10 + fehler;
  }


  /* -------------------- print results ----------------------------- */

  printf("\n\n"
         "Output data:\n"
         "------------\n"
         "error code from awp():              %24d\n"
         "final step size:                    %24.15"LZP"e\n"
         "number of evaluations of r.h. side: %24ld\n"
         "final x-value x =                   %24.15"LZP"e\n\n",
         fehler, h, aufrufe, x0);

  for (i = 0; i < n; i++)
    printf("approximate solution y%d(x) = %24.15"LZP"e\n",
           i + 1, y0[i]);

  if (beispiel->exakte_loesung != NULL)  /* exact solution available? */
  {
    (*beispiel->exakte_loesung)(x0, yex);
    printf("\n");
    for (i = 0; i < n; i++)
      printf("'exact' solution  y%d(x)    = %24.15"LZP"e\n",
             i + 1, yex[i]);
    printf("\ndifference  approximate solution - 'exact' solution:\n");
    for (i = 0; i < n; i++)
      printf("%24.15"LZP"g\n", y0[i] - yex[i]);
  }


  return 0;
}

/* -------------------------- END  m_awp.c -------------------------- */
