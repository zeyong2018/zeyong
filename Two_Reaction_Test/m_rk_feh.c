#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>              /*  for  umleiten, fprintf, stderr,   */
                                /*       scanf, printf, NULL, REAL,   */
                                /*       LZS, LZP, fehler_melden      */
#include <vmblock.h>            /*  for  vmalloc, vmcomplete, vmfree, */
                                /*       vminit, VEKTOR               */
#include <rk_fehl.h>            /*  for  rk_fehl                      */
#include <t_dgls.h>             /*  for  bsptyp, dgls_waehlen         */



/***********************************************************************
*                                                                      *
* Test program for the function rk_fehl() from the module  rk_fehl to  *
* solve a first order ordinary system of differential equations using  *
* the Runge-Kutta-Fehlberg method of order O(h^6) with local error     *
* estimation and step size control.                                    *
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
* the initial value solver is executed and the computed results are    *
* put out.                                                             *
*                                                                      *
* To solve a differential equation, please proceed as in example       *
* t_dgls.c.                                                            *
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
* xdiff      Length of interval of integration                         *
*                                                                      *
* The number n of differential equations follows from the number of    *
* the DE system chosen; it is stored in conjunction with the right     *
* hand side function in t_dgls.c.                                      *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL   x0,
         *y0,
         xdiff,
         h,
         hmax,
         epsabs,
         epsrel;
  int    fehler,
         i,
         n,
         bspnummer;
  bsptyp *beispiel;
  void   *vmblock;                     /* List of dynamic allocations */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard ones           */
    return fehler;                                          /* 1 or 2 */

  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nExample:                                  ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = dgls_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("non defined example number",
                  0, __FILE__, __LINE__);
    return 3;
  }

  n = beispiel->n;
  vmblock = vminit();                           /* initialize storage */
  y0 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of available memory", 0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "absolute error bound epsabs:              ");
#endif
  scanf("%"LZS"f", &epsabs);

#ifdef INTERAKTIV
  fprintf(stderr, "relative error bound epsrel:              ");
#endif
  scanf("%"LZS"f", &epsrel);

#ifdef INTERAKTIV
  fprintf(stderr, "initial x-value  x0:                      ");
#endif
  scanf("%"LZS"f", &x0);

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "initial y-value  y0[%d] at x0:             ",
                    i);
#endif
    scanf("%"LZS"f", y0 + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "starting step size h:                     ");
#endif
  scanf("%"LZS"f", &h);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal step size hmax:                   ");
#endif
  scanf("%"LZS"f", &hmax);

#ifdef INTERAKTIV
  fprintf(stderr, "Length of interval of integration xdiff:  ");
#endif
  scanf("%"LZS"f", &xdiff);


  /* ------------ put out input data -------------------------------- */

  printf("\n"
         "Solve a first order ordinary system of DEs\n"
         "==========================================\n"
         "using the Runge-Kutta-Fehlberg method [O(h^6)]\n"
         "==============================================\n\n\n"
         "given DE system:\n"
         "----------------\n"
         "%s\n\n"
         "Input data:\n"
         "-----------\n"
         "Example  = %24d\n"
         "n        = %24d\n"
         "epsabs   = %24.15"LZP"e\n"
         "epsrel   = %24.15"LZP"e\n"
         "x0       = %24.15"LZP"e\n"
         "xdiff    = %24.15"LZP"e\n"
         "h        = %24.15"LZP"e\n"
         "hmax     = %24.15"LZP"e\n",
         (*beispiel->dgl_text)(), bspnummer, n, epsabs, epsrel,
         x0, xdiff, h, hmax);

  for (i = 0; i < n; i++)
    printf("y0[%d]    = %24.15"LZP"e\n", i, y0[i]);


  /* ------------ solve the DE system ------------------------------- */

  fehler = rk_fehl(&x0, xdiff, n, y0, beispiel->rechte_seite, &h, hmax,
                   epsabs, epsrel);

  if (fehler != 0 && fehler != 3)
  {
    fehler_melden("rk_fehl()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }


  /* --------------------- put out the results  --------------------- */

  printf("\n\n"
         "Output data:\n"
         "------------\n"
         "error code from rk_fehl():   %24d\n"
         "final step size:             %24.15"LZP"e\n"
         "final x-value  x =           %24.15"LZP"e\n"
         "\n", fehler, h, x0);

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
