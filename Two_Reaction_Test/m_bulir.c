#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, fprintf, stderr, printf,  */
                         /*       scanf, NULL, REAL, LZS, LZP,        */
                         /*       fehler_melden                       */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vmfree,        */
                         /*       vminit, VEKTOR                      */
#include <bulirsch.h>    /*  for  bul_stoe                            */
#include <t_dgls.h>      /*  for  bsptyp, dgls_waehlen                */
#ifdef __TURBOC__
#ifdef __MSDOS__
#include <alloc.h>       /*  for  coreleft                            */
#endif
#endif



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
* Test program for the function bul_stoe() from the module bulirsch    *
* for solving a first order ordinary system of differential equations  *
* using the extrapolation method of  Bulirsch-Stoer-Gragg.             *
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
* the extrapolations method is used as often as needed to reach the    *
* desired x-value, except for a fatal error and the computed results   *
* are put out.                                                         *
*                                                                      *
* To solve a differential equation, please proceed as in example       *
* t_dgls.c.                                                            *
*                                                                      *
* Construction of input data files:                                    *
* =================================                                    *
* bspnummer  Number of DE system in  t_dgls.c                          *
* x0         initial x-value                                           *
* xend       final desired x-value                                     *
* y0[0]   \  initial y-value at x0                                     *
* ...      >                                                           *
* y0[n-1] /                                                            *
* epsabs     desired absolute error bound                              *
* epsrel     desired relative error bound                              *
* fmax       upper bound for the number of allowed function            *
*            evaluations of the right hand side of the DE system       *
* h          starting step size                                        *
* hmax       maximal step size                                         *
*                                                                      *
* The number n of differential equations follows from the number of    *
* the DE system chosen; it is stored in conjunction with the right     *
* hand side function in t_dgls.c.                                      *
***********************************************************************/

{
  REAL   x0,
         xend,
         *y0,
         *yex,
         epsabs,
         epsrel,
         h,
         hmax;
  long   fmax,           /* upper bound for the number of allowed     */
                         /* function evaluations of the right hand    */
                         /* side of the DE system                     */
         aufrufe;        /* number of calls of dgl()                  */
  int    bspnummer,
         n,
         fehler,         /* error code from bul_stoe()                */
         i;              /* loop counter                              */
  bsptyp *beispiel;      /* pointer to the structure that describes   */
                         /* the DE system at hand                     */
  void   *vmblock;       /* List of dynamic allocations               */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard ones           */
    return fehler;  /* 1 or 2 */

  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "Example:                                     ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = dgls_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("non existing example number",
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
  fprintf(stderr, "initial x-value  x0:    ");
#endif
  scanf("%"LZS"f", &x0);

#ifdef INTERAKTIV
  fprintf(stderr, "desired final x-value xend: ");
#endif
  scanf("%"LZS"f", &xend);

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "Function value y0[%d] at x0:            ",
                    i);
#endif
    scanf("%"LZS"f", y0 + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "absolute error bound epsabs:               ");
#endif
  scanf("%"LZS"f", &epsabs);

#ifdef INTERAKTIV
  fprintf(stderr, "relative error bound epsrel:               ");
#endif
  scanf("%"LZS"f", &epsrel);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of calls of right side:     ");
#endif
  scanf("%ld", &fmax);

#ifdef INTERAKTIV
  fprintf(stderr, "starting step size h:                  ");
#endif
  scanf("%"LZS"f", &h);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal step size hmax:       ");
#endif
  scanf("%"LZS"f", &hmax);


  /* ------------ put out input ------------------------------------- */

  printf("\n"
         "Solve a first order ordinary system of differential "
         "equations\n"
         "===================================================="
         "=========\n"
         "using the extrapolation method of Bulirsch-Stoer-Gragg\n"
         "======================================================\n"
         "\n\n"
         "DE system:\n"
         "----------\n"
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


  /* ------------ solve the system of DEs -------------------------- */
#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "\nfree before:  %u\n", coreleft());
#endif
#endif
#endif
  fehler = TRUE;  /* report that no old values are used              */
                  /* (they do not exist yet)                         */
  do
  {
    fehler = bul_stoe(&x0, xend, n, beispiel->rechte_seite, y0, epsabs,
                      epsrel, &h, hmax, fehler, fmax, &aufrufe);
  }
  while (fehler == 1);   /* repeat call with old values as xend was  */
                         /* not reached due to excessive calls of    */
                         /*  dgl()                                   */
#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "free after: %u\n", coreleft());
#endif
#endif
#endif

  if (fehler != 0)
  {
    fehler_melden("bul_stoe()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }

  /*  The following call of "bul_stoe()" is used to free the allocated*/
  /*  storage now by calling with an improper value of n = 0.         */
  (void)bul_stoe(&x0, xend, 0, beispiel->rechte_seite, y0, epsabs,
                 epsrel, &h, hmax, fehler, fmax, &aufrufe);
#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "free after: %u\n", coreleft());
#endif
#endif
#endif


  /* -------------------- put out results -------------------------- */

  printf("\n\n"
         "Output data:\n"
         "------------\n"
         "error code of bul_stoe():                 %24d\n"
         "final step size:                          %24.15"LZP"e\n"
         "calls of dgl():                           %24ld\n"
         "final x-value x =                         %24.15"LZP"e\n\n",
         fehler, h, aufrufe, x0);

  for (i = 0; i < n; i++)
    printf("approximate solution y%d(x) = %24.15"LZP"e\n",
           i + 1, y0[i]);

  if (beispiel->exakte_loesung != NULL)       /* exact solution       */
  {                                           /* available?           */
    (*beispiel->exakte_loesung)(x0, yex);
    printf("\n");
    for (i = 0; i < n; i++)
      printf("'exact' solution  y%d(x)    = %24.15"LZP"e\n",
             i + 1, yex[i]);
    printf("\ndifference approx. solution - 'exact' solution:\n");
    for (i = 0; i < n; i++)
      printf("%24.15"LZP"g\n", y0[i] - yex[i]);
  }


  return 0;
}
