#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  REAL, umleiten, fprintf, stderr,    */
                         /*       scanf, NULL, fehler_melden, readln, */
                         /*       LZS, printf, LZP                    */
#include <vmblock.h>     /*  for  vminit, vmalloc, VEKTOR, vmcomplete */
#include <brown.h>       /*  for  brown                               */
#include <nlglstst.h>    /*  for  bsptyp1, nlgls_waehlen              */



/***********************************************************************
* This program tests the function brown() which is designed to solve   *
* a nonlinear system of equations                                      *
*             f0  (x0,x1,...,xn-1)  =  0                               *
*             f1  (x0,x1,...,xn-1)  =  0                               *
*             ...                                                      *
*             fn-1(x0,x1,...,xn-1)  =  0                               *
* using Brown's method                                                 *
*                                                                      *
* Scope of program :                                                   *
* ==================                                                   *
* The program reads the input from stdin and writes output onto stdout.*
* If the first command line parameter is not void, it is taken  as     *
* a file name and associated with stdin, similarly for the second      *
* one and stdout.                                                      *
* Calls for input are written onto stderr which also receives error    *
* messages. If before compiling the macro INTERAKTIV has not been set  *
* up, the input calls are not recorded.                                *
*                                                                      *
* After reading, the input is printed out for control purposes. Then   *
* Brown's method is performed and, barring premature errors, the end   *
* results are printed.                                                 *
*                                                                      *
* If the user wants to test his own systems, he should follow the      *
* example in nlglstst.c.                                               *
*                                                                      *
* Form of the input file :                                             *
* ========================                                             *
* bspnummer  Number of the example from  nlglstst.c                    *
* eps        desired accuracy                                          *
* prot       Protocol is kept if prot = 1, not kept for prot = 0       *
* maxit      maximal number of steps                                   *
* x0[0]  \                                                             *
* x0[1]   \  Starting vector for iteration                             *
* ...     /                                                            *
* x0[n-1]/                                                             *
*                                                                      *
* The order of the system is prescribed via the number of the test     *
* example.                                                             *
***********************************************************************/

int main(int  argc,
         char *argv[]
        )

{
  int     n,             /* size of system                            */
          prot,          /* Protocol flag                             */
          maxit,         /* maximal number of iterations              */
          fehler,        /* error code of brown()                     */
          bspnummer,     /* Number of chosen example                  */
          itanz,         /* number of iterations performed            */
          i;             /* Loop variable                             */
  REAL    eps,           /* desired accuracy                          */
          *x0,           /* [0..n-1] starting vector                  */
          *x1;           /* [0..n-1] approximate solution             */
  bsptyp1 *beispiel;     /* pointer to the system of current example  */
  void    *vmblock;      /* List of dynamically allocated vectors     */


  if ((fehler = umleiten(argc, argv))   /* assign in/ouput files to   */
      != 0)                             /* standard in/output files   */
    return fehler;  /* 1 or 2 */

  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nExample:                               ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = nlgls_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("not one of the listed examples",
                  0, __FILE__, __LINE__);
    return 3;
  }
  readln();

  n = beispiel->n;
  vmblock = vminit();                 /* initialize files             */
  x0 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  x1 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))   /* allocation successful ?             */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "desired accuracy:                         ");
#endif
  scanf("%"LZS"f", &eps);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr, "with (1) or without (0) protocol:         ");
#endif
  scanf("%d", &prot);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of iterations:             ");
#endif
  scanf("%d", &maxit);
  readln();
  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "Component %3d of the starting vector x0:  ", i);
#endif
    scanf("%"LZS"f", x0 + i);
  }


  /* ------------ print input for checking purposes ----------------- */

  printf("\n\n"
         "Brown's method for nonlinear systems of equations\n"
         "=================================================\n"
         "\n"
         "system to be solved\n"
         "%s\n"
         "Starting vector\n",
         beispiel->fkt_text()
        );

  for (i = 0; i < n; i++)
    printf(" %9.3"LZP"e", x0[i]);
  printf("\n\n"
         "error bound                  = %"LZP"e\n"
         "maximal number of iterations = %d\n",
          eps, maxit
         );
  if (prot)
    printf("Intermediate results are saved.\n");
  else
    printf("Intermediate results are not kept.\n");


  /* ------------ solve nonlinear system ---------------------------- */

  fehler = brown(beispiel->fkt, n, x0, eps, prot, maxit, x1, &itanz);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("brown(): too many steps",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("brown(): linearized system singular",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("brown(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 4:
      fehler_melden("brown(): wrong input parameter: "
                    "fkt = NULL or n < 1 or maxit < 1",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 5:
      fehler_melden("brown(): error calling fkt()",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("brown(): other error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* --------------------- print solution --------------------------- */

  printf("\nsolution vector\n");
  for (i = 0; i < n; i++)
    printf("%16"LZP"e ", x1[i]);
  printf("\n\n%d number of iterations\n", itanz);


  return 0;
}
