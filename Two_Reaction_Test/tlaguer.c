#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>    /* for   REAL, umleiten, scanf, printf, FOUR,   */
                      /*       FORMAT_2010LF, MACH_EPS, WriteHead,    */
                      /*       LogError, ReadVec, WriteEnd            */
#include <vmblock.h>  /* for   vminit, vmalloc, VEKTOR, VVEKTOR,      */
                      /*       vmcomplete                             */
#include <flaguer.h>  /* for   laguerre                               */



/* ------------------------------------------------------------------ */

int main
        (
         int  argc,
         char *argv[]
        )

/***********************************************************************
* This is a test programm for Laguerre's method to find all real roots *
* of a real polynomial.                                                *
*                                                                      *
* Way of operation:                                                    *
* =================                                                    *
* The program reads the input from file stdin and writes the output    *
* onto file stdout. If the first command line parameter is not void,   *
* it is interpreted as the input file and is associated with stdin.    *
* Analogously, the second command line parameter is associated with    *
* the output file stdout.                                              *
*                                                                      *
* After reading the input files, their content is printed as a         *
* safeguard. Then the Laguerre method is applied and its output is     *
* stored.                                                              *
*                                                                      *
* Form of the input files:                                             *
* ========================                                             *
* n       degree of the polynomial p                                   *
* a[0] \  coefficients (in ascending order) of the polynomial p:       *
* ...   >     p(x)  =  a[0] + a[1] * x + ... + a[n] * x^n              *
* a[n] /                                                               *
***********************************************************************/

{
  REAL *a,
       *x;
  int  *iter,
       i,
       n,
       nulanz,
       fehler;
  void *vmblock;


  if ((fehler = umleiten(argc, argv))   /* if necessary assign the    */
      != 0)                             /* input and/or output files  */
    return fehler;  /* 1 or 2 */        /* to the standard input or   */
                                        /* output files               */


  WriteHead("Method of Laguerre");

  if (scanf("%d", &n) <= 0)
  {
    LogError("wrong input", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError("degree of the polynomial must exceed zero.",
             0,  __FILE__, __LINE__);
    return 1;
  }

  vmblock = vminit();
  a    = (REAL *)vmalloc(vmblock, VEKTOR,  n + 1, 0);
  x    = (REAL *)vmalloc(vmblock, VEKTOR,  n,     0);
  iter = (int  *)vmalloc(vmblock, VVEKTOR, n,     sizeof(*iter));
  if (! vmcomplete(vmblock))
  {
    LogError("Memory overflow", 0, __FILE__, __LINE__);
    return 1;
  }

  fehler = ReadVec(n + 1, a);
  if (fehler)
  {
    LogError("wrong input", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf("coefficients (in ascending order):\n");

  for (i = 0; i <= n; i++)
  {
    printf("\n a[%2d] = ", i);
    printf(FORMAT_2010LF, a[i]);
  }

  fehler = laguerre(n, a, FOUR * MACH_EPS, FOUR * MACH_EPS, 500,
                    x, iter, &nulanz);
  if (! fehler)
  {
    printf("\n\n");
    printf("Number  root                      step count   "
           "function value\n");
    for (i = 0; i < nulanz; i++)
      printf(" \n %5d  "FORMAT_2016LE"       %5d  "FORMAT_LE,
             i, x[i], iter[i], horner(n, a, x[i]));
    printf("\n");
  }
  else
  {
    LogError("laguerre()", fehler,  __FILE__, __LINE__);
    return 1;
  }


  WriteEnd();

  return 0;
}
