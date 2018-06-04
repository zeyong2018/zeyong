#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*--------------------------------------------------------------------*/
/* Test programm for fzeroin                                          */
/* The module tfunc1 must be linked.                                  */
/*--------------------------------------------------------------------*/

#include <basis.h>    /* for   REAL, MACH_EPS, printf, ZERO, TWO,     */
                      /*       NULL, ONE, HALF, FORMAT_LF, FORMAT_LE, */
                      /*       WriteHead, WriteEnd                    */
#include <tfunc1.h>   /* for   f1, f2, f3, f4                         */
#include <fzeroin.h>  /* for   zeroin                                 */



int main(int  argc,
         char *argv[]
        )

{
  REAL x1,
       x2,
       ff2,
       abserr,
       relerr;
  int  rc,
       iter,
       fno;


  /* --- assign a potential output file to the standard output    --- */

  if (argc >= 2)                           /* at least one argument?  */
    if (freopen(argv[1], "w", stdout) == NULL)        /* open output  */
    {                                                 /* file         */
      fprintf(stderr, "Error when opening %s!\n", argv[1]);
      return 1;
    }


  WriteHead("Zeroin method");

  abserr = relerr = FOUR * MACH_EPS;

  for (fno = 1; fno <= 6; fno++)
  {
    printf("Function number = %d\n", fno);
    switch (fno)
    {
      case 1:   x1 = ZERO; x2 = TWO;
                rc = zeroin(f1, &abserr, &relerr, 300, NULL,
                            x1, &x2, &ff2, &iter);
                break;
      case 2:   x1 = ZERO; x2 = TWO;
                rc = zeroin(f2, &abserr, &relerr, 300, NULL,
                            x1, &x2, &ff2, &iter);
                break;
      case 3:   x1 = -ONE; x2 = TWO;
                rc = zeroin(f3, &abserr, &relerr, 300, NULL,
                            x1, &x2, &ff2, &iter);
                break;
      case 4:   x1 = ZERO; x2 = TWO;
                rc = zeroin(f4, &abserr, &relerr, 300, NULL,
                            x1, &x2, &ff2, &iter);
                break;
      case 5:   x1 = -ONE; x2 = TWO;
                rc = zeroin(f4, &abserr, &relerr, 300, NULL,
                            x1, &x2, &ff2, &iter);
                break;
      case 6:   x1 = HALF; x2 = TWO;
                rc = zeroin(f4, &abserr, &relerr, 300, NULL,
                            x1, &x2, &ff2, &iter);
                break;
      default:  return 1;
    }

    printf("Output value:             % d\n"
           "Root:                     "FORMAT_LF
           "\nFunction value:           "FORMAT_LE
           "\nIterations:               % d\n\n",
           rc, x2, ff2, iter);
  }


  WriteEnd();
  return 0;
}
