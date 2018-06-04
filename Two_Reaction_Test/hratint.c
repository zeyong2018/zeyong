#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test program for rational interpolation              Egg, 11.17.1991 *
***********************************************************************/

#include <basis.h>
#include <stdio.h>
#include <ratint.h>

#define TestFunktion  2
#define INTERVMIN (-ONE)
#define INTERVMAX ONE
#define Nmax      (20)
#define MAXI      (Nmax+1)


#if TestFunktion == 1

   REAL f (REAL x) { return ( (x*x + (REAL)2.2) * x - (REAL)45.0); }
   char  *fname = "(x*x + 2.2) * x - 45.0";

#elif TestFunktion == 2

   REAL f (REAL x) { return ONE / ((x-THREE) * (x-THREE)); }
   char  *fname = "1 / ((x-3) * (x-3))";

#else /* TestFunktion == 3 */

   REAL f (REAL x) { return log (x+THREE); }
   char  *fname = "log (x+3)";

#endif

int main (int argc, char *argv[])
{
  int    n, i, error, numgrad, inf_c[MAXI];
  REAL   x[MAXI], y[MAXI];
  REAL   eps, x0;

  /* --- assign input file to standard input file ------------------- */
  if (argc >= 2)                               /* at least one entry? */
    if (freopen(argv[1], "r", stdin) == NULL)      /* open input file */
    {
      fprintf(stderr, "error opening %s!\n", argv[1]);
      return 2;
    }

  printf ("Test program:             Rational interpolation\n");
  printf ("Test function:            %s\n", fname);
  printf ("Interval:                 %"LZP"g .. %"LZP"g\n",
          INTERVMIN, INTERVMAX);

  eps = (REAL)1.0e-10;

  printf ("number of nodes : ");
  scanf  ("%d", &n);
  n--;                                               /* use n - 1 !  */
  for (i = 0; i <= n; i++)                           /* compute nodes*/
  {
    x[i] =  INTERVMIN + (INTERVMAX-INTERVMIN) * i / n;
    y[i] =  f (x[i]);
  }

  printf ("degree of numerator polynomial:  ");
  scanf  ("%d", &numgrad);

  error = ratint (n, numgrad, x, y, inf_c, eps);

  if (error == 0)
  {
    printf ("%-15s  %-15s  %-15s  %-15s\n",
            "place of evaluation","computed value",
            "Original value",     "Difference");
    for (i = 0; i <= Nmax; i++)
    {
      REAL rv, fv;
      x0 = INTERVMIN + (INTERVMAX-INTERVMIN) * i / Nmax;
      rv = ratval (n, x0, x, y, inf_c);
      fv = f (x0);
      printf ("%15.5"LZP"f  %15.5"LZP"f  %15.5"LZP"f  %15.10"LZP"f\n",
              x0, rv, fv, rv-fv);
    }
  }
  else printf ("error number = %d\n",error);

  return 0;
}

/***********************************************************************
* The accuracy is easily checked: the original and computed values     *
* must be nearly identical, or their difference must be close to zero. *
***********************************************************************/
