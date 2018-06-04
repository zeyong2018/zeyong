#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>         /*  for  REAL, SIN, freopen, stdout, NULL, */
                           /*       fprintf, stderr, PI, TEN, ZERO,   */
                           /*       spline, spwert, COS, FABS         */
#include <kubsplin.h>      /*  for  spline                            */
#include <spliwert.h>      /*  for  spwert                            */



/* ------------------------------------------------------------------ */

static REAL fkt(REAL x)       /* Test function for differentiation    */
{
  return SIN(x);
}



/* ------------------------------------------------------------------ */

int main
        (
          int  argc,
          char *argv[]
        )

/***********************************************************************
* Test program for numerical differentiation using cubic splines.      *
* For the test function   f(x) = sin(x)  we select several equidistant *
* nodes. For this data set we compute a natural cubic spline using the *
* function spline().                                                   *
* Its derivative serves as an approximation for the first derivative   *
* of the test function.                                                *
* For comparison, we compute the derivatives of the spline at several  *
* locations and compare the results with the cosine function at these  *
* locations.                                                           *
***********************************************************************/

{
#define NMAX     20                                /* number of nodes */

  REAL x[NMAX],
       y[NMAX],
       b[NMAX],
       c[NMAX],
       d[NMAX],
       h,
       x0,
       ausg[3];
  int  n,
       fehler;


  /* --- assign input/output files to standard files ---------------- */

  if (argc >= 2)                            /* at least one argument? */
    if (freopen(argv[1], "w", stdout) == NULL)    /* open output file */
    {
      fprintf(stderr, "error opening %s!\n", argv[1]);
      return 1;
    }

  /* ---- select some equidistant nodes  ---------------------------- */

  h = PI / TEN;
  for (n = 0, x0 = ZERO; x0 <= PI; x0 += h, n++)
    x[n] = x0,
    y[n] = fkt(x0);


  /* ---- compute coefficients for natural spline through nodes ----- */

  fehler = spline(n, x, y, 2, ZERO, ZERO, 0, b, c, d);

  if (fehler)
  {
    fprintf(stderr, "*** error in spline(): error = %d\n", fehler);
    return 10 + fehler;
  }

  /* --- compare spline derivative with exact derivative ------------ */

  printf("Numerical derivation with splines\n"
         "   X         S'(X)         COS(X)         error\n"
         "-----------------------------------------------------\n"
        );

  for (x0 = ZERO; x0 <= PI; x0 += PI / (REAL)15.0)
  {
    spwert(n - 1, x0, y, b, c, d, x, ausg);
    printf("%8.4"LZP"f  %12.8"LZP"f  %12.8"LZP"f  %14.8"LZP"g\n",
           x0, ausg[0], COS(x0), FABS(ausg[0] - COS(x0)));
  }


  return 0;
}
