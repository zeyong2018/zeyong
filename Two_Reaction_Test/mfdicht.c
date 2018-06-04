#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>     /*  for  PI, TRUE, FALSE, umleiten, SIN, COS,  */
                       /*       REAL, printf, scanf, NULL, fprintf,   */
                       /*       stderr, FABS, REAL, TWO,              */
                       /*       fehler_melden, LZP                    */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vminit, VVEKTOR  */
#include <fft.h>       /*  for  complex, fdicht                       */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
*         This program tests the function fdicht()   with functional   *
*         values of the sum of the complex basic functions             *
*                                                                      *
*              2 * EXP(I*(-3)*OMEGA*X) + I * EXP(I*7*OMEGA*X)          *
*                                                                      *
*         ( I : imaginary unit, OMEGA : 2*PI/P, P : period )           *
*         for different sets of equidistant nodes from the interval of *
*         periodicity  [0, P). It then compares the function values at *
*         the shifted nodes with the exact ones, which should differ   *
*         only the up to errors on the level of the machine constant.  *
***********************************************************************/

{
#define ANZAHLBEISPIELE     5
#define maxM             4096

  static int  M[ANZAHLBEISPIELE] = { 16, 81, 263, 1101, maxM };
  static REAL p[ANZAHLBEISPIELE] = { (REAL)6.283185308, (REAL)7.0,
                                     (REAL)3.14159, TWO, TEN
                                   };

  complex     *F;       /* [0..maxM-1] vector; table of values of test*/
                        /* function (1)                               */
  complex     fj;       /* "exact" value for F[j] (from fdicht())     */
  complex     diff;     /* distance between `exact' and computed value*/
  REAL        diffmax;  /* accuracy of result from fdicht()           */
  REAL        fakt;     /* step size                                  */
  REAL        xj;       /* current node                               */
  REAL        weg;
  REAL        xjweg;
  int         fehler;   /* error code of  fdicht()                    */
  int         i;        /* variable                                   */
  int         j;        /* variable                                   */
  void        *vmblock; /* List of dynamically allocated vectors      */


  /* --- assign output file to standard output                     -- */

  if (argc >= 2)                           /* at least one argument?  */
    if (freopen(argv[1], "w", stdout) == NULL)        /* open output  */
    {                                                 /* file      */
      fprintf(stderr, "Error in opening %s!\n", argv[1]);
      return 1;
    }


  vmblock = vminit();                 /* initialize storage           */
  F = (complex *)vmalloc(vmblock, VVEKTOR, maxM, sizeof(*F));
  if (! vmcomplete(vmblock))   /* Storage demand not                  */
  {                            /*          satisfied?                 */
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }


  for (i = 0; i < ANZAHLBEISPIELE; i++)
  {
    /* find real and imaginary parts of functional values */
    fakt = TWO * PI / (REAL)M[i];
    for (j = 0; j < M[i]; j++)
    {
      xj     = (REAL)j * fakt;
      F[j].x =  TWO * COS(THREE * xj) - SIN((REAL)7.0 * xj);
      F[j].y = -TWO * SIN(THREE * xj) + COS((REAL)7.0 * xj);
      F[j].y = ZERO;
    }
    weg = HALF * p[i] / (REAL)M[i];

    fehler = fdicht(M[i], F, p[i], weg);

    switch (fehler)
    {
      case 0:
        break;
      case 1:
        fehler_melden("fdicht(): M < 1", 10 + fehler,
                      __FILE__, __LINE__);
        break;
      case 2:
        fehler_melden("fdicht(): M too large",
                      10 + fehler, __FILE__, __LINE__);
        break;
      case 3:
        fehler_melden("fdicht(): memory exhausted",
                      10 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("fdicht(): unknown eror",
                      10 + fehler, __FILE__, __LINE__);
    }

    /* Check computed functioanl values                    */
    for (diffmax = ZERO, j = 0; j < M[i]; j++)
    {
      xjweg   = (REAL)j * fakt + TWO * PI / p[i] * weg;
      fj.x    =  TWO * COS(THREE * xjweg) - SIN((REAL)7.0 * xjweg);
      fj.y    = -TWO * SIN(THREE * xjweg) + COS((REAL)7.0 * xjweg);
      fj.y    = ZERO;
      diff.x  = sqr(F[j].x - fj.x);
      diff.y  = sqr(F[j].y - fj.y);
      diffmax = max(diffmax, SQRT(diff.x + diff.y));
    }
    printf("number = %4d  period = %10.5"LZP"f"
           "  error = %12.5"LZP"e\n", M[i], p[i], diffmax);
  }


  return 0;
}
