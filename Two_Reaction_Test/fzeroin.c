#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE fzeroin.c ----------------------- */

/***********************************************************************
*                                                                      *
* Compute a zero of a continuous real valued function with the         *
* ------------------------------------------------------------         *
* Zeroin method.                                                       *
* --------------                                                       *
*                                                                      *
* exported funktion:                                                   *
*   - zeroin():  Zeroin- method for real functions                     *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 mit 80387                          *
* Author:               Siegmar Neuner                                 *
* Modifications:        Juergen Dietel, Rechenzentrum, RWTH Aachen     *
* Source:               FORTRAN source code                            *
* Date:                 11. 27. 1992                                   *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  REAL, REALFCT, FOUR, MACH_EPS, ZERO, */
                        /*       FABS, HALF, NULL, TWO, ONE, THREE,   */
                        /*       FILE, fprintf, REALFCT, sign         */
#include <fzeroin.h>    /*  for  zeroin                               */



/* ------------------------------------------------------------------ */
/*.BA*/

int zeroin                       /* Find roots with the Zeroin method */
/*.IX{zeroin}*/
          (
           REALFCT fkt,         /* Function ..........................*/
           REAL    *abserr,     /* absolute error bound ..............*/
           REAL    *relerr,     /* relative error bound ..............*/
           int     fmax,        /* maximal number of calls for fkt() .*/
           char    *protnam,    /* Name of the log file ..............*/
           REAL    a,           /* [a,b] = inclusion interval ........*/
           REAL    *b,          /* right endpoint or zero ............*/
           REAL    *fb,         /* Function value at the root b ......*/
           int     *fanz        /* number of actual function calls ...*/
          )                     /* error code ........................*/

/***********************************************************************
* Given a real valued function f on an interval [a,b] with             *
* f(a) * f(b) < 0, we compute a root of f in [a,b].                    *
* The Zeroin method combines bisection and secant methods with inverse *
* quadratic interpolation.                                             *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* fkt      Function f, whose root we want to compute                   *
* abserr\  Error bounds with   abserr >= 0  and  relerr >= 0. Their    *
* relerr/  sum must exceed zero. For break-off we use the test         *
*              |xm|  <=  0.5 * (abserr + |b| * relerr),                *
*          where xm denotes half the length of the final inclusion     *
*          interval. For relerr = 0 we test only for absolute error,   *
*          while for  abserr = 0, we test the relative error only.     *
*          abserr and relerr are used as put in only when both exceed  *
*          four times the machine constant. Otherwise we set them to   *
*          this value.                                                 *
* fmax     upper limit of calls of fkt()                               *
* protnam  Name of a file used for intermediate results. If the pointer*
*          is set to zero, we do not use this file.                    *
* a,b      end points of the interval, that includes a root            *
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* abserr\  actual error bounds used                                    *
* relerr/                                                              *
* b        approximate root                                            *
* fb       value of f at the root b                                    *
* fanz     number of calls of  fkt()                                   *
*                                                                      *
* Function value:                                                      *
* ==============                                                       *
* = -2: abserr or relerr is negative, or both are zero, or fmax < 1.   *
* = -1: The necessary assumption  f(a) * f(b) < 0  is violated.        *
* =  0: Desired accuracy has been reache :                             *
*           |xm|  <=  0.5 * (abserr + |b| * relerr).                   *
* =  1: b on output is a root with  fkt(b) = 0.                        *
* =  2: Either a or b is a root at input.                              *
* =  3: After fmax calls of fkt() the algorithm has not found a root   *
* =  4: Error when opening intermediate result file                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, REALFCT, FOUR, MACH_EPS, ZERO, FABS, HALF, NULL, bi, TWO, ONE, *
* THREE, FILE, fprintf                                                 *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL fa,                /* Function value fkt(a)                    */
       fc,                /* Function value fkt(c)                    */
       eps,               /* minimal value for error bounds abserr    */
                          /* and relerr                               */
       tol1,              /* auxiliary variable for mixed error bound */
       xm,                /* half of the current interval length      */
       c = ZERO,          /* value inside [a,b]                       */
       d = ZERO,          /* Distance to the nearest approximate root */
       e = ZERO,          /* previous value of d                      */
       p,
       q,
       r,
       s,
       tmp;               /* auxiliary variable to check inclusion    */
                          /*   f(a) * f(b) < 0                        */
  int  fehler;            /* error code of this function              */
  FILE *prodat;           /* intermediate data file                   */


  /* ----------------- initialize variables ------------------------- */

  eps = FOUR * MACH_EPS;

  fa    = fkt(a);                   /* evaluate fkt() at the end      */
  *fb   = fkt(*b);                  /* points a and b                 */
  *fanz = 2;

  /* ---------- check   f(a) * f(b) < 0  ---------------------------- */

  if ((tmp = fa * *fb) > ZERO)
    return -1;
  else if (tmp == ZERO)
    return 2;

  /* ---------- check usability of given error bounds -------------- */

  if (*abserr < ZERO || *relerr < ZERO || *abserr + *relerr <= ZERO ||
      fmax < 1)
    return -2;
  if (*relerr == ZERO && *abserr < eps)
    *abserr = eps;
  else if (*abserr == ZERO && *relerr < eps)
    *relerr = eps;
  else
  {
    if (*abserr < eps)
      *abserr = eps;
    if (*relerr < eps)
      *relerr = eps;
  }


  prodat = NULL;
  if (protnam != NULL &&                      /* intermediate results?*/
      (prodat = fopen(protnam, "w")) == NULL) /* Error opening file ? */
    return 4;


  for (fc = *fb; TRUE; )                         /* start iteration  */
  {

    if (*fb * (fc / FABS(fc)) > ZERO)     /* no inclusion of a root   */
      c  = a,                             /* between b and c ?        */
      fc = fa,                            /* alter c so that b and c  */
      e  = d   = *b - a;                  /* include the root of f    */


    if (FABS(fc) < FABS(*fb))        /* If fc has the smaller modulus */
      a   = *b,                      /* interchange interval end      */
      *b  = c,                       /* points.                       */
      c   = a,
      fa  = *fb,
      *fb = fc,
      fc  = fa;

    if (prodat != NULL)                /* Want intermediate results? */
    {
      fprintf(prodat, "a = %20.14"LZP"f  b = %20.14"LZP"f  "
                      "c = %20.14"LZP"f\n", a, *b, c);
      fprintf(prodat, "fa= %20.14"LZP"f  fb= %20.14"LZP"f  "
                      "fc= %20.14"LZP"f\n", fa, *fb, fc);
    }

    tol1 = HALF * (*abserr + *relerr * FABS(*b));
    xm   = HALF * (c - *b);

    if (FABS(xm) <= tol1)               /* reached desired accuracy ? */
    {
      fehler = 0;
      break;
    }
    if (*fb == ZERO)                /* Is the best approximate root b */
    {                               /* a precise root of f ?          */
      fehler = 1;
      break;
    }

    r = ZERO;

    if (FABS(e) < tol1 || FABS(fa) <= FABS(*fb))
    {
      e = d = xm;
      if (prodat != NULL)
        fprintf(prodat, "Bisection\n");
    }

    else
    {

      if (a != c)            /* if a is not equal to c :              */
        q =  fa / fc,        /* With a, b and c we have 3 points for  */
        r = *fb / fc,        /* an inverse quadratic interpolation    */
        s = *fb / fa,
        p = s * (TWO * xm * q * (q - r) - (*b - a) * (r - ONE)),
        q = (q - ONE) * (r - ONE) * (s - ONE);

      else                   /* If a equals  c :                      */
        s = *fb / fa,        /* Use the secant method or linear       */
        p = TWO * xm * s,    /* interpolation                         */
        q = ONE - s;

      if (p > ZERO)          /* alter the sign of  p/q for the        */
        q = -q;              /* subsequent division                   */
      else
        p = FABS(p);

      if ((TWO * p  >=  THREE * xm * q - FABS(tol1 * q)) ||
          (p        >=  FABS(HALF * e * q))
         )
      {
        e = d = xm;
        if (prodat != NULL)
          fprintf(prodat, "Bisection\n");
      }

      else
      {
        e = d;        /* Compute the quotient p/q for both iterations */
        d = p / q;    /* which will be used to modify b               */

        if (prodat != NULL)            /* Want intermediate data ?    */
          if (r == ZERO)
            fprintf(prodat, "Secant method\n");
          else
            fprintf(prodat, "Inverse quadratic interpolation\n");
      }
    }

    a  = *b;         /* store the best approximation b and its        */
    fa = *fb;        /* function value fb in a and fa.                */

    if (FABS(d) > tol1)                         /* d large enough?    */
    {
      *b += d;                                  /* add d to b         */

      if (prodat != NULL)              /* Want intermediate data?     */
        fprintf(prodat, "Difference d from new b:  "
                        "d = %20.14"LZP"f", d);
    }

    else                                        /* d too small?       */
    {                                           /* add tol1 to b      */
      *b += sign(tol1, xm);

      if (prodat != NULL)              /* Want intermediate data ?    */
        if (xm < ZERO)
          fprintf(prodat, "Subtract error bound: "
                          "d = %20.14"LZP"f\n", -tol1);
        else
          fprintf(prodat, "Add error bound:     "
                          "d = %20.14"LZP"f\n", tol1);
    }

    *fb = fkt(*b);                     /* compute function value at b */
    ++*fanz;                           /*      up evaluation counter  */

    if (prodat != NULL)                /* Want intermediate data?     */
      fprintf(prodat, "b = %20.14"LZP"f  fb= %20.14"LZP"f\n"
                      "Number of functional evaluations = %4d",
                      *b, *fb, *fanz);

    if (*fanz > fmax)               /* too many function evaluations? */
    {
      fehler = 3;
      break;
    }
  }                                                  /* end iteration */


  if (protnam != NULL)          /* If intermediate data file was in   */
    fclose(prodat);             /* use, close it now.                 */

  return fehler;
}

/* ------------------------- END fzeroin.c -------------------------- */
