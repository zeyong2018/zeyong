#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/
/*.KA{C 2}{Nonlinear Equations in One Variable}
          {Nonlinear Equations in One Variable}*/
/*.FE{C 2.5.1}
     {The Newton Method for Simple Roots}
     {The Newton Method for Simple Roots}*/

/*.BE*/
/* ------------------------- MODULE fnewton.c ----------------------- */

#include <basis.h>
#include <u_proto.h>

#define ITERMAX 300                    /* Maximal number of iterations*/

#define ABSERR ZERO                    /* acceptable absolute errorr  */

                                       /* acceptable relative error   */
#define RELERR (REAL)((REAL)128.0 * MACH_EPS)

                                       /* Maximal functional error    */
#define FCTERR (REAL)((REAL)4.0 * MACH_EPS)


/*.BA*/

int newton              /* Newton method in one dimension    .........*/
/*.IX{newton}*/
           (
            REALFCT   fct,        /* Function ........................*/
            REALFCT   fderv,      /* 1st derivative ..................*/
            REAL *    x,          /* Starting value / solution .......*/
            REAL *    fval,       /* Functional value at the root ... */
            int  *    iter        /* Number of iterations ............*/
           )
/*====================================================================*
 *                                                                    *
 *  The function newton uses Newton's iterative method to solve       *
 *  the equation fct(x) = 0.                                          *
 *  The function fct and its first derivative must be available as    *
 *  given input parameters.                                           *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   ============                                                     *
 *      To find roots of a continously differentiable function fct.   *
 *      The user must supply a sufficiently good approximation to     *
 *      the root.                                                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      fct      REAL fct (REAL);                                     *
 *               Function, whose root is to be found.                 *
 *               fct has the form:                                    *
 *               REAL fct (REAL x)                                    *
 *               { REAL f;                                            *
 *                 f = ...;                                           *
 *                 return (f);                                        *
 *               }                                                    *
 *      fderv    REAL fderv (REAL);                                   *
 *               1st derivative of fct; the function fderv has the    *
 *               same form as fct.                                    *
 *      x        REAL   *x;                                           *
 *               Starting value for the iteration.                    *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      x        REAL   *x;                                           *
 *               Computed approximate root of fct.                    *
 *      fval     REAL   *fval;                                        *
 *               Functional value of fct at x; this must be near zero.*
 *      iter     int *iter;                                           *
 *               Number of iteration steps.                           *
 *                                                                    *
 *   Return values:                                                   *
 *   =============                                                    *
 *      = 0      root x has been found with abs(fct(x)) < FCTERR      *
 *      = 1      stop with abs(xnew-xold) < ABSERR + xnew * RELERR    *
 *      = 2      maximal number of iterations reached                 *
 *      = 3      unadmissable input fumctions fct or fderv            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *                                                                    *
 *      REAL   fct ():  Function, whose zero is to be computed        *
 *      REAL   fderv(): 1st derivative of fct                         *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used: ABSERR, RELERR, MACH_EPS, FCTERR, ITERMAX,       *
 *   ==============  EPSROOT                                          *
 *                                                                    *
 *====================================================================*/
{
  REAL fs, diff, y = ONE;
  int rc;

  *iter = 0;

  if (fct == NULL || fderv == NULL) return (3);

  rc = 2;
  while (*iter < ITERMAX)                       /* Newton iteration   */
  {
    y = (*fct) (*x);                           /* Function value at x */

    if (ABS (y) < FCTERR)     /* Function value < FCTERR > acceptable */
    {
      rc = 0;
      break;
    }

    fs = (*fderv) (*x);                        /* 1st derivative at x */

    if (ABS (fs) < EPSROOT) fs = EPSROOT;

    (*iter)++;

    diff = y / fs;
    *x -= diff;                            /* xnew = xold - fval / fs */

    if (ABS (diff) < ABS (*x) * RELERR + ABSERR)
    {
      rc = 1;
      break;
    }
  }

  *fval = y;
  return (rc);

}


/*.BA*/



/*.FE{}
     {The Newton Method for Polynomials}
     {The Newton Method for Polynomials}*/

int newpoly              /* Newton method for polynomials  ...........*/
/*.IX{newpoly}*/
            (
             int       n,         /* degree of polynomial ............*/
             REAL      coeff[],   /* vector of coefficients ..........*/
             REAL *    x,         /* Starting value / solution .......*/
             REAL *    fval,      /* Functional value at x ...........*/
             int  *    iter       /* Number of iterations ............*/
            )
/*====================================================================*
 *                                                                    *
 *  The function newpoly performs Newton iteration to solve the       *
 *  polynomial equation fct(x) = 0, if fct is a given polynomial.     *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   ============                                                     *
 *    To find a root of the nth degree polynom:                       *
 *                                                n-1              n  *
 *    coeff[0] + coeff[1] * x +...+ coeff[n-1] * x   + coeff[n] * x   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n; (n > 0)                                       *
 *               degree of the polynom.                               *
 *      coeff    REAL   *coeff;                                       *
 *               vector of polynomial coefficients                    *
 *               coeff[0],..,coeff[n];                                *
 *      x        REAL   *x;                                           *
 *               Starting value for the iteration.                    *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      x        REAL   *x;                                           *
 *               Computed approximate root of fct                     *
 *      fval     REAL   *fval;                                        *
 *               Functional value at the computed root, which ought   *
 *               to be nearly zero.                                   *
 *      iter     int *iter;                                           *
 *               Number of iterations performed.                      *
 *                                                                    *
 *   Return value:                                                    *
 *   ============                                                     *
 *      = 0      root found with abs(fct) < FCTERR                    *
 *      = 1      stop with abs(xnew-xold) < ABSERR + xnew * RELERR    *
 *      = 2      maximal number of iterations reached                 *
 *      = 3      disallowed input parameter                           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Used functions:                                                  *
 *   ==============                                                   *
 *                                                                    *
 *      REAL   polval (): Horner scheme to evaluate a polynomial and  *
 *                        its derivative.                             *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Used constants: ABSERR, RELERR, MACH_EPS, FCTERR, ITERMAX,       *
 *   ==============  EPSROOT                                          *
 *                                                                    *
 *====================================================================*/
{
  REAL fs, diff, y;

  if (n < 1 || coeff[n] == ZERO) return (3); /* wrong parameter       */
  if (coeff == NULL) return (3);

  *iter = 0;                      /* Initialize iterations counter    */

  while (*iter < ITERMAX)
  {
    polval (n, coeff, *x, &y, &fs);       /* Compute value of the     */
                                          /* polynomial and its       */
                                          /* derivative               */
    if (ABS (y) <= FCTERR)                /* Accurate enough ?        */
    {
      *fval = y;
      return (0);
    }

    if (ABS (fs) < EPSROOT)       /* To avoid overflow, bound the     */
      fs = EPSROOT;               /* value of the derivative to       */
                                  /* EPSROOT                          */

    (*iter)++;                    /* Up iterations counter            */

    diff = y / fs;                /* new step size                    */
    *x -= diff;                   /* xnew =                           */
                      /*   = xold - function value/1st derivative     */

    if (ABS (diff) <= ABS (*x) * RELERR + ABSERR)
    {                             /* step size sufficiently small ?   */
      *fval = y;
      return (1);
    }
  }

  *fval = y;
  return (2);
}


int polval                /* Newton method for polynomials ...........*/
/*.IX{polval}*/
           (
            int     n,            /* degree of the polynomial ........*/
            REAL    coeff[],      /* Vector of coefficients ..........*/
            REAL    x,            /* Value for evaluation ............*/
            REAL *  val,          /* value of the polynomial at x ....*/
            REAL *  dval          /* 1st derivative at x .............*/
           )
/*====================================================================*
 *                                                                    *
 *  polval evaluates a polynomial, that is given by its coefficient   *
 *  vector coeff in ascending order, for a given x. Moreover it finds *
 *  the value of the 1st derivative at x.                             *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               degree of the polynomial.                            *
 *      coeff    REAL   *coeff;                                       *
 *               coefficient vector for the polynomial consisting of  *
 *               n + 1 components coeff[0],...,coeff[n].              *
 *      x        REAL   *x;                                           *
 *               value where we want to evaluate the polynomial.      *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      val      REAL   *val;                                         *
 *               Polynomial value at x.                               *
 *      dval     REAL   *dval;                                        *
 *               value of the 1st derivative at x.                    *
 *                                                                    *
 *   Return values:                                                   *
 *   =============                                                    *
 *      = 0      ok                                                   *
 *      = 1      error : n < 0 at the input                           *
 *                                                                    *
 *====================================================================*/
{
  int i;

  if (n < 0) return (1);

  if (n == 0)
  {
    *val = coeff[0];
    *dval = ZERO;
    return (0);
  }
                                         /* Horner scheme             */
  for (*val = *dval = ZERO, i = n; i >= 1; i--)
  {
    *val = *val * x + coeff[i];          /* val  : Polynomial value   */
    *dval = *dval * x + *val;            /* dval : 1st derivative     */
  }

  *val = *val * x + coeff[0];

  return (0);
}
/*.BA*/



/*.FE{C 2.5.3}
     {Newton's Method for Multiple Zeros}
     {Newton's Method for Multiple Zeros;
      a Modified Newton's Method}*/

int newmod              /* Modified Newton Method ....................*/
/*.IX{newmod}*/
           (
            REALFCT   fct,        /* Function ........................*/
            REALFCT   fderv1,     /* 1st derivative ..................*/
            REALFCT   fderv2,     /* 2nd derivative ..................*/
            REAL *    x,          /* Starting value / solution .......*/
            REAL *    fval,       /* Functional value at x ...........*/
            int *     iter,       /* Number of iterations ............*/
            int *     mul         /* Multiplicity of the root ........*/
           )
/*====================================================================*
 *                                                                    *
 *  newmod computes a rooot of a twice continuously differentiable    *
 *  function fct.                                                     *
 *  The function fct, its 1st and 2nd derivative must be explicitly   *
 *  known. If there is some indication that the root of fct may be    *
 *  multiple, we recommend to use this more expensive procedure newmod*
 *  instead of the simple Newton method.                              *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   ============                                                     *
 *      Find roots of the at least twice continuously differentiable  *
 *      function fct, especially designed for multiple roots.         *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      fct      REAL fct (REAL);                                     *
 *               Function, whose root is wanted.                      *
 *               fct has the same form as in newton.                  *
 *      fderv1   REAL fderv (REAL);                                   *
 *               1st derivative of fct of the same form as fct.       *
 *      fderv2   REAL fderv (REAL);                                   *
 *               ditto, the second derivative of fct.                 *
 *      x        REAL   *x;                                           *
 *               Starting value for the iteration.                    *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      x        REAL   *x;                                           *
 *               Computed approximate root of fct.                    *
 *      fval     REAL   *fval;                                        *
 *               Functional value at x. This must be near zero.       *
 *      iter     int *iter;                                           *
 *               Number of iterations performed.                      *
 *      mul      int *mul;                                            *
 *               Multiplicity of the root.                            *
 *                                                                    *
 *   Return value:                                                    *
 *   ============                                                     *
 *      = 0      root found with  abs(fct) < FCTERR                   *
 *      = 1      stop with abs(xnew-xold) < ABSERR + xnew * RELERR    *
 *      = 2      maximum number of iterations reached                 *
 *      = 3      unadmissable function calls                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *                                                                    *
 *      REAL   fct ():   Function, whose root is desired              *
 *      REAL   fderv1(): 1st derivative of fct                        *
 *      REAL   fderv2(): 2nd derivative of fct                        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used: ABSERR, RELERR, MACH_EPS, FCTERR, ITERMAX,       *
 *   ==============  EPSROOT                                          *
 *                                                                    *
 *====================================================================*/
{
  REAL fs, fss, diff, xj = ZERO, y;
  int  rc = 2;

  if (fct == NULL || fderv1 == NULL || fderv2 == NULL)
    return (3);

  *iter = 0;                     /* Initialize iteration counter      */

  while (*iter < ITERMAX)    /* Limit number of iterations to ITERMAX */
  {
    y = (*fct) (*x);                                /* Function value */

    if (ABS (y) < FCTERR)
    {
      rc = 0;
      break;
    }

    fs = (*fderv1) (*x);                            /* 1st derivative */
    fss = (*fderv2) (*x);                           /* 2nd derivative */

    (*iter)++;

    if (ABS (fs) < EPSROOT) fs = SIGN (EPSROOT, fs);  /* divide by 0? */

    xj = ONE / (ONE - y * fss / (fs * fs));

    diff = xj * y / fs;
    *x -= diff;                 /* xnew = xold - xj * fs / 1st deriv. */

    if (ABS (diff) < ABS (*x) * RELERR + ABSERR)
    {
      rc = 1;
      break;
    }
  }

  *fval = y;
  *mul = (int) (xj + HALF);

  return (rc);
}

/* ------------------------- END fnewton.c -------------------------- */
