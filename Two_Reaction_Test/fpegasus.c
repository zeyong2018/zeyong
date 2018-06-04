#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 2.8.2}
     {Pegasus Method}
     {Pegasus Method}*/

/*.BE*/
/* ----------------------- MODULE fpegasus.c ------------------------ */

#include <basis.h>
#include <u_proto.h>

#define ITERMAX 300                   /* Maximal number of iterations */

#define ABSERR ZERO                   /* Admissable absolute error    */
#define RELERR ((REAL)4.0 * MACH_EPS) /* Admissable relative error    */
#define FCTERR EPSQUAD                /* Maximal function errror      */


/*.BA*/

int pegasus             /* Pegasus Method    .........................*/
/*.IX{pegasus}*/
            (
             REALFCT   fct,       /* Function ........................*/
             REAL *    x1,        /* Starting value 1 ................*/
             REAL *    x2,        /* Starting value 2 / solution .....*/
             REAL *    f2,        /* Function value at x2 ............*/
             int *     iter       /* Number of iterations ............*/
            )
/*====================================================================*
 *                                                                    *
 *  pegasus computes one zero of the continuous function fct,         *
 *  provided that the two starting values x1 and x2 satisfy:          *
 *                fct(x1) * fct(x2) <= 0  .                           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Applications:                                                    *
 *   ============                                                     *
 *      Determine one root of the continuous function fct, if an      *
 *      inclusion interval [x1, x2] is known for the root.            *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      fct      REAL   fct (REAL);                                   *
 *               Function, whose root is to be found.                 *
 *               fct has the form:                                    *
 *                   REAL   fct (REAL x)                              *
 *                   { REAL f;                                        *
 *                     f = ...;                                       *
 *                     return(f);                                     *
 *                   }                                                *
 *      x1,x2    REAL   *x1, *x2;                                     *
 *               Starting values with fct(x1) * fct(x2) <= 0.         *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      x2       REAL   *x2;                                          *
 *               Computed approximation for the root of fct           *
 *                                                                    *
 *      f2       REAL   *f2;                                          *
 *               Functional value at the approximate root, this       *
 *               must be nearly zero.                                 *
 *      iter     int *iter;                                           *
 *               Number of iterations that were performed.            *
 *                                                                    *
 *   Return values:                                                   *
 *   =============                                                    *
 *      = -1     No inclusion: fct(x2) * fct(x1) > 0                  *
 *      =  0     Root has been found with ABS(f2) < FCTERR            *
 *      =  1     Break-off with                                       *
 *                   ABS(xnew-xold) < ABSERR + xnew * RELERR,         *
 *               check functional value                               *
 *      =  2     Iteration limit reached                              *
 *      =  3     Improper input parameters                            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used:   ABSERR, RELERR, MACH_EPS, EPSROOT, ITERMAX     *
 *   ==============                                                   *
 *                                                                    *
 *====================================================================*/
{
  REAL f1, x3, f3, s12;
  int  rc = 2;

  *iter = 0;                          /* Initialize iteration counter */

  if (fct == NULL) return (3);

   f1 = (*fct)(*x1);                   /* Function values at *x1, *x2 */
  *f2 = (*fct)(*x2);

  if (f1 * (*f2) > ZERO) return (-1);   /* No inclusion -> Error  */

  if (f1 * (*f2) == ZERO)             /* One starting value is a root */
  {
    if (f1 == ZERO)
    {
      *x2 = *x1;
      *f2 = ZERO;
    }
    return (0);
  }

  while (*iter <= ITERMAX)              /* Pegasus iteration          */
  {
    (*iter)++;

    s12 = (*f2 - f1) / (*x2 - *x1);     /* Secant slope          */

    x3  = *x2 - *f2 / s12;              /* new approximation     */
    f3  = (*fct)(x3);

    if (*f2 * f3 <= ZERO)               /* new inclusion interval */
    {
      *x1 = *x2;
       f1 = *f2;
    }
    else
      f1 *= *f2 / ( *f2 + f3 );

    *x2 = x3;
    *f2 = f3;

    if ( ABS(*f2) < FCTERR )            /* Root found        */
    {
      rc = 0;
      break;
    }
                                    /* Break-off with small step size */
    if ( ABS(*x2 - *x1) <= ABS(*x2) * RELERR + ABSERR )
    {
      rc = 1;
      break;
    }
  }

  if ( ABS(f1) < ABS(*f2) )      /* Choose approximate root with    */
  {                              /* least magnitude function value  */
    *x2 = *x1;
    *f2 = f1;
  }

  return (rc);
}

/* ------------------------- ENDE fpegasus.c ------------------------ */
