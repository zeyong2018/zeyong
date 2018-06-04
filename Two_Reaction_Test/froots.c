#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 2.8.4}
     {The King and the Anderson-Bj"orck-King Method}
     {The King and the Anderson-Bj"orck-King Methods,
      the Illinois Method}*/

/*.BE*/
/* ------------------------ MODULE froots.c ------------------------- */

#include <basis.h>
#include <u_proto.h>

#define ITERMAX 300                /* Maximal number of function      */
                                   /* evaluations allowed             */
#define ABSERR (REAL)0.0           /* Allowable absolute error        */
#define RELERR ((REAL)4.0 * MACH_EPS)   /* Allowable relative error   */
#define FCTERR EPSQUAD             /* Maximal function error          */


static void swap2 (REAL * x, REAL * y, REAL * fx, REAL * fy);

/*.BA*/

int roots               /* Pegasus, Anderson-Bjoerck-King method    ..*/
/*.IX{roots}*/
          (
           int       method,      /* method    .......................*/
           REALFCT   fct,         /* function ........................*/
           int       quadex,      /* quadratic extrapolation   .......*/
           REAL *    x1,          /* Starting value 1 ................*/
           REAL *    x2,          /* Starting value 2 / solution .....*/
           REAL *    fx2,         /* Function value at x2 ............*/
           int *     iter         /* Iteration number.................*/
          )
/*====================================================================*
 *                                                                    *
 *  The function roots computes a zero of the continuous function fct.*
 *  If the starting values  x1 and x2 satisfy inclusion:              *
 *  fct(x1) * fct(x2) <= 0.0, then each of the four methods           *
 *  (Pegasus, Pegasus-King, Anderson-Bjoerck, Anderson-Bjoerck-King)  *
 *  is convergent.                                                    *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   ============                                                     *
 *      Find a root of the continuous function fct. If an interval    *
 *      of inclusion [x1, x2] is known for the root, then the methods *
 *      are surely and always convergent.                             *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   INPUT PARAMETERS:                                                *
 *   ================                                                 *
 *      method   int method;                                          *
 *               Chosen method:                                       *
 *          =1   Pegasus method                                       *
 *          =2   Pegasus-King method                                  *
 *          =3   Anderson-Bjoerck method                              *
 *       sonst   Anderson-Bjoerck-King method                         *
 *      fct      REAL   fct (REAL);                                   *
 *               Function, whose root is desired.                     *
 *               fct has the form:                                    *
 *                   REAL   fct (REAL x)                              *
 *                   { REAL   f;                                      *
 *                     f = ...;                                       *
 *                     return(f);                                     *
 *                   }                                                *
 *      quadex   int quadex;                                          *
 *          =0   we allow only linear extrapolation, if then starting *
 *               values do not include the root;                      *
 *               this is advisable, if a multiple root is suspected   *
 *       sonst   quadratic extrapolation.                             *
 *      x1,x2    REAL   *x1, *x2;                                     *
 *               Starting values for the iteration; if only one       *
 *               starting value is known, one may set *x1 = *x2.      *
 *               The method will automatically create a second        *
 *               starting value.                                      *
 *                                                                    *
 *   OUTPUT PARAMETERS:                                               *
 *   =================                                                *
 *      x2       REAL   *x2;                                          *
 *               Computed approximation for the root of fct.          *
 *      f2       REAL   *f2;                                          *
 *               Functional value at x2. This must be close to zero.  *
 *      iter     int *iter;                                           *
 *               Number of iterations carried out.                    *
 *                                                                    *
 *   Return values:                                                   *
 *   =============                                                    *
 *      =  0     Zero found with ABS(f2) < FCTERR                     *
 *      =  1     Break-off with                                       *
 *                    ABS(xnew-xold) < ABSERR + xnew * RELERR,        *
 *               check function value                                 *
 *      =  2     Maximum number of iterations has been reached        *
 *      =  3     Improper input function                              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *                                                                    *
 *      swap2 (): Interchanges two  x  and two function values        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used: ABSERR, RELERR, MACH_EPS, EPSROOT, FCTERR,       *
 *   ==============  ITERMAX                                          *
 *                                                                    *
 *====================================================================*/
{
  REAL  f1, x3, f3, fquot, q, x3new;
  int   sec = 0, neg = 0, incl, rc = 2;

  if (fct == NULL) return (3);

  if (*x1 == *x2)                       /* If *x1 = *x2, change *x2  */
    *x2 = (REAL) ((ONE + EPSROOT) * (*x2) + EPSROOT);

    f1 = (*fct) (*x1);       /* Function value at the starting points */
  *fx2 = (*fct) (*x2);
                             /* Initialize counter for number of      */
  *iter = 2;                 /* function evaluations.                 */


  if ( ABS (f1) < ABS (*fx2) )       /* Store function value of least */
    swap2 (x1, x2, &f1, fx2);        /* magnitude in f2               */

  if ( ABS (*fx2) < FCTERR )         /* Is *x2 a root ?               */
    return (0);

  if (f1 * (*fx2) > ZERO)            /* if incl = 0: no inclusion,    */
  {
    incl = 0;
    f3 = *fx2;
  }
  else                               /* else: inclusion               */
  {
    incl = 1;                        /* if sec = 1: Next step is a    */
    sec = 1;                         /* secant step                   */
  }


  while ( *iter <= ITERMAX )         /* Iteration start  .............*/
  {
    if (!incl)                       /* If *x1, *x2 do not include    */
    {
      fquot = f1 / *fx2;             /* the root                     */
      if (fquot > ONE)
      {
        if ( quadex && (fquot - f1 / f3 > ONE) )
          f1 *= ONE - *fx2 / f3;
      }
      else
      if (fquot < ZERO)
      {
        incl = 1;
        if ( ABS (*x1 - *x2) <= ABS (*x2) * RELERR + ABSERR )
        {
          rc = 1;
          break;
        }
        else sec = 1;
      }
      else return (2);                   /* No root found             */
    }

    q = *fx2 / (*fx2 - f1);
    x3 = *x2 + q * (*x1 - *x2);

    if (!incl)             /* If there is no inclusion, we construct  */
      if ( *x2 == x3 )     /* a x3, different from both  *x2 and *x1. */
                           /* If there is no such x3, then the root   */
      {                    /* is at *x2.                              */
        x3new = *x2 + (*x1 - *x2) / THREE;
        if (x3new == *x2)
        {
          rc = 1;
          break;
        }
        else
        {
          do
          {
            q += q;
            x3new = *x2 + q * (*x1 - *x2);
          }
          while (x3new == *x2);

          if (x3new == *x1)
          {
            rc = 1;
            break;
          }
          else x3 = x3new;
        }
      }
      else
      if ( *x1 == x3 )
      {
        x3new = *x1 + (*x2 - *x1) / THREE;
        if (x3new == *x1)
        {
          rc = 1;
          break;
        }
        else
        {
          q = f1 / (f1 - *fx2);
          do
          {
            q += q;
            x3new = *x1 + q * (*x2 - *x1);
          }
          while (x3new == *x1);

          if (x3new == *x2)
          {
            rc = 1;
            break;
          }
          else x3 = x3new;
        }
      }
                             /* now we have: x1 != x3 and *x2 != x3.  */
      f3 = (*fct) (x3);      /* Compute the function value at  x3.    */
      (*iter)++;             /* Up the function evaluation counter.   */

      if ( ABS(f3) < FCTERR )                     /* Is x3 a root  ?  */
      {
        rc = 0;
        swap2 (x2, &x3, fx2, &f3);                /* swap and done    */
        return (0);
      }

      if (!incl) swap2 (x1, x2, &f1, fx2);
      else
      if ( *fx2 * f3 < ZERO )
      {
        neg = 1;                       /* If neg = 1: *x2 and x3      */
        swap2 (x1, x2, &f1, fx2);      /* include the root            */
      }
      else neg = 0;

      swap2 (x2, &x3, fx2, &f3);

      if (!incl) continue;     /* no inclusion => continue with while */
                               /* else: Are *x1 or *x2 roots ?        */

      if ( ABS (*x1 - *x2) <= ABS (*x2) * RELERR + ABSERR )
      {
        rc = 1;
        break;
      }

      /* Otherwise use the chosen method to find a new value for f1:  */

      switch (method)
      {
        case 1:  if (!neg)         /* Pegasus method                  */
                   f1 *= f3 / (*fx2 + f3);
                 break;

        case 2:  if (sec)          /* Pegasus-King method             */
                 {
                   f1 *= f3 / (*fx2 + f3);
                   sec = 0;
                 }
                 else if ( !neg )
                        f1 *= f3 / (*fx2 + f3);
                      else sec = 1;
                 break;

        case 3:  if (!neg)         /* Anderson-Bjoerck method         */
                 {
                   q = ONE - *fx2 / f3;
                   if (q <= ZERO) q = HALF;
                   f1 *= q;
                 }
                 break;

        default: if (sec)          /* Anderson-Bjoerck-King method    */
                 {
                   q = ONE - *fx2 / f3;
                   f1 *= (q > ZERO) ? q : HALF;
                   sec = 0;
                 }
                 else if ( !neg )
                 {
                   q = ONE - *fx2 / f3;
                   f1 *= (q > ZERO) ? q : HALF;
                 }
                 else sec = 1;
                 break;
      }

  }  /* ende while ( *iter < ITERMAX ) */

  if ( ABS (f1) < ABS (*fx2) )    /* Store least magnitude function   */
    swap2 (x1, x2, &f1, fx2);     /* value and corresponding x value  */
                                  /* in *fx2 and *x2.                 */
  return (rc);
}


static void swap2 (REAL * x,
                   REAL * y,
                   REAL * fx,
                   REAL * fy)
/*====================================================================*
 *                                                                    *
 *  The procedure swap2 exchanges the values stored in x,y and those  *
 *  stored in fx,fy simultaneously.                                   *
 *                                                                    *
 *====================================================================*/
{
  REAL tmp;

  tmp = *x;    *x = *y;    *y = tmp;
  tmp = *fx;  *fx = *fy;  *fy = tmp;
}

/* -------------------------- END froots.c -------------------------- */
