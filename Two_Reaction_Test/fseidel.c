#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/
/*.KA{C 5}{Iterative Methods for Linear Systems}
          {Iterative Methods for Linear Systems}*/
/*.FE{C 5.4}
     {The Gau"s-Seidel Iteration}{The Gau"s-Seidel Iteration}*/

/*.BE*/
/* ------------------------ MODULE fseidel.c ------------------------ */

#include <basis.h>
#include <u_proto.h>

#define ITERMAX 300                   /* Maximal number of iterations */

/*.BA*/

int seidel              /* Gauss Seidel Method with relaxation .......*/
/*.IX{seidel}*/
           (
            int     crit,         /* crit = 0, 1, 2, 3 ...............*/
            int     n,            /* size of matrix ..................*/
            REAL *  mat[],        /* matrix ..........................*/
            REAL    b[],          /* Right hand side .................*/
            REAL    omega,        /* Relaxaktion coefficient .........*/
            REAL    x[],          /* solution ........................*/
            REAL    residu[],     /* Residuum vector .................*/
            int *   iter          /* # of iterations .................*/
           )
/*====================================================================*
 *                                                                    *
 *  seidel solves the linear system  mat * x = b  iteratively.        *
 *  Here  mat  is a nonsingular  n x n  matrix, b is the right hand   *
 *  side for teh linear system and x is the solution.                 *
 *                                                                    *
 *  seidel uses the Gauss Seidel Method with relaxation for a given   *
 *  relaxation coefficient 0 < omega < 2.                             *
 *  If  omega = 1, the standard Gauss Seidel method (without          *
 *  relaxation) is performed.                                         *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   =============                                                    *
 *      Solve linear systems with nonsingular system matrices that    *
 *      satisfy one of the following criteria: row sum criterion,     *
 *      column sum criterion or the criterion of Schmidt and v. Mises.*
 *      Only if at least one of these criteria is satisfied for mat,  *
 *      convergence of the scheme is guaranteed.                      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      crit     int crit;                                            *
 *               select criterion                                     *
 *               =1 : row sum criterion                               *
 *               =2 : column sum criterion                            *
 *               =3 : criterion of Schmidt-v.Mises                    *
 *               other : no check                                     *
 *      n        int n;  ( n > 0 )                                    *
 *               size of mat, b and x                                 *
 *      mat      REAL   *mat[n];                                      *
 *               Matrix of the liear system                           *
 *      b        REAL   b[n];                                         *
 *               Right hand side                                      *
 *      omega    REAL   omega; ( 0.0 < omega < 2.0 )                  *
 *               Relaxation coefficient.                              *
 *      x        REAL   x[n];                                         *
 *               Starting vector for iteration                        *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      x        REAL   x[n];                                         *
 *               solution vector                                      *
 *      residu   REAL   residu[n];                                    *
 *               residual vector  b - mat * x; close to zero vector   *
 *      iter     int *iter;                                           *
 *               Number of iterations performed                       *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      =  0     solution has been found                              *
 *      =  1     n < 1  or omega <= 0 or omega >= 2                   *
 *      =  2     improper mat or b or x                               *
 *      =  3     one diagonal element of mat vanishes                 *
 *      =  4     Iteration number exceeded                            *
 *      = 11     column sum criterion violated                        *
 *      = 12     row sum criterion violated                           *
 *      = 13     Schmidt-v.Mises criterion violated                   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :     NULL, MACH_EPS                              *
 *   ================                                                 *
 *                                                                    *
 *   Macros :         ABS, SQR, SQRT                                  *
 *   ========                                                         *
 *                                                                    *
 *====================================================================*/
{
  int   i, j, rc = 0;
  REAL  tmp, eps;

  *iter = 0;                          /* Initialize iteration counter */

  if ( n < 1 ||                           /* Check omega              */
       omega <= ZERO || omega >= TWO ) return (1);

  if (mat == NULL) return (2);

  for (i = 0; i < n; i++)
    if (mat[i] == NULL) return (2);

  if (x == NULL || b == NULL || residu == NULL) return (2);

  eps = (REAL) (MACH_EPS * (REAL) 128.0);

  for (i = 0; i < n; i++)                 /* transform mat so that all*/
  {                                       /* diagonals equal 1        */
    if (mat[i][i] == ZERO) return (3);
    tmp = ONE / mat[i][i];
    for (j = 0; j < n; j++)
      mat[i][j] *= tmp;
    b[i] *= tmp;                 /* adjust right hand side  b         */
  }

  switch (crit)                  /* check convergence criteria        */
  {

    case 1: for (i = 0; i < n; i++)         /* row sum criterion      */
            {
              for (tmp = ZERO, j = 0; j < n; j++)
                tmp += ABS (mat[i][j]);
              if (tmp >= TWO) return (11);
            }
            break;

    case 2: for (j = 0; j < n; j++)         /* column sum criterion   */
            {
              for (tmp = ZERO, i = 0; i < n; i++)
                tmp += ABS (mat[i][j]);
              if (tmp >= TWO) return (12);
            }
            break;

    case 3: for (tmp = ZERO, i = 0; i < n; i++)
              for (j = 0; j < n; j++)        /* criterion of Schmidt, */
                tmp += SQR (mat[i][j]);      /* v. Mises              */
            tmp = SQRT (tmp - ONE);
            if (tmp >= ONE) return (13);
            break;

    default: break;                          /* do not check          */
  }

  for (i = 0; i < n; i++) residu[i] = x[i];  /* store x in residu     */

  while (*iter <= ITERMAX)                   /* Begin iteration       */
  {
    (*iter)++;

    for (i = 0; i < n; i++)
    {
      for (tmp = b[i], j = 0; j < n; j++)
        tmp -= mat[i][j] * residu[j];
      residu[i] += omega * tmp;
    }

    for (i = 0; i < n; i++)          /* check break-off criterion     */
    {
      tmp = x[i] - residu[i];
      if (ABS (tmp) <= eps)
      {
        x[i] = residu[i];            /* If rc = 0 at end of loop      */
        rc = 0;                      /*  -> stop iteration            */
      }
      else
      {
        for (j = 0; j < n; j++) x[j] = residu[j];
        rc = 4;
        break;
      }
    }
    if (rc == 0) break;                        /* solution found      */
  }                                            /* End iteration       */

  for (i = 0; i < n; i++)                     /* find residual vector */
  {
    for (tmp = b[i], j = 0; j < n; j++)
      tmp -= mat[i][j] * x[j];
    residu[i] = tmp;
  }

  return (rc);
}

/* ------------------------- END fseidel.c -------------------------- */
