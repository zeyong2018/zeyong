#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.7.1}
     {The Cholesky Decomposition}
     {The Cholesky Decomposition}*/

/*.BE*/
/* ------------------------- MODULE fcholy.c ------------------------ */

#include <basis.h>
#include <u_proto.h>

/*.BA*/

int choly               /* Cholesky Method ...........................*/
/*.IX{choly}*/
          (
           int     mod,           /* Modus: 0, 1, 2 ..................*/
           int     n,             /* Dimension of matrix .............*/
           REAL *  mat[],         /* matrix ..........................*/
           REAL    b[],           /* Right hand side of system .......*/
           REAL    x[]            /* solution vector .................*/
          )
/*====================================================================*
 *                                                                    *
 *  The function cholesky solves linear systems :  mat * x = b        *
 *  for positive definite symmetric n x n matrices mat using the      *
 *  Cholesky method.                                                  *
 *                                                                    *
 *  mat must be symmetric and positive definite, or the method will   *
 *  fail. I.e. for all nonzero vectors y we must have y' * mat * y > 0*
 *                                                                    *
 *  cholesky uses only the lower triangle of mat.                     *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *                                                                    *
 *      Solve linear systems with positive definite system matrices   *
 *      efficeiently.                                                 *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mod      int mod;                                             *
 *       = 0     Factor mat and solve linear system                   *
 *       = 1     Compute factorization only                           *
 *       = 2     Solve system only; for this call the factorization   *
 *               must be stored in mat from chodec. Used to solve for *
 *               many right hand sides.                               *
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat, size of b and x                    *
 *      mat      REAL   *mat[n];                                      *
 *                 mod = 0, 1: Matrix of the linear system            *
 *                 mod = 2   : Cholesky factor                        *
 *      b        REAL   b[n];          ( for mod = 0, 2 )             *
 *               Right hand side of system of equations               *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      mat      REAL   *mat[n];       ( for mod = 0, 1 )             *
 *               Cholesky decomposition of input matrix mat           *
 *      x        REAL   x[n];          ( for mod = 0, 2 )             *
 *               solution vector                                      *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or false input parameter                       *
 *      = 2      Matrix nit positive definite                         *
 *      = 3      wrong call                                           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use  :                                              *
 *   ===================                                              *
 *      int chodec()  : Factorization                                 *
 *      int chosol()  : Solving system                                *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constant used :   NULL                                           *
 *   ===============                                                  *
 *                                                                    *
 *====================================================================*/
{
  int rc;

  if (mat == NULL || n < 1) return (1);

  switch (mod)
  {
    case 0: /* Factor matrix and solve system ........................*/
            rc = chodec (n, mat);
            if (rc == 0)
              return (chosol (n, mat, b, x));
            else
              return (rc);

    case 1: /* factor only ...........................................*/
            return (chodec (n, mat));

    case 2: /* solve only ............................................*/
            return (chosol (n, mat, b, x));
  }

  return (3);                                    /* Wrong call        */
}


int chodec              /* Cholesky decomposition ....................*/
/*.IX{chodec}*/
           (
            int     n,            /* size of matrix ..................*/
            REAL *  mat[]         /* input matrix/Cholesky factor ....*/
           )
/*====================================================================*
 *                                                                    *
 *  chodec decomposes the symmetric positive definite matrix mat.     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat                                     *
 *      mat      REAL   *mat[n];                                      *
 *               Matrix                                               *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      mat      REAL   *mat[n];                                      *
 *               Cholesky decomposition in the lower triangle         *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 1                                                *
 *      = 2      Matrix not  positive definite                        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use   :                                             *
 *   ===================                                              *
 *                                                                    *
 *   From the C library: sqrt()                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants in use :    EPSQUAD                                    *
 *   ===================                                              *
 *                                                                    *
 *====================================================================*/
{
  register j, k, i;
  REAL   sum;

  if (n < 1) return (1);                         /* n < 1  error  !   */
  if (mat == NULL) return (1);
  for (j = 0; j < n; j++)
    if (mat[j] == NULL) return (1);

  if (mat[0][0] < EPSQUAD) return (2);           /* mat not positive  */
                                                 /* definite          */
  mat[0][0] = SQRT (mat[0][0]);
  for (j = 1; j < n; j++) mat[j][0] /= mat[0][0];

  for (i = 1; i < n; i++)
  {
    sum = mat[i][i];
    for (j = 0; j < i; j++)  sum -= SQR (mat[i][j]);

    if (sum < EPSQUAD) return(2);            /* not positive definite */
    mat[i][i] = SQRT (sum);
    for (j = i + 1; j < n; j++)
    {
      sum = mat[j][i];
      for (k = 0; k < i; k++)
        sum -= mat[i][k] * mat[j][k];
      mat[j][i] = sum / mat[i][i];
    }
  }

  return(0);
}


int chosol              /* Cholesky solver ...........................*/
/*.IX{chosol}*/
           (
            int     n,            /* Dimension of matrix .............*/
            REAL *  lmat[],       /* Cholesky matrix .................*/
            REAL    b[],          /* Right hand side of system .......*/
            REAL    x[]           /* solution vector .................*/
           )
/*====================================================================*
 *                                                                    *
 *  chosol finds the solution  x  of the linear system B' *  B * x = b*
 *  for a lower triangular nonsingular matrix B as supplied in chodec.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of lmat, size of b and x                   *
 *      lmat     REAL   *lmat[n];                                     *
 *               lower triangular matrix as supplied by  chodec       *
 *      b        REAL   b[n];                                         *
 *               Right hand side                                      *
 *                                                                    *
 *   Output parameter:                                                *
 *   =================                                                *
 *      x        REAL   x[n];                                         *
 *               solution vector                                      *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      improper lwer triangular matrix or  n < 1            *
 *                                                                    *
 *====================================================================*/
{
  register j, k;
  REAL     sum;

  if (n < 1) return (1);
  if (lmat == NULL || b == NULL || x == NULL) return (1);

  for (j = 0; j < n; j++)
    if (lmat[j] == NULL) return (1);

  if (lmat[0][0] == ZERO) return (1); /* improper factor matrix       */

  x[0] = b[0] / lmat[0][0];          /* update right hand side        */
  for (k = 1; k < n; k++)
  {
    for (sum = ZERO, j = 0; j < k; j++)
      sum += lmat[k][j] * x[j];
    if (lmat[k][k] == ZERO) return (1);
    x[k] = (b[k] - sum) / lmat[k][k];
  }

  x[n-1] /= lmat[n-1][n-1];          /* back substitution             */
  for (k = n - 2; k >= 0; k--)
  {
    for (sum = ZERO, j = k + 1; j < n; j++)
      sum += lmat[j][k] * x[j];
    x[k] = (x[k] - sum) / lmat[k][k];
  }

  return (0);
}

/* --------------------------- END fcholy.c ------------------------- */
