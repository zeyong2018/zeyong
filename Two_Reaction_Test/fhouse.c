#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.14}
     {Solving Linear Systems via Householder Transformations}
     {Solving Linear Systems via Householder Transformations}*/

/*.BE*/
/* ------------------------- MODULE fhouse.c ------------------------ */

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>

/*.BA*/

int house               /* Householder Method ........................*/
/*.IX{house}*/
          (
           int     m,             /* # of rows .......................*/
           int     n,             /* # of columns ....................*/
           REAL *  mat[],         /* Input matrix ....................*/
           REAL    b[]            /* righ thand side/solution ........*/
          )
/*====================================================================*
 *  The function house solves a linear system of equations            *
 *            mat * x = b.                                            *
 *  Here m is the number of rows in the matrix mat, n the number of   *
 *  its columns with  m >= n and  rank (mat) = n.                     *
 *  b is the right hand side, an m vector, and x is the solution      *
 *  n vector.                                                         *
 *                                                                    *
 *  house uses Householder transformations to solve the overdetermined*
 *  linear system. x is the solution of the minimization problem for  *
 *   mat * x - b in the euclidean norm.                               *
 *  This solution need not solve the linear system too well.          *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      m        int m;  ( m > 0 )                                    *
 *               number of rows in  mat                               *
 *      n        int n;  ( n > 0 )                                    *
 *               number of columns in mat                             *
 *      mat      REAL   *mat[]:                                       *
 *               system matrix mat:                                   *
 *                    mat[i][j], i = 0, ..., m-1, j = 0, ..., n-1.    *
 *      b        REAL   b[m];                                         *
 *               Right hand side                                      *
 *                                                                    *
 *    mat and b are over-written during the computations !            *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      b        REAL   b[n];                                         *
 *               solution vector                                      *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      m or n < 1 or m < n                                  *
 *      = 2      mat does not have maximal rank n.                    *
 *      = 3      Matrix numerically has rank < n                      *
 *      = 4      lack of memory                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use :                                               *
 *   ===================                                              *
 *                                                                    *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants in use  :  NULL, MACH_EPS                              *
 *   ===================                                              *
 *                                                                    *
 *   Macros: SQRT                                                     *
 *   =======                                                          *
 *                                                                    *
 *====================================================================*/
{
  int i, j, k;
  REAL f, r, alpha, ak, eps, tmp, sum,
         norm, maxnorm;
  REAL *d;
  void *vmblock;

  if ((m < 1) || (n < 1) || (m < n)) return (1);

  eps = (REAL)(TWO * MACH_EPS);

  vmblock = vminit();
  d = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
    return 4;

  for (i = 0; i < n; i++)            /*  Householder transformation   */
  {
    r = ZERO;
    for (k = i; k < m; k++)
      r += mat[k][i] * mat[k][i];

    if (r == ZERO)                   /*  Matrix not of full rank      */
    {
      vmfree(vmblock);
      return (2);
    }

    if (mat[i][i] >= ZERO)
      alpha = SQRT (r);
    else
      alpha = - SQRT (r);

    ak = ONE / (r + alpha * mat[i][i]);
    mat[i][i] += alpha;

    d[i] = - alpha;

    maxnorm = ZERO;
    for (k = i + 1; k < n; k++)
    {
      norm = f = ZERO;
      for (j = i; j < m; j++)
      {
        tmp = mat[j][k];
        f += tmp * mat[j][i];
        norm += tmp * tmp;
      }

      if (norm > maxnorm)
        maxnorm = norm;

      f *= ak;
      for (j = i; j < m; j++)
        mat[j][k] -= f * mat[j][i];
    }

    if (ABS(alpha) < eps * SQRT (maxnorm))     /* solvable ?          */
    {
      vmfree(vmblock);
      return (3);
    }

    for (f = ZERO, j = i; j < m; j++)      /* update right hand side  */
      f += b[j] * mat[j][i];

    f *= ak;
    for (j = i; j < m; j++)
      b[j] -= f * mat[j][i];

  }  /*  for i  */

  for (i = n - 1; i >= 0; i--)             /* back substitution       */
  {
    sum = ZERO;
    for (k = i + 1; k < n; k++)
      sum += mat[i][k] * b[k];

    b[i] = (b[i] - sum) / d[i];
  }

  vmfree(vmblock);

  return (0);
}


int mhouse              /* Householder method for m right hand sides .*/
/*.IX{mhouse}*/
           (
            int     m,            /* # of rows .......................*/
            int     n,            /* # of columns ....................*/
            int     k,            /* # right hand sides ..............*/
            REAL *  mat[],        /* matrix ..........................*/
            REAL *  xmat[]        /* Right hand sides/solutions ......*/
           )
/*====================================================================*
 *  The function mhouse solves a linear system for several right hand *
 *  sides:  mat * xmat = rs. (for details see house)                  *
 *  k denotes the number of right hand sides, rs is m x k; xmat is    *
 *  n x k, where each column denotes one solution vector.             *
 *                                                                    *
 *  mhouse uses house to solve each linear system in the minimization *
 *  sense. The individual solution vectors may not be exact solutions *
 *  to the linear overdetermined system.                              *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int m;  ( m > 0 )                                    *
 *               number of rows in mat                                *
 *      n        int n;  ( n > 0 )                                    *
 *               number of columns in mat                             *
 *      k        int k;  ( k > 0 )                                    *
 *               number of right hand sides in rs                     *
 *      mat      REAL   *mat[];                                       *
 *               Matrix :                                             *
 *                 mat[i][j], i = 0,...,m-1, j = 0,...,n-1.           *
 *      xmat     REAL   *xmat[];                                      *
 *               matrix of right hand sides :                         *
 *                 xmat[i][j], i = 0,...,m-1, j = 0,...,k-1.          *
 *                                                                    *
 *    mat and xmat are over-written during computations !             *
 *                                                                    *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      xmat     REAL   *xmat[];                                      *
 *               matrix of solution vectors                           *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      m or n < 1 or m < n or k < 1                         *
 *      = 2      mat not of maximal rank                              *
 *      = 3      Matrix not of maximal rank numerically               *
 *      = 4      lack of memory space                                 *
 *                                                                    *
 *                                                                    *
 *   Functions in use  :                                              *
 *   ===================                                              *
 *                                                                    *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :  NULL, MACH_EPS                                 *
 *   ================                                                 *
 *                                                                    *
 *   Macros: SQRT                                                     *
 *   =======                                                          *
 *                                                                    *
 *====================================================================*/
{
  int  i, j, k0;
  REAL f, r, alpha, ak, eps, tmp, sum, norm, maxnorm;
  REAL *d;
  void *vmblock;

  if ((m < 1) || (n < 1) || (k < 1) || (m < n))
    return (1);

  eps = (REAL)(TWO * MACH_EPS);

  vmblock = vminit();
  d = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
    return 4;

  for (i = 0; i < n; i++)            /* Householder transformation    */
  {
    r = ZERO;
    for (k0 = i; k0 < m; k0++)
      r += mat[k0][i] * mat[k0][i];

    if (r == ZERO)                   /* Matrix not of maximal rank    */
    {
      vmfree(vmblock);
      return (2);
    }

    if (mat[i][i] >= ZERO)
      alpha = SQRT (r);
    else
      alpha = - SQRT (r);

    ak = ONE / (r + alpha * mat[i][i]);
    mat[i][i] += alpha;

    d[i] = - alpha;

    maxnorm = ZERO;                                 /* check rank     */
    for (k0 = i + 1; k0 < n; k0++)
    {
      norm = f = ZERO;
      for (j = i; j < m; j++)
      {
        tmp = mat[j][k0];
        f += tmp * mat[j][i];
        norm += tmp * tmp;
      }
      if (norm > maxnorm)
        maxnorm = norm;

      f *= ak;
      for (j = i; j < m; j++)
        mat[j][k0] -= f * mat[j][i];
    }

    if (ABS(alpha) < eps * SQRT (maxnorm))        /*  no solution ?  */
    {
      vmfree(vmblock);
      return (3);
    }

    for (k0 = 0; k0 < k; k0++)           /* update right hand sides  */
    {
      f = ZERO;
      for (j = i; j < m; j++)
        f += xmat[j][k0] * mat[j][i];
      f *= ak;
      for (j = i; j < m; j++)
        xmat[j][k0] -= f * mat[j][i];
    }
  }  /*  for i  */

  for (j = 0; j < k; j++)            /*  Loop over k right hand sides */
    for (i = n - 1; i >= 0; i--)               /*  back substitution  */
    {
      for (sum = ZERO, k0 = i + 1; k0 < n; k0++)
        sum += mat[i][k0] * xmat[k0][j];

      xmat[i][j] = (xmat[i][j] - sum) / d[i];
    }

  vmfree(vmblock);

  return (0);
}

/* -------------------------- END fhouse.c -------------------------- */
