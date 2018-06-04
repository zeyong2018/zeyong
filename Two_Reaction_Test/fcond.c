#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.15.1}{Error and the Condition Number}{Error and the Condition
                Number}*/

/*.BE*/
/* ------------------------- MODULE fcond.c ------------------------- */

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>

/*.BA*/

REAL hcond              /* Hadamard condition number .................*/
/*.IX{hcond}*/
             (
              int     n,          /* size of matrix ..................*/
              REAL *  mat[]       /* matrix ..........................*/
             )
/*====================================================================*
 *                                                                    *
 *  hcond computes the Hadamard conditions number of an n x n matrix. *
 *  If the return value is order of megnitudes less than 1, then the  *
 *  matrix is ill conditioned. Any solution to a linear system with   *
 *  this matrix will be error prone.                                  *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat                                     *
 *      mat      REAL   *mat[n];                                      *
 *               n x n matrix                                         *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      REAL     < 0.0: error                                         *
 *                 = -1.0 :  n < 1                                    *
 *                 = -2.0 :  lack of memory space                     *
 *                 = -3.0 :  Matrix is singular (det (A) = 0.0)       *
 *                                                                    *
 *               >= 0.0 :                                             *
 *               Hadamard condition number of mat                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use :                                               *
 *   ==================                                               *
 *                                                                    *
 *      int gaudec ():    Factors mat into L*U                        *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :     NULL                                        *
 *   ================                                                 *
 *                                                                    *
 *   Macros: ABS, SQRT                                                *
 *   ======                                                           *
 *====================================================================*/
{
  register j, i;
  REAL     temp, cond, **lu;
  int      rc, signd, *perm;
  void *vmblock;

  if (n < 1) return (-ONE);

                                              /* allocate buffer for  */
  vmblock = vminit();                         /* L * U factorization  */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));

  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return (-TWO);
  }

  rc = gaudec (n, mat, lu, perm, &signd);       /* Compute            */
  if (rc != 0 || signd == 0)                    /* decomposition      */
  {
    vmfree(vmblock);
    return (-THREE);
  }

  cond = ONE;                                   /* compute condition  */
  for (i = 0; i < n; i++)                       /* number             */
  {
    for (temp = ZERO, j = 0; j < n; j++)
      temp += SQR (mat[i][j]);
    cond *= lu[i][i] / SQRT (temp);
  }

  vmfree(vmblock);                                     /* free buffer */

  return (ABS (cond));
}
/*.BA*/



/*.FE{C 4.15.2}{Condition Estimates}{Condition Estimates}*/

/*.FE{}{Condition Estimate according to Cline}
       {Condition Estimate according to Cline}*/

REAL ccond              /* Conditions estimate according to Cline ....*/
/*.IX{ccond}*/
             (
              int     n,          /* Dimension of matrix .............*/
              REAL *  mat[]       /* matrix ..........................*/
             )
/*====================================================================*
 *                                                                    *
 *  ccond estimates the condition number cond (mat) of an  n x n      *
 *  matrix mat according to Cline.                                    *
 *                        -1                                          *
 *  cond (A) = | A | * | A  |, where we choose | | to be the maximum  *
 *  norm.                                                             *
 *                                                                    *
 *  A large value for cond (A) means ill conditioning of A.           *
 *  Consequently solutions of linear systems with A cannot be found   *
 *  with high accuracy.                                               *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat                                     *
 *      mat      REAL   *mat[n];                                      *
 *               n x n matrix                                         *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      REAL     < 0.0: error                                         *
 *                 = -1.0 :  n < 1                                    *
 *                 = -2.0 :  lack of memory                           *
 *                 = -3.0 :  Matrix is singular (det (A) = 0.0)       *
 *                                                                    *
 *               >= 0.0 :                                             *
 *               Cline's estimate of the condition number of mat      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use   :                                             *
 *   ===================                                              *
 *                                                                    *
 *      int gaudec ():    LU factorization of  mat                    *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :  NULL                                           *
 *   ================                                                 *
 *                                                                    *
 *   Macros: ABS                                                      *
 *   ======                                                           *
 *====================================================================*/
{
  register i, j, k;
  REAL     **lu, *x, *y, *z,
           v, smi, spl, sum, xnorm, znorm, matnorm;
  int      rc, signd, *perm;
  void *vmblock;

  if (n < 1) return (-ONE);
                                              /* allocate buffer for  */
  vmblock = vminit();                         /* LU factorization     */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));

  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return (-TWO);
  }

  rc = gaudec (n, mat, lu, perm, &signd);       /* Compute LU         */
  if (rc != 0 || signd == 0)                    /* decomposition      */
  {
    vmfree(vmblock);                            /* free storage       */
    return (-THREE);
  }

  x = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  y = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  z = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))                    /* not enough room ...*/
  {
    vmfree(vmblock);
    return (-THREE);
  }
                                    /* Determine x = (+-1,+-1...,+-1) */
  x[0] = ONE;                       /* so that y is "maximal"         */
  y[0] = ONE / lu[0][0];
  for (i = 1; i < n; i++)
    y[i] = - lu[0][i] * y[0] / lu[i][i];

  for (k = 1; k < n; k++)
  {
    v = ONE / lu[k][k];
    x[k] = y[k] - v;
    y[k] += v;
    smi = ABS (x[k]);
    spl = ABS (y[k]);
    for (i = k + 1; i < n; i++)
    {
      v = lu[k][i] / lu[i][i];
      x[i] = y[i] - v * x[k];
      y[i] -= v * y[k];
      smi += ABS (x[i]);
      spl += ABS (y[i]);
    }

    if (smi > spl)
    {
      for (i = k; i < n; i++) y[i] = x[i];
      x[k] = -ONE;
    }
    else
      x[k] = ONE;
  }

  for (i = n - 1; i >= 0; i--)         /* back substitution ..........*/
  {
    z[i] = y[i];
    for (j = i + 1; j < n; j++)
      z[i] -= lu[j][i] * y[j];
  }

  znorm = ZERO;                        /* find norms .................*/
  xnorm = ZERO;
  matnorm = ZERO;
  for (i = 0; i < n; i++)
  {
    if (ABS (z[i]) > znorm)            /* Max norm of z ..........*/
      znorm = ABS (z[i]);
    if (ABS (x[i]) > xnorm)            /* Max norm of x ..........*/
      xnorm = ABS (x[i]);

    sum = ZERO;
    for (j = 0; j < n; j++)            /* Max norm of mat ........*/
      sum += ABS (mat[i][j]);
    if (sum > matnorm)
      matnorm = sum;
  }

  vmfree(vmblock);                              /* free storage   */

  return (matnorm * (znorm / xnorm));
}
/*.BA*/



/*.FE{}{Condition Estimate according to Forsythe/Moler}
       {Condition Estimate according to Forsythe/Moler}*/

REAL fcond              /* Condition estimate of Forsythe/Moler ......*/
/*.IX{fcond}*/
             (
              int     n,          /* Dimension of matrix .............*/
              REAL *  mat[]       /* matrix ..........................*/
             )
/*====================================================================*
 *                                                                    *
 *  fcond estimates the condition number cond (mat) of an  n x n      *
 *  matrix mat according to Forsythe and Moler.                       *
 *                        -1                                          *
 *  cond (A) = | A | * | A  |, where we choose | | to be the maximum  *
 *  norm.                                                             *
 *                                                                    *
 *  A large value for cond (A) means ill conditioning of A.           *
 *  Consequently solutions of linear systems with A cannot be found   *
 *  with high accuracy.                                               *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat                                     *
 *      mat      REAL   *mat[n];                                      *
 *               n x n matrix                                         *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      REAL     < 0.0: error                                         *
 *                 = -1.0 :  n < 1                                    *
 *                 = -2.0 :  lack of memory                           *
 *                 = -3.0 :  Matrix is singular (det (A) = 0.0)       *
 *                                                                    *
 *               >= 0.0 :                                             *
 *               Forsythe/Moler estimate of condition number of mat   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use   :                                             *
 *   ===================                                              *
 *                                                                    *
 *      int gaudec ():    LU factorization of  mat                    *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :  NULL                                           *
 *   ================                                                 *
 *                                                                    *
 *   Macros: ABS, MACH_EPS                                            *
 *   ======                                                           *
 *====================================================================*/
{
  register i, j;
  REAL     **lu, *x, *b, *r, nom, denom;
  int      rc, signd, *perm;
  LONG_REAL sum;
  void *vmblock;

  if (n < 1) return (-ONE);
                                              /* allocate storage for */
  vmblock = vminit();                         /* LU decomposition     */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));
  x = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  b = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  r = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);

  if (! vmcomplete(vmblock))               /* not enough space        */
  {
    vmfree(vmblock);
    return (-TWO);
  }

  for (i = 0; i < n; i++) b[i] = ONE;      /* b = vector of ones      */

  /* Solve mat * x = b ...............................................*/
  rc = gauss (0, n, mat, lu, perm, b, x, &signd);

  if (rc)                                  /* Matrix singular         */
  {
    vmfree(vmblock);
    return (-THREE);
  }

  for (i = 0; i < n; i++)                  /* find residual vector in */
  {                                        /* double precision        */
    sum = (LONG_REAL) b[i];
    for (j = 0; j < n; j++)
      sum -= (LONG_REAL) mat[i][j] * (LONG_REAL) x[j];
    r[i] = (REAL) sum;
  }

  /* Solve mat * b = r, i.e. perform one step of iterative refinement */
  rc = gauss (2, n, mat, lu, perm, r, b, &signd);

  if (rc != 0 || signd == 0)               /* Should not occur        */
  {
    vmfree(vmblock);
    return (-THREE);
  }

  denom = nom = ZERO;                      /* find max norms          */
  for (i = 0; i < n; i++)
  {
    if (denom < ABS (b[i])) denom = ABS (b[i]);
    if (nom   < ABS (x[i])) nom   = ABS (x[i]);
  }

  vmfree(vmblock);

  return (denom / nom / MACH_EPS);
}

#define MAXITER 30      /* Maximal number of iterations in iterative .*/
                        /* refinement ................................*/
/*.BA*/



/*.FE{C 4.15.4}{Iterative Refinement}{Iterative Refinement}*/

int gausoli              /* Gauss with iterative refinement ..........*/
/*.IX{gausoli}*/
            (
             int     n,            /* Dimension of matrix ............*/
             REAL *  mat[],        /* original matrix ................*/
             REAL *  lumat[],      /* LU decomposition ...............*/
             int     perm[],       /* row interchange vector .........*/
             REAL    b[],          /* Right hand side ................*/
             REAL    x[]           /* solution .......................*/
            )
/*====================================================================*
 *                                                                    *
 *  gausoli solves a linear system lumat * x = b with lumat from      *
 *  gaudec by using gausol and iterative refinement.                  *
 *  Iterative refinement is terminated when the relative improvement  *
 *  does no longer  exceeds  2*MACH_EPS, when the norm of the residual*
 *  vector begins to increase or when the max iteration number MAXITER*
 *  has been reached.                                                 *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of lumat                                   *
 *      mat      REAL   *mat[];                                       *
 *               original system matrix in vector form                *
 *      lumat    REAL   *lumat[];                                     *
 *               LU factorization as supplied by gaudec               *
 *               NOTE: mat and lumat must be stored separately here ! *
 *      perm     int perm[];                                          *
 *               row permutation vector of lumat                      *
 *      b        REAL   b[];                                          *
 *               Right hand side                                      *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      x        REAL   x[];                                          *
 *               solution vector                                      *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      =-1      MAXITER reached, check result carefully              *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or invalid input parameter                     *
 *      = 2      lack of memory                                       *
 *      = 3      invalid LU decomposition (zero diagonal entry)       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Funstions used :                                                 *
 *   ================                                                 *
 *                                                                    *
 *      int gausol ():    determine first solution                    *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants in use  :  NULL, MAXROOT, MACH_EPS, MAXITER            *
 *   ===================                                              *
 *                                                                    *
 *====================================================================*/
{
  int       i, j, k, rc;
  REAL      *r, *z, maxx, maxz, oldmaxz, eps;
  LONG_REAL sumld;
  void *vmblock;

  if (n < 1) return (1);
  if (mat == lumat) return (1);
                             /* Solve system via gaussian elimination */
  if ((rc = gausol (n, lumat, perm, b, x)) != 0)
    return rc;

  eps = (REAL) (TWO * MACH_EPS);
  oldmaxz = MAXROOT;

  vmblock = vminit();
  z = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  r = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
    return 2;

  for (k = 1; ; k++)
  {
    maxx = ZERO;
    for (i = 0; i < n; i++)   /* Double precision for residual vector */
    {
      sumld = (LONG_REAL) b[i];
      for (j = 0; j < n; j++)
        sumld -= (LONG_REAL) mat[i][j] * (LONG_REAL) x[j];

      r[i] = (REAL) sumld;
      if (ABS (x[i]) > maxx) maxx = ABS (x[i]);
    }

    rc = gausol (n, lumat, perm, r, z);    /* solve mat * z = r ......*/
    if (rc) break;

    maxz = ZERO;                 /* improve x, find  max (ABS(z[i]))  */
    for (i = 0; i < n; i++)
    {
      x[i] += z[i];
      if (ABS (z[i]) > maxz) maxz = ABS (z[i]);
    }

    if (maxz < eps * maxx)    /* Check stopping criterion ............*/
    {
      rc = 0;
      break;
    }

    if (oldmaxz < maxz)       /* Increasing residual vector norm ? ...*/
    {
      rc = 0;
      break;
    }

    if (k >= MAXITER)         /* Maximal number of iterations reached */
    {
      rc = -1;
      break;
    }

    oldmaxz = maxz;           /* record maximum norm of z ............*/
  }  /* end of k */

  vmfree(vmblock);            /* free storage ........................*/

  return rc;
}

/* -------------------------- END fcond.c --------------------------- */
