#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/
/*.KA{C 7}{Eigenvalues and Eigenvectors of Matrices}
          {Eigenvalues and Eigenvectors of Matrices}*/
/*.FE{C 7.3.1}
     {The Dominant Eigenvalue and the Associated Eigenvector}
     {The Dominant Eigenvalue and the Associated Eigenvector
      of a Matrix}*/

/*.BE*/
/* ------------------------ MODULE fmises.c ------------------------- */

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>

#define MAXIT 200
#define eps (REAL)256.0 * MACH_EPS

/*.BA*/

int mises               /* Vector iteration for max modulus eigenvalue*/
/*.IX{mises}*/
          (
           int     n,             /* Dimension of matrix .............*/
           REAL *  mat[],         /* matrix ..........................*/
           REAL    x[],           /* Eigenvector .....................*/
           REAL *  ew             /* maximum modulus eigenvalue ......*/
          )
/*====================================================================*
 *                                                                    *
 *  The function mises determines the maximal modulus eigenvalue of a *
 *  matrix and the corresponding eigenvector by vector iteration.     *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *      For real  n x n matrices, with a single dominant eigenvalue   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 1 )                                    *
 *               Dimension of  mat                                    *
 *      mat      REAL   *mat[n];                                      *
 *               input matrix                                         *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      x        REAL   x[n];                                         *
 *               Eigenvector of mat for the maximal modulus evalue    *
 *      ew       REAL   *ew;                                          *
 *               maximum modulus eigenvalue of mat                    *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all is ok                                            *
 *      = 1      n < 2                                                *
 *      = 2      lack of memory space                                 *
 *      = 3      maximal number of iterations reached                 *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   functions in use  :                                              *
 *   ===================                                              *
 *                                                                    *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors or matrices            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants in use  :  NULL, MACH_EPS, MAXIT                       *
 *   ===================                                              *
 *                                                                    *
 *   Macros: SQRT                                                     *
 *   ======                                                           *
 *====================================================================*/
{
  int  i, j, iter;
  REAL *y, s, tmp;
  void *vmblock;

  if (n < 2) return (1);                    /*  n must be  > 1        */
  vmblock = vminit();                       /*  allocate buffers      */
  y = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
    return 2;

  s = ONE / SQRT ((double) n);              /*  initialize x with     */
  for (i = 0; i < n; i++)                   /*  norm 1                */
    x[i] = s;

  for (iter = 1; iter <= MAXIT; iter++)     /*  Iteration             */
  {
    *ew = ZERO;
    tmp = ZERO;
    for (i = 0; i < n; i++)                 /*  compute mat * x       */
    {
      for (y[i] = ZERO, j = 0; j < n; j++)
        y[i] += mat[i][j] * x[j];
      *ew += x[i] * y[i];                  /*  compute  x * mat * x   */
      tmp += y[i] * y[i];                  /*  and norm y             */
    }

    if (tmp == ZERO)                       /*  mat = zero matrix      */
    {
      vmfree(vmblock);
      return (0);
    }

    tmp = SQRT (tmp);
    for (s = ZERO, i = 0; i < n; i++)      /*  compute norm of        */
    {                                      /*        mat * x - ew * x */
      s += (*ew * x[i] - y[i]) * (*ew * x[i] - y[i]);
      x[i] = y[i] / tmp;                   /*  set up x for next      */
    }                                      /*  iteration              */
    if (SQRT (s) < (REAL) (eps * (*ew)))   /*  accurate enough ?      */
    {
      vmfree(vmblock);
      return (0);
    }
  }

  vmfree(vmblock);
  return (3);
}

/* --------------------------- END fmises.c ------------------------- */
