#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.10.1}{Systems with Tridiagonal Matrices}
               {Systems with Tridiagonal Matrices}*/

/*.BE*/
/* ------------------------- MODULE ftdiag.c ------------------------ */

#include <basis.h>
#include <u_proto.h>

/*.BA*/

int trdiag              /* Tridiagonal linear systems ................*/
/*.IX{trdiag}*/
           (
            int     n,            /* size of system matrix ...........*/
            REAL    lower[],      /* lower co-diagonal ...............*/
            REAL    diag[],       /* Diagonal ........................*/
            REAL    upper[],      /* upper co-diagonal ...............*/
            REAL    b[],          /* Right hand side / solution ......*/
            int     rep           /* rep = 0, 1 ......................*/
           )
/*====================================================================*
 *                                                                    *
 *  trdiag solves the linear system  A * x = b for x where A is a     *
 *  tridiagonal n x n matrix, which is given by 3 vectors lower,      *
 *  upper and diag as follows :                                       *
 *                                                                    *
 *       ( diag[0]  upper[0]    0        0  .   .     .   0      )    *
 *       ( lower[1] diag[1]   upper[1]   0      .     .   .      )    *
 *       (   0      lower[2]  diag[2]  upper[2]   0       .      )    *
 *  A =  (   .        0       lower[3]  .     .       .          )    *
 *       (   .          .           .        .     .      0      )    *
 *       (   .             .            .        .      .        )    *
 *       (                   .             .        . upper[n-2] )    *
 *       (   0 .   .    .       0        lower[n-1]   diag[n-1]  )    *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications :                                                   *
 *   ==============                                                   *
 *      Preferred for diagonally dominant tridiagonal matrices, as    *
 *      they occur in spline interpolation for example.               *
 *      Diagonally dominant  matrices always allow an LU factorization*
 *      For non diagonally dominant tridiagonal matrices we recommend *
 *      the funktion band, which uses pivot search and is more stable.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        size of system matrix ( > 1 )  int n                 *
 *                                                                    *
 *      lower    lower co-diagonal         REAL   lower[n]            *
 *      diag     main diagonal             REAL   diag[n]             *
 *      upper    upper co-diagonal         REAL   upper[n]            *
 *                                                                    *
 *               For rep != 0  lower, diag and upper contain the      *
 *               LU factorization of the matrix already.              *
 *                                                                    *
 *      b        right hand side           REAL   b[n]                *
 *      rep      = 0  first call           int rep                    *
 *               !=0  repeated call for same system with new          *
 *                    right hand side                                 *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      b        solution for the system;  REAL   b[n]                *
 *               (the right hand side is lost)                        *
 *                                                                    *
 *      lower    ) if  rep = 0 , after the run, these contain the     *
 *      diag     ) LU factors                                         *
 *      upper    )                                                    *
 *                                                                    *
 *   The determinant of the matrix is given as                        *
 *      det A = diag[0] * ... * diag[n-1] .                           *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 2                                                *
 *      = 2      LU factorization does not exist, matrix singular     *
 *                                                                    *
 *====================================================================*/
{
  register i;

  if (n < 2) return (1);                    /*  n at least 2          */

  if (lower == NULL || diag == NULL || upper == NULL ||
      b == NULL) return (1);

  if (rep == 0)                             /*  for rep = 0, determine*/
  {                                         /*  LU factorization      */
    for (i = 1; i < n; i++)
    {
      if (ABS(diag[i-1]) < MACH_EPS)        /*  if one   diag[i] = 0  */
        return (2);                         /*  we have no LU fact.   */
      lower[i] /= diag[i-1];
      diag[i] -= lower[i] * upper[i-1];
    }
  }

  if (ABS(diag[n-1]) < MACH_EPS) return (2);

  for (i = 1; i < n; i++)                   /*  update  b            */
    b[i] -= lower[i] * b[i-1];

  b[n-1] /= diag[n-1];                      /* backsubsitute         */
  for (i = n-2; i >= 0; i--)
    b[i] = ( b[i] - upper[i] * b[i+1] ) / diag[i];

  return (0);
}

/* --------------------------- END ftdiag.c ------------------------- */
