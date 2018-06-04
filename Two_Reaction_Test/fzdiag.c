#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.11.1}{Systems with a Cyclically Tridiagonal Matrix}
               {Systems with a Cyclically Tridiagonal Matrix}*/

/*.BE*/
/* ------------------------- MODULE fzdiag.c ------------------------ */

#include <basis.h>
#include <u_proto.h>

/*.BA*/

int tzdiag              /* cyclic tridiagonal linear systems .........*/
/*.IX{tzdiag}*/
           (
            int   n,              /* size of matrix ..................*/
            REAL  lower[],        /* Sub-diagonal ....................*/
            REAL  diag[],         /* Diagonal ........................*/
            REAL  upper[],        /* Super-diagonal ..................*/
            REAL  lowrow[],       /* row below .......................*/
            REAL  ricol[],        /* column to the right .............*/
            REAL  b[],            /* right hand side, or solution ....*/
            int   rep             /* rep = 0, 1 ......................*/
           )
/*====================================================================*
 *                                                                    *
 *  tzdiag computes the solution  x  of the linear syatem             *
 *  A * x = b where A is cyclic and tridiagonal. A is given via 5     *
 *  the vectors lower, upper, diag, lowrow and ricol :                *
 *                                                                    *
 *       ( diag[0]  upper[0]    0        0  .   . 0   ricol[0]   )    *
 *       ( lower[1] diag[1]   upper[1]   0      .     .   0      )    *
 *       (   0      lower[2]  diag[2]  upper[2]   0       .      )    *
 *  A =  (   .        0       lower[3]  .     .       .   .      )    *
 *       (   .          .           .        .     .      0      )    *
 *       (   .             .            .        .      .        )    *
 *       (   0               .             .        . upper[n-2] )    *
 *       ( lowrow[0]  0 .  .  .  0       lower[n-1]   diag[n-1]  )    *
 *                                                                    *
 *  Additional storage for lowrow[1],..,lowrow[n-3] and ricol[1],..., *
 *  ricol[n-3] must be supplied, as these values will appear in the   *
 *  decomposition which will be stored in these 5 vectors.            *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   =============                                                    *
 *      For diagonally dominant cyclic tridiagonal linear systems as  *
 *      they appear in spline interpolation.                          *
 *      Diagonally dominant matrices always allow an LU factorization.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        Dimension of matrix  ( > 2 )  int n                  *
 *      lower    lower co-diagonal             REAL   lower[n]        *
 *      diag     main diagonal                 REAL   diag[n]         *
 *      upper    upper co-diagonal             REAL   upper[n]        *
 *      b        right hand side               REAL   b[n]            *
 *      rep      = 0  first call               int rep                *
 *               !=0  subsequant call for same matrix, different b    *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      b        solution of system            REAL   b[n]            *
 *                                                                    *
 *      lower    ) if  rep = 0 ,these vectors contain the             *
 *      diag     ) factorization. The original matrix is lost.        *
 *      upper    )                                                    *
 *      lowrow   )                             REAL   lowrow[n-2]     *
 *      ricol    )                             REAL   ricol[n-2]      *
 *                                                                    *
 *   If  rep = 0, the determinant of A is computable as :             *
 *      det A = diag[0] * ... * diag[n-1]  .                          *
 *                                                                    *
 *   Return values:                                                   *
 *   ==============                                                   *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other invalid input parameter               *
 *      = 2      LU decomposition does not exist                      *
 *                                                                    *
 *====================================================================*/
{
  REAL     tmp;
  register i;

  if (n < 3) return (1);
  if (lower == NULL || diag == NULL || upper == NULL ||
      lowrow == NULL || ricol == NULL) return (1);

  if (rep == 0)                             /*  if rep = 0,           */
  {                                         /*  decompose matrix      */
    lower[0] = upper[n-1] = ZERO;

    if (ABS (diag[0]) < MACH_EPS) return (2);
                                          /* If a diagonal entry is   */
    tmp = ONE / diag[0];                  /* too close to zero, stop  */
    upper[0] *= tmp;
    ricol[0] *= tmp;

    for (i = 1; i < n-2; i++)
    {
      diag[i] -= lower[i] * upper[i-1];
      if (ABS(diag[i]) < MACH_EPS) return (2);
      tmp = ONE / diag[i];
      upper[i] *= tmp;
      ricol[i] = -lower[i] * ricol[i-1] * tmp;
    }

    diag[n-2] -= lower[n-2] * upper[n-3];
    if (ABS(diag[n-2]) < MACH_EPS) return (2);

    for (i = 1; i < n-2; i++)
      lowrow[i] = -lowrow[i-1] * upper[i-1];

    lower[n-1] -= lowrow[n-3] * upper[n-3];
    upper[n-2] = ( upper[n-2] - lower[n-2] * ricol[n-3] ) / diag[n-2];

    for (tmp = ZERO, i = 0; i < n-2; i++)
      tmp -= lowrow[i] * ricol[i];
    diag[n-1] += tmp - lower[n-1] * upper[n-2];

    if (ABS (diag[n-1]) < MACH_EPS) return (2);
  }  /* end if ( rep == 0 ) */

  b[0] /= diag[0];                                     /* update b   */
  for (i = 1; i < n - 1; i++)
    b[i] = ( b[i] - b[i-1] * lower[i] ) / diag[i];

  for (tmp = ZERO, i = 0; i < n - 2; i++)
    tmp -= lowrow[i] * b[i];

  b[n-1] = ( b[n-1] + tmp - lower[n-1] * b[n-2] ) / diag[n-1];

  b[n-2] -= b[n-1] * upper[n-2];           /* back substitute         */
  for (i = n - 3; i >= 0; i--)
    b[i] -= upper[i] * b[i+1] + ricol[i] * b[n-1];

  return (0);
}

/* -------------------------- END fzdiag.c -------------------------- */
