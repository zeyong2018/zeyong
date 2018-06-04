#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{}
     {Gau"s Algorithm for Band Matrices without Pivots}
     {Gau"s Algorithm for Band Matrices without Pivots}*/

/*.BE*/
/* ------------------------- MODULE fbando.c ------------------------ */

#include <basis.h>
#include <u_proto.h>

/*.BA*/

int bando               /* Linear banded system without using pivots .*/
/*.IX{bando}*/
          (
          int    mod,             /* Modus: 0, 1, 2 ..................*/
          int    n,               /* size of system ..................*/
          int    ld,              /* # of lower co-diagonals .........*/
          int    ud,              /* # of upper co-diagonals .........*/
          REAL * pmat[],          /* condensed input matrix ..........*/
          REAL   b[]              /* right hand side .................*/
         )
/*====================================================================*
 *                                                                    *
 *  The function bando solves a linear banded system:  pmat * x = b.  *
 *  Here pmat is a nonsingular n x n matrix in condensed form, i.e.   *
 *  represented in an ld+ud+1 x n matrix for its ld lower and ud upper*
 *  co-diagonals. b denotes the right hand side of the system, and x  *
 *  is the solution.                                                  *
 *                                                                    *
 *  bando uses the Gauss algorithm without column pivot search.       *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   =============                                                    *
 *      Solve linear systems with nonsingular banded system matrices. *
 *      Particularly useful for large sparse and banded and diagonally*
 *      dominant matrices.                                            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mod      int mod;                                             *
 *               calling modus for band:                              *
 *       = 0     factor matrix and solve linear system                *
 *       = 1     factor only, store factors in pmat                   *
 *       = 2     solve linear system only; for this call, the factors *
 *               must already be available in pmat, such as when      *
 *               many systems are solved for differing right hand     *
 *               sides and the same system matrix.                    *
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 2 )                                    *
 *               Dimension of pmat, size of b                         *
 *      ld       int ld; ( ld >= 0 )                                  *
 *               number of lower co-diagonals                         *
 *      ud       int ud; ( ud >= 0 )                                  *
 *               number of upper co-diagonals                         *
 *      pmat     REAL   *pmat[n];                                     *
 *               mod = 0, 1:                                          *
 *               Matrix of the system in comndensed form.             *
 *               Each row has length at least ld + 1 + ud + min(ld,ud)*
 *               where the columns 0, .., ld-1 denote the lower       *
 *               co-diagonals, column ld stores the diagonal and the  *
 *               columns ld+1, .., ld+ud contain the upper            *
 *               co-diagonals.                                        *
 *               If A is the original uncondensed band matrix, then : *
 *               A[i][k] = pmat[i][ld+k-i],                           *
 *                           for k,i inside the band                  *
 *               mod = 2:                                             *
 *               LU factors in condensed form.                        *
 *      b        REAL   b[n];          ( for mod = 0, 2 )             *
 *               Right hand side                                      *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      pmat     REAL   *pmat[n];      ( for mod = 0, 1 )             *
 *               LU factorization in condensed form                   *
 *      b        REAL   b[n];          ( for mod = 0, 2 )             *
 *               solution vector for the system                       *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *      = 2      lack of space                                        *
 *      = 3      Matrix is numerically singular                       *
 *      = 4      wrong calling modus                                  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used :                                                 *
 *   ================                                                 *
 *      int banodec() : Factor matrix                                 *
 *      int banosol() : solve linear system                           *
 *                                                                    *
 *====================================================================*/
{
  int rc;

  if (n < 3 || ld < 0 || ud < 0 || n < ld + 1 + ud) return (1);

  switch (mod)
  {
    case 0: /* factor matrix and solve system ........................*/
            rc = banodec (n, ld, ud, pmat);
            if (rc == 0)
              return (banosol (n, ld, ud, pmat, b));
            else
              return (rc);

    case 1: /* factor only ...........................................*/
            return (banodec (n, ld, ud, pmat));

    case 2: /* solve only ............................................*/
            return (banosol (n, ld, ud, pmat, b));
  }

  return (3);
}


int banodec             /* Decompose a banded matrix .................*/
/*.IX{banodec}*/
            (
             int    n,            /* size of system ..................*/
             int    ld,           /* # of lower co-diagonals .........*/
             int    ud,           /* # of upper co-diagonals .........*/
             REAL * pmat[]        /* condensed input matrix ..........*/
            )
/*====================================================================*
 *                                                                    *
 *  The function banodec factors a condensed banded matrix pmat.      *
 *  banddec uses the Gauss algorithm without column pivot search.     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 2 )                                    *
 *               Dimension of pmat, size of b                         *
 *      ld       int ld; ( ld >= 0 )                                  *
 *               number of lower co-diagonals                         *
 *      ud       int ud; ( ud >= 0 )                                  *
 *               number of upper co-diagonals                         *
 *      pmat     REAL   *pmat[n];                                     *
 *               Matrix of the system in comndensed form.             *
 *               Each row has length at least ld + 1 + ud + min(ld,ud)*
 *               where the columns 0, .., ld-1 denote the lower       *
 *               co-diagonals, column ld stores the diagonal and the  *
 *               columns ld+1, .., ld+ud contain the upper            *
 *               co-diagonals.                                        *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      pmat     REAL   *pmat[n];                                     *
 *               LU factorization in condensed form                   *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *      = 2      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :     NULL, MACH_EPS                              *
 *   ================                                                 *
 *                                                                    *
 *   Macros: min                                                      *
 *   =======                                                          *
 *                                                                    *
 *====================================================================*/
{
  int      kend, kjend, jm, jk;
  register k, j, i;

  if (ld == 0) return (0);           /* Matrix already in upper       */
                                     /* triangular form               */
  if (pmat == NULL) return (1);

  for (k = 0; k < n; k++)
    if (pmat[k] == NULL) return (1);

  for (i = 0; i < n - 1; i++)               /* loop over all rws      */
  {
    kend =  min (ld + 1, n - i);
    kjend = min (ud + 1, n - i);

    if (ABS(pmat[i][ld]) < MACH_EPS)        /* LU decompsition does   */
      return (2);                           /* not exist              */

    for (k = 1; k != kend; k++)             /* loop over all rows     */
    {                                       /* below row i            */
      pmat[k+i][ld-k] /= pmat[i][ld];

      for (j = 1; j != kjend; j++)
      {
        jk = j + ld - k;
        jm = j + ld;
        pmat[k+i][jk] -= pmat[k+i][ld-k] * pmat[i][jm];
      }
    } /* end k */
  }  /* end i */

  return (0);
}


int banosol             /* Solve a banded system .....................*/
/*.IX{banosol}*/
            (
             int    n,            /* size of system ..................*/
             int    ld,           /* # of lower co-diagonals .........*/
             int    ud,           /* # of upper co-diagonals .........*/
             REAL * pmat[],       /* condensed input matrix ..........*/
             REAL   b[]           /* right hand side .................*/
            )
/*====================================================================*
 *                                                                    *
 *  The function banosol solves a factored linear banded system in    *
 *  condensed form using banodec.                                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 2 )                                    *
 *               Dimension of pmat, size of b                         *
 *      ld       int ld; ( ld >= 0 )                                  *
 *               number of lower co-diagonals                         *
 *      ud       int ud; ( ud >= 0 )                                  *
 *               number of upper co-diagonals                         *
 *      pmat     REAL   *pmat[n];                                     *
 *               Matrices for the factored system in comndensed form. *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      b        REAL   b[n];                                         *
 *               solution vector of linear system                     *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :     NULL                                        *
 *   ================                                                 *
 *                                                                    *
 *   Macros: min                                                      *
 *   =======                                                          *
 *                                                                    *
 *====================================================================*/
{
  int      kend;
  register i, k;
                                        /*  Invalid input parameter   */
  if (n < 3 || ld < 0 || ud < 0 || ud + ld + 1 > n) return (1);

  if (pmat == NULL || b == NULL) return (1);

  for (i = 0; i < n; i++)
    if (pmat[i] == NULL) return (1);

  for (i = 0; i < n - 1; i++)
  {
    kend = min (ld + 1, n - i);
    for (k = 1; k != kend; k++)
      b[k+i] -= pmat[k+i][ld-k] * b[i];
  }

  for (i = n - 1; i >= 0 ; i--)          /* back substitution        */
  {
    kend = min (ud + 1, n - i);
    for (k = 1; k < kend; k++)
      b[i] -= pmat[i][k+ld] * b[i+k];
    b[i] /= pmat[i][ld];
  }

  return (0);
}

/* -------------------------- END fbando.c -------------------------- */
