#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.13}{Linear Systems with Band Matrices}
             {Linear Systems with Band Matrices}*/

/*.FE{}
     {Gau"s Algorithm for Band Matrices using Pivots}
     {Gau"s Algorithm for Band Matrices using Pivots}*/

/*.BE*/
/* -------------------------- MODULE fband.c ------------------------ */

#include <basis.h>
#include <u_proto.h>
/*.BA*/

int band                /* Linear systems with banded matrices .......*/
/*.IX{band}*/
         (
          int    mod,             /* Modus: 0, 1, 2 ..................*/
          int    n,               /* size of system ..................*/
          int    ld,              /* # of lower co-diagonals .........*/
          int    ud,              /* # of upper co-diagonals .........*/
          REAL * pmat[],          /* condensed input matrix ..........*/
          REAL   b[],             /* right hand side .................*/
          int    perm[],          /* row permutation vector ..........*/
          int *  signd            /* sign of perm ....................*/
         )
/*====================================================================*
 *                                                                    *
 *  The function band solves a linear banded system:  pmat * x = b.   *
 *  Here pmat is a nonsingular n x n matrix in condensed form, i.e.   *
 *  represented in an ld+ud+1 x n matrix for its ld lower and ud upper*
 *  co-diagonals. b denotes the right hand side of the system, and x  *
 *  is the solution.                                                  *
 *                                                                    *
 *  band uses the Gauss algorithm with column pivot search.           *
 *  The result of pivoting are min( ud, ld) additional columns,       *
 *  so that pmat needs all in all a n x (ld+1+ud+min(ld,ud)) matrix.  *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Applications:                                                    *
 *   =============                                                    *
 *      Solve linear systems with nonsingular banded system matrices. *
 *      Particularly useful for large sparse and banded matrices with *
 *      n >> ld+1+ud.                                                 *
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
 *      perm     int perm[n];          ( for mod = 2 )                *
 *               row permutation vector                               *
 *      signd    int *signd;           ( for mod = 2 )                *
 *               sign of perm. The determinant of A can be computed   *
 *               as the product of the diagonal entries times signd.  *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      pmat     REAL   *pmat[n];      ( for mod = 0, 1 )             *
 *               LU factorization in condensed form                   *
 *      perm     int perm[n];          ( for mod = 0, 1 )             *
 *               row permutation vector                               *
 *      b        REAL   b[n];          ( for mod = 0, 2 )             *
 *               solution vector for the system                       *
 *      signd  int *signd;             ( for mod = 0, 1 )             *
 *               sign of perm                                         *
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
 *      int banddec() : Factor matrix                                 *
 *      int bandsol() : solve linear system                           *
 *                                                                    *
 *====================================================================*/
{
  int rc;

  if (n < 1 || ld < 0 || ud < 0) return(1);

  switch (mod)
  {
    case 0: /* factor and solve system ...............................*/
            rc = banddec (n, ld, ud, pmat, perm, signd);
            if (rc == 0)
              return (bandsol (n, ld, ud, pmat, b, perm));
            else
              return (rc);

    case 1: /* factor only ...........................................*/
            return (banddec (n, ld, ud, pmat, perm, signd));

    case 2: /* solve only ............................................*/
            return (bandsol (n, ld, ud, pmat, b, perm));
  }

  return (4);
}


int banddec             /* Factor a banded matrix ....................*/
/*.IX{banddec}*/
            (
             int    n,            /* size of system ..................*/
             int    ld,           /* # of lower co-diagonals .........*/
             int    ud,           /* # of upper co-diagonals .........*/
             REAL * pmat[],       /* condensed input matrix ..........*/
             int    perm[],       /* row permutation vector ..........*/
             int *  signd         /* sign of perm ....................*/
            )
/*====================================================================*
 *                                                                    *
 *  The function banddec factors a condensed banded matrix pmat.      *
 *  banddec uses the Gauss algorithm with column pivot search.        *
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
 *      perm     int perm[n];                                         *
 *               row permutation vector                               *
 *      signd    int *signd;                                          *
 *               sign of perm. The determinant of A can be computed   *
 *               as the product of the diagonal entries times signd.  *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      pmat     REAL   *pmat[n];                                     *
 *               LU factorization in condensed form                   *
 *      perm     int perm[n];                                         *
 *               row permutation vector                               *
 *      signd  int *signd;                                            *
 *               sign of perm                                         *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *      = 2      lack of space                                        *
 *      = 3      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :     NULL, MACH_EPS                              *
 *   ================                                                 *
 *                                                                    *
 *   Macros: min, max, SWAP                                           *
 *   =======                                                          *
 *                                                                    *
 *====================================================================*/
{
  int      j0, mm, up, istart, iend, step, kstart,
           kend, kjend, km, jm, jk;

  register k, j, i;
  REAL   piv;

  if (ld < 0 || ud < 0 || n < 1)
   return (1);                         /* Invalid input parameter    */

 if (pmat == NULL) return (1);

 for (i = 0; i < n; i++)
   if (pmat[i] == NULL) return (1);
  mm = ld + 1 + ud + min (ld, ud);

  up = ld <= ud;                       /* up = 0 ==> transform into   */
                                       /* lower triangular matrix     */
  for (i = 0; i < n; i++)
    for (k = ld + ud + 1; k < mm; k++) /* initialize needed extra     */
      pmat[i][k] = ZERO;               /* columns                     */

  *signd = 1;                               /* initialize signd       */
  if (up)
  {
    istart = 0; iend = n-1; step = 1;       /* Find start, end and    */
    kstart = 1;                             /* directions of loops    */
  }                                         /* depending of up        */
  else
  {
    istart = n-1; iend = 0; step = -1;
    kstart = -1;
  }

  for (i = istart; i != iend; i += step)    /* loop over all rows     */
  {
    kend = (up ? min (ld+1, n-i) : max (-i-1, -ud-1));
    j0 = 0;
    piv = ABS (pmat[i][ld]);                /* choose pivot           */
    for (k = kstart; k != kend; k += step)
      if ( ABS (pmat[k+i][ld-k]) > piv )
      {
        piv = ABS (pmat[k+i][ld-k]);
        j0 = k;
      }

    if (piv < MACH_EPS) return (3);        /* If piv = 0, matrix is   */
                                           /* singular                */
    perm[i] = j0;
    kjend =  (up ? min (j0+ud+1, n-i) : max ( -i-1, j0-ld-1));

    if (j0 != 0)
    {
      *signd = - *signd;                   /* swap rows               */
      for (k = 0; k != kjend; k += step)
      {
        km = k + ld;
        if (km < 0) km += mm;
        SWAP (REAL, pmat[i][km], pmat[i+j0][k+ld-j0])
      }
    }

    for (k = kstart; k != kend; k +=step)      /* loop over all rows  */
    {                                          /* below row i         */
      pmat[k+i][ld-k] /= pmat[i][ld];
      for (j = kstart; j != kjend; j += step)
      {                                        /* loop over all       */
        jk = j + ld - k;                       /* columns right of i  */
        jm = j + ld;
                                  /* additional columns from pivoting */
        if (jk < 0) jk += mm;     /* are stored to right of column    */
        if (jm < 0) jm += mm;     /*  ud+ld+1                         */

        pmat[k+i][jk] -= pmat[k+i][ld-k] * pmat[i][jm];
      }
    }   /* end k */
  }   /* end i */

  piv = ABS (pmat[iend][ld]);              /* choose pivot            */
  if (piv < MACH_EPS) return(3);           /* If piv = 0, matrix is   */
                                           /* singular                */
  perm[iend] = 0;
  return (0);
}


int bandsol             /* Solve a banded system .....................*/
/*.IX{bandsol}*/
            (
             int    n,            /* size of system ..................*/
             int    ld,           /* # of lower co-diagonals .........*/
             int    ud,           /* # of upper co-diagonals .........*/
             REAL * pmat[],       /* condensed input matrix ..........*/
             REAL   b[],          /* right hand side .................*/
             int    perm[]        /* row permutation vector ..........*/
            )
/*====================================================================*
 *                                                                    *
 *  The function bandsol solves a factored linear banded system in    *
 *  condensed form using banddec.                                     *
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
 *      perm     int perm[n];                                         *
 *               row permutation vector                               *
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
 *   Macros: min, max, SWAP                                           *
 *   =======                                                          *
 *                                                                    *
 *====================================================================*/
{
  register i, k;
  int      s, mm, up, istart, iend, step, kstart,
           kend, km;
                                       /* invalid input               */
 if (ld < 0 || ud < 0 || n < 1)
    return (1);

  if (pmat == NULL || b == NULL || perm == NULL) return (1);

  for (i = 0; i < n; i++)
    if (pmat[i] == NULL) return (1);

  mm = ld + ud + 1 + min (ld, ud);     /* mm = max. column number     */

  up = ld <= ud;                       /* up = 0 ==> Tranform to      */
                                       /* lower triangular matrix     */
  if (up)
  {
    istart = 0; iend = n-1; step = 1;    /* determine bounds and      */
    kstart = 1; s = -1;                  /* direction of loop         */
  }                                      /* depending on up           */
  else
  {
    istart = n-1; iend = 0; step = -1;
    kstart = -1; s = 1;
  }

  for (i = istart; i != iend; i += step)    /* update right hand side */
  {                                         /* according to perm      */
    if (perm[i] != 0)
      SWAP (REAL, b[i], b[i+perm[i]])

    kend = (up ? min (ld+1, n-i) : max (-i-1, -ud-1));

    for (k = kstart; k != kend; k += step)
      b[k+i] -= pmat[k+i][ld-k] * b[i];
  }

  for (i = iend; i != istart + s ; i -= step)
  {
    kend =  (up ? min (ld+ud+1, n-i) : max (-i-1, -ud-ld));

    for (k = kstart; k != kend; k += step)
    {
      km = k + ld;                     /* update and                  */
      if (km < 0) km += mm;            /* back substitute             */
        b[i] -= pmat[i][km] * b[i+k];
    }
    b[i] /= pmat[i][ld];
  }

  return (0);
}

/* --------------------------- END fband.c -------------------------- */
