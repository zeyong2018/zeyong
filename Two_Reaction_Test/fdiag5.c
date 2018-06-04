#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.12.1}{Systems with Five-Diagonal Matrices}
               {Systems with Five-Diagonal Matrices}*/

/*.BE*/
/* ------------------------- MODULE fdiag5.c ------------------------ */

#include <basis.h>
#include <u_proto.h>


#define MACH4_EPS (FOUR * MACH_EPS)

/*.BA*/

int diag5               /* 5 diagonal linear systems .................*/
/*.IX{diag5}*/
          (
           int   mod,             /* Modus: 0, 1, 2 ..................*/
           int   n,               /* size of matrix ..................*/
           REAL  ld2[],           /* 2. lower co-diagonal ............*/
           REAL  ld1[],           /* 1. lower co-diagonal ............*/
           REAL  d[],             /* main diagonal ...................*/
           REAL  ud1[],           /* 1. upper co-diagonal ............*/
           REAL  ud2[],           /* 2. upper co-diagonal ............*/
           REAL  b[]              /* right hand side/solution ........*/
          )
/*====================================================================*
 *                                                                    *
 *  diag5 solves the linear system  A * x = b  for a five diagonal    *
 *  n x n matrix A, which is given by the five vectors ld2, ld1, d,   *
 *  ud1, ud2 as follows:                                              *
 *                                                                    *
 *       ( d[0]     ud1[0]    ud2[0]   0     .  .   .        0      ) *
 *       ( ld1[1]   d[1]      ud1[1]   ud2[1]   0   .  .   . 0      ) *
 *       ( ld2[2]   ld1[2]    d[2]     ud1[2]   ud2[2]  0 .. 0      ) *
 *  A =  (   0      ld2[3]    ld1[3]   d[3]    .    .      . .      ) *
 *       (   .  .        .           .         .   .    .    0      ) *
 *       (   .    .         .            .        .    .   ud2[n-3] ) *
 *       (           .        .               .         .  ud1[n-2  ) *
 *       (   0 .   .   0       ld2[n-1]        ld1[n-1]    d[n-1]   ) *
 *                                                                    *
 *  A must be nonsingular.                                            *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *      Mainly used for diagonally dominant 5 diagonal matrices, as   *
 *      they occur in spline interpolation.                           *
 *      Diagonally dominant matrices always admit an unpivoted LU     *
 *      factorization. If A is not diagonally dominant, then the      *
 *      function band is recommended which uses column pivot search.  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mod      int mod;                                             *
 *               calling mode of diag5:                               *
 *       = 0     factor matrix and solve linear system                *
 *       = 1     factor only; factors are stored in ld2, ld1, d, ud1, *
 *               and ud2                                              *
 *       = 2     solve linear system only; for this call the LU       *
 *               decomposition must already be available in the above *
 *               five vectors. This is used to solve several systems  *
 *               for the same A and different right hand sides.       *
 *                                                                    *
 *   Input parameters:                                                *
 *   =================                                                *
 *      n        size of matrix  ( >= 3 )      int n                  *
 *                                                                    *
 *      ld2      2. lower co-diagonal          REAL   ld2[n]          *
 *      ld1      1. lower co-diagonal          REAL   ld1[n]          *
 *      d        main diagonal                 REAL   d  [n]          *
 *      ud1      1. upper co-diagonal          REAL   ud1[n]          *
 *      ud2      2. upper co-diagonal          REAL   ud2[n]          *
 *                                                                    *
 *               If  mod = 2, the vectors ld2, ld1, d, ud1, ud2       *
 *               contain the data of the LU factorization.            *
 *                                                                    *
 *      b        right hand side               REAL   b[n]            *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      ld1      )                                                    *
 *      ld2      ) if  mod = 0, 1, these contain the data for the     *
 *      d        ) factorization of A                                 *
 *      ud1      )                                                    *
 *      ud2      )                                                    *
 *                                                                    *
 *      b        solution vector               REAL   b[n]            *
 *                                                                    *
 *   The determinant of A is given by                                 *
 *      det A = d[0] * ... * d[n-1]                                   *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3                                                *
 *      = 2      Matrix is numerically singular                       *
 *      = 3      wrong call                                           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use:                                                *
 *   =================                                                *
 *                                                                    *
 *        int diag5dec (): forms LU factorization                     *
 *        int diag5sol (): solves the linear system                   *
 *                                                                    *
 *====================================================================*/
{
  int rc;

  switch (mod)
  {
    case 0:                                /* factor and solve .......*/
            rc = diag5dec (n, ld2, ld1, d, ud1, ud2);
            if (rc == 0)
              return (diag5sol (n, ld2, ld1, d, ud1, ud2, b));
            else
              return (rc);

    case 1:                                /* factor only ............*/
            return (diag5dec (n, ld2, ld1, d, ud1, ud2));

    case 2:                                /* solve only .............*/
            return (diag5sol (n, ld2, ld1, d, ud1, ud2, b));
  }

  return (3);
}


int diag5dec            /* LU factorization of a 5 diagonal matrix ...*/
/*.IX{diag5dec}*/
             (
              int   n,            /* size of matrix  .................*/
              REAL  ld2[],        /* 2. lower co-diagonal ............*/
              REAL  ld1[],        /* 1. lower co-diagonal ............*/
              REAL  d[],          /* main diagonal ...................*/
              REAL  ud1[],        /* 1. upper co-diagonal ............*/
              REAL  ud2[]         /* 2. upper co-diagonal ............*/
             )
/*====================================================================*
 *                                                                    *
 *  diag5dec factors a five diagonal n x n matrix A, which is given   *
 *  by the five vectors ld2, ld1, d, ud1, ud2 as follows:             *
 *                                                                    *
 *       ( d[0]     ud1[0]    ud2[0]   0     .  .   .        0      ) *
 *       ( ld1[1]   d[1]      ud1[1]   ud2[1]   0   .  .   . 0      ) *
 *       ( ld2[2]   ld1[2]    d[2]     ud1[2]   ud2[2]  0 .. 0      ) *
 *  A =  (   0      ld2[3]    ld1[3]   d[3]    .    .      . .      ) *
 *       (   .  .        .           .         .   .    .    0      ) *
 *       (   .    .         .            .        .    .   ud2[n-3] ) *
 *       (           .        .               .         .  ud1[n-2  ) *
 *       (   0 .   .   0       ld2[n-1]        ld1[n-1]    d[n-1]   ) *
 *                                                                    *
 *  A must be nonsingular.                                            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   =================                                                *
 *      n        size of matrix  ( >= 3 )      int n                  *
 *                                                                    *
 *      ld2      2. lower co-diagonal          REAL   ld2[n]          *
 *      ld1      1. lower co-diagonal          REAL   ld1[n]          *
 *      d        main diagonal                 REAL   d  [n]          *
 *      ud1      1. upper co-diagonal          REAL   ud1[n]          *
 *      ud2      2. upper co-diagonal          REAL   ud2[n]          *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      ld1      )                                                    *
 *      ld2      ) these contain the data for the LU                  *
 *      d        ) factorization of A                                 *
 *      ud1      )                                                    *
 *      ud2      )                                                    *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other wrong input                           *
 *      = 2      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*/
{
  register int i;
  REAL     row,      /* 1-norm of a row of matrix A                   */
           D,        /* reciprocal value of row (serves as reference  */
                     /* in comparing the diagonal elements with zero) */
           ud1i,     /* ud1[i] for i=2,...,n-2 and zero for i=n-1     */
           ud2i;     /* ud2[i] for i=2,...,n-3 and zero for i=n-2,n-1 */

  if (n < 3) return (1);     /* Dimension too small ..................*/

  if (ld2 == NULL || ld1 == NULL || d == NULL || ud1 == NULL ||
      ud2 == NULL) return (1);

  row = FABS(d[0]) + FABS(ud1[0]) + FABS(ud2[0]);
  if (row == ZERO)                        /* numerically  singular ...*/
    return 2;
  D = ONE / row;
  if (FABS(d[0]) * D <= MACH4_EPS)        /* numerically  singular ...*/
    return 2;

  ud1[0] /= d[0];
  ud2[0] /= d[0];
  row = FABS(ld1[1]) + FABS(d[1]) + FABS(ud1[1]) + FABS(ud2[1]);
  if (row == ZERO)                        /* numerically  singular ...*/
    return 2;
  D = ONE / row;

  d[1] -= ld1[1] * ud1[0];

  if (FABS(d[1]) * D <= MACH4_EPS)        /* numerically  singular ...*/
    return 2;

  ud1[1] = (ud1[1] - ld1[1] * ud2[0]) / d[1];
  ud2[1] /= d[1];

  ud1i = ud1[2];
  ud2i = ud2[2];
  for (i = 2; i < n; i++)
  {
    row = FABS(ld2[i]) + FABS(ld1[i]) + FABS(d[i]) +
          FABS(ud1i) + FABS(ud2i);
    if (row == ZERO)                      /* numerically  singular ...*/
      return 2;
    D = ONE / row;
    ld1[i] -= ld2[i] * ud1[i-2];
    d[i]   -= (ld2[i] * ud2[i-2] + ld1[i] * ud1[i-1]);

    if (FABS(d[i]) * D <= MACH4_EPS)      /* numerically  singular ...*/
      return 2;

    if (i < n - 1)
      ud1[i]  = (ud1[i] - ld1[i] * ud2[i-1]) / d[i];

    if (i < n - 2)
      ud2[i] /= d[i],
      ud1i   =  ud1[i + 1];
    else
      ud1i   =  ZERO;

    if (i < n - 3)
      ud2i = ud2[i + 1];
    else
      ud2i = ZERO;
  }

  return (0);
}


int diag5sol            /* solve a five diagonal linear system .......*/
/*.IX{diag5sol}*/
             (
              int   n,            /* size of matrix ..................*/
              REAL  ld2[],        /* 2. lower co-diagonal ............*/
              REAL  ld1[],        /* 1. lower co-diagonal ............*/
              REAL  d[],          /* main diagonal ...................*/
              REAL  ud1[],        /* 1. upper co-diagonal ............*/
              REAL  ud2[],        /* 2. upper co-diagonal ............*/
              REAL  b[]           /* right hand side/solution ........*/
             )
/*====================================================================*
 *                                                                    *
 *  diag5sol solves the linear system  L * U * x = b  for a given five*
 *  diagonal LU factorization given by the five vectors ld2, ld1, d,  *
 *  ud1, ud2 in diag5dec                                              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   =================                                                *
 *      n        size of matrix  ( >= 3 )      int n                  *
 *                                                                    *
 *      ld2      2. lower co-diagonal          REAL   ld2[n]          *
 *      ld1      1. lower co-diagonal          REAL   ld1[n]          *
 *      d        main diagonal                 REAL   d  [n]          *
 *      ud1      1. upper co-diagonal          REAL   ud1[n]          *
 *      ud2      2. upper co-diagonal          REAL   ud2[n]          *
 *                                                                    *
 *      The vectors ld2, ld1, d, ud1, ud2 contain the data of the     *
 *      LU factorization.                                             *
 *                                                                    *
 *      b        right hand side               REAL   b[n]            *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      b        solution vector               REAL   b[n]            *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3                                                *
 *      = 2      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*/
{
  register int i;

  if (n < 3) return (1);

  if (ld2 == NULL || ld1 == NULL || d == NULL || ud1 == NULL ||
      ud2 == NULL || b == NULL) return (1);

  if (ABS(d[0]) < MACH_EPS) return (2);

                                                /* update b ..........*/
  b[0] /= d[0];

  if (ABS(d[1]) < MACH_EPS) return (2);
  b[1] = (b[1] - ld1[1] * b[0]) / d[1];

  for (i = 2; i < n; i++)
  {
    if (ABS (d[i]) < MACH_EPS) return (2);
    b[i] = (b[i] - ld2[i] * b[i-2] - ld1[i] * b[i-1]) / d[i];
  }

                                             /* back substitution ....*/
  b[n-2] -= ud1[n-2] * b[n-1];

  for (i = n - 3; i >= 0; i--)
    b[i] -= (ud1[i] * b[i+1] + ud2[i] * b[i+2]);

  return (0);
}

/* --------------------------- END fdiag5.c ------------------------- */
