#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.12.2}
     {Systems with Five-Diagonal Symmetric Matrices}
     {Systems with Five-Diagonal Symmetric Matrices}*/

/*.BE*/
/* ----------------------- MODULE fdiag5pd.c ------------------------ */

#include <basis.h>
#include <u_proto.h>

/*.BA*/

int diag5pd             /* 5 diagonal symmetric strongly nonsingular .*/
/*.IX{diag5pd}*/
            (
             int   mod,           /* Modus: 0, 1, 2 ..................*/
             int   n,             /* # matrix rows ...................*/
             REAL  d[],           /* main diagonal ...................*/
             REAL  ud1[],         /* first co-diagonal ...............*/
             REAL  ud2[],         /* second co-diagonal ..............*/
             REAL  b[]            /* Right hand side .................*/
            )
/*====================================================================*
 *                                                                    *
 *  diag5pd determine the solution  x  of the linear system           *
 *  A * x = b with a  5-diagonal strongly nonsingular symmetric       *
 *  n x n coefficient matrix A, given by the 3 vectors d, ud1, ud2:   *
 *                                                                    *
 *       ( d[0]     ud1[0]    ud2[0]   0     .  .   .        0      ) *
 *       ( ud1[0]   d[1]      ud1[1]   ud2[1]   0   .  .   . 0      ) *
 *       ( ud2[0]   ud1[1]    d[2]     ud1[2]   ud2[2]  0 .. 0      ) *
 *  A =  (   0      ud2[1]    ud1[2]   d[3]    .    .      . .      ) *
 *       (   .  .        .           .       .     .    .    0      ) *
 *       (   .    .         .            .        .     .  ud2[n-3] ) *
 *       (           .        .               .         .  ud1[n-2] ) *
 *       (   0 .   .   0       ud2[n-3]        ud1[n-2]    d[n-1]   ) *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *      Such matrices occur during spline interpolation.              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mod      int mod;                                             *
 *               kind of call of diag5pd;                             *
 *       = 0     factor matrix and solve linear system                *
 *       = 1     factor only, write answer onto  d, ud1, and ud2      *
 *       = 2     Solve only. Here the factorization must be available *
 *               This modus is useful to solve several sysstems of    *
 *               equations for the same system matrix.                *
 *                                                                    *
 *   Input parameters:                                                *
 *   =================                                                *
 *      n        Dimension of matrix  ( >= 3 ) int n                  *
 *                                                                    *
 *      d        main diagonal                 REAL   d  [n]          *
 *      ud1      1. co-diagonal                REAL   ud1[n]          *
 *      ud1      2. co-diagonal                REAL   ud2[n]          *
 *                                                                    *
 *               If  mod = 2 ,  d, ud1, ud2  contain the R'*D*R       *
 *               factorization                                        *
 *                                                                    *
 *      b        right hand side of system     REAL   b[n]            *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      d        ) If  mod = 0, 1  these contain the factors of the   *
 *      ud1      ) matrix in place of the original matrix in condensed*
 *      ud2      ) form.                                              *
 *                                                                    *
 *      b        solution vector               REAL   b[n]            *
 *                                                                    *
 *   The determinant of the matrix is given ad if rep = 0 or 3        *
 *      det A = d[0] * ... * d[n-1]                                   *
 *                                                                    *
 *   Return value:                                                    *
 *   ============                                                     *
 *      = 0      all ok                                               *
 *      = 1      n < 3                                                *
 *      = 2      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   functions in use  :                                              *
 *   ===================                                              *
 *                                                                    *
 *        int diag5pddec (): Computes R'*D*R factorization            *
 *        int diag5pdsol (): solves linear system                     *
 *                                                                    *
 *====================================================================*/
{
  int rc;

  if (n < 3) return (1);  /*  Dimension  wrong .......................*/

  switch (mod)
  {
    case 0:                      /* Decompose and solve system .......*/
            rc = diag5pddec (n, d, ud1, ud2);
            if (rc == 0)
              return (diag5pdsol (n, d, ud1, ud2, b));
            else
              return (rc);

    case 1:             /* Factor only ...............................*/
            return (diag5pddec (n, d, ud1, ud2));

    case 2:            /* Solve only .................................*/
            return (diag5pdsol (n, d, ud1, ud2, b));
  }

  return (3);
}


int diag5pddec       /* Factor 5 diagonal strongly nonsingular matrix */
/*.IX{diag5pddec}*/
               (
                int   n,          /* # Matrix rows ...................*/
                REAL  d[],        /* main diagonal ...................*/
                REAL  ud1[],      /* 1. co-diagonal ..................*/
                REAL  ud2[]       /* 2. co-diagonal ..................*/
               )
/*====================================================================*
 *                                                                    *
 *  diag5pddec decomposes the  5-diagonal  symmetric and strongly     *
 *  nonsingular matrix A, which is given by the 3 vectors d, ud1, ud2 *
 *                                                                    *
 *       ( d[0]     ud1[0]    ud2[0]   0     .  .   .        0      ) *
 *       ( ud1[0]   d[1]      ud1[1]   ud2[1]   0   .  .   . 0      ) *
 *       ( ud2[0]   ud1[1]    d[2]     ud1[2]   ud2[2]  0 .. 0      ) *
 *  A =  (   0      ud2[1]    ud1[2]   d[3]    .    .      . .      ) *
 *       (   .  .        .           .       .     .    .    0      ) *
 *       (   .    .         .            .        .     .  ud2[n-3] ) *
 *       (           .        .               .         .  ud1[n-2] ) *
 *       (   0 .   .   0       ud2[n-3]        ud1[n-2]    d[n-1]   ) *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        Dimension of A   ( >= 3 )     int n                  *
 *                                                                    *
 *      d        main diagonal                 REAL   d  [n]          *
 *      ud1      1. co-diagonal                REAL   ud1[n]          *
 *      ud2      2. co-diagonal                REAL   ud2[n]          *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      d        ) contain the Cholesky factors   R' * D * R          *
 *      ud1      )                                                    *
 *      ud2      )                                                    *
 *                                                                    *
 *   Retrurn values :                                                 *
 *   ================                                                 *
 *      = 0      all ok                                               *
 *      = 1      n < 3                                                *
 *      = 2      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*/
{
  register int i;
  REAL e_1, e_2, tmp, sum;

  if (n < 3) return (1);                       /* n too small ?       */

  if (d == NULL || ud1 == NULL || ud2 == NULL) return (1);

  ud1[n-1] = ZERO;
  ud2[n-1] = ZERO;
  ud2[n-2] = ZERO;

  sum = ABS(d[0]) + ABS(ud1[0]) + ABS(ud2[0]);
  if (sum == ZERO) return (3);

  if (d[0] < sum * MACH_EPS) return (2);

  tmp     = ud1[0];
  ud1[0] /= d[0];
  e_2     = ud2[0];
  ud2[0] /= d[0];

  sum = ABS(tmp) + ABS(d[1]) + ABS(ud1[1]) + ABS(ud2[1]);
  if (sum == ZERO) return (2);            /* singular                 */

  d[1] -= tmp * ud1[0];

  if (ABS(d[1]) < sum * MACH_EPS) return (2); /* not str. nonsingular */

  tmp     = ud1[1];
  ud1[1]  = (ud1[1] - e_2 * ud1[0]) / d[1];
  e_1     = ud2[1];
  ud2[1] /= d[1];

  for (i = 2; i < n; i++)
  {
    sum = ABS(e_2) + ABS(tmp) + ABS(d[i]) + ABS(ud1[i]) + ABS(ud2[i]);
    if (sum == ZERO) return (2);

    d[i]   -= e_2 * ud2[i-2] + d[i-1] * SQR (ud1[i-1]);

    if (ABS(d[i]) < sum * MACH_EPS) return (2); /* not str. nonsing.? */

    if (i < n - 1)
    {
      tmp = ud1[i];
      ud1[i]  = (ud1[i] - e_1 * ud1[i-1]) / d[i];
    }

    if (i < n - 2)
    {
      e_2     = e_1;
      e_1     = ud2[i];
      ud2[i] /= d[i];
    }
  }

  return (0);
}


int diag5pdsol          /* Solve systems for 5 diagonal symmetric m. .*/
/*.IX{diag5pdsol}*/
               (
                int   n,          /* size of matrix ..................*/
                REAL  d[],        /* main diagonal ...................*/
                REAL  ud1[],      /* 1. co-diagonal ..................*/
                REAL  ud2[],      /* 2. co-diagonal ..................*/
                REAL  b[]         /* Right hand side .................*/
               )
/*====================================================================*
 *                                                                    *
 *  diag5pdsol solves the linear system                               *
 *  R'  * D * R * x = b, where the Cholesky decomposition is stored   *
 *  in the 3 vectors d, ud1 and ud2 from diag5pddec                   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        Dimension of matrix  ( >= 3 ) int n                  *
 *                                                                    *
 *      d        main diagonal                 REAL   d  [n]          *
 *      ud1      1. co-diagonal                REAL   ud1[n]          *
 *      ud2      2. co-diagonal                REAL   ud2[n]          *
 *                                                                    *
 *      containing the Cholesky factors of a 5-diagonal matrix        *
 *                                                                    *
 *      b        right hand side of system                            *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      b        solution vector               REAL   b[n]            *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 3                                                *
 *      = 2      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*/
{
  register int i;

  if (n < 3) return (1);

  if (d == NULL || ud1 == NULL || ud2 == NULL || b == NULL) return (1);

  /* updating ........................................................*/

  if (ABS(d[1]) < MACH_EPS) return (2);
  b[1] -= ud1[0] * b[0];

  for (i = 2; i < n; i++)
    b[i] -= ud1[i-1] * b[i-1] + ud2[i-2] * b[i-2];

  for (i = 0; i < n; i++)
  {
    if (ABS (d[i]) < MACH_EPS) return (2);
    b[i] /= d[i];
  }

  /* back substitution ...............................................*/

  b[n-2] -= ud1[n-2] * b[n-1];

  for (i = n - 3; i >= 0; i--)
    b[i] -= (ud1[i] * b[i+1] + ud2[i] * b[i+2]);

  return (0);
}

/* -------------------------- END fdiag5pd.c ------------------------ */
