#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE fzytrdsy.c ------------------------ */

/***********************************************************************
*                                                                      *
* Solve a linear system with a cyclic tridiagonal symmetric and        *
* -------------------------------------------------------------        *
* strongly nonsingular matrix                                          *
* ---------------------------                                          *
*                                                                      *
* exported function:                                                   *
*   - zytrdsy():  Solution                                             *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Gisela Engeln-Muellges (FORTRAN)               *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               FORTRAN code                                   *
* Datum:                1.12.1993                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>         /*  for  REAL, FOUR, MACH_EPS, FABS, ZERO, */
                           /*       ONE, sqr                          */
#include <fzytrdsy.h>      /*  for  zytrdsy                           */



/* ------------------------------------------------------------------ */

static int zzytrdsy
/*.IX{zzytrdsy}*/
                   (
                    int  n,
                    REAL diag[],
                    REAL oben[],
                    REAL rechts[]
                   )


/***********************************************************************
* The cyclic tridiagonal symmetric ans strongly nonsingular matrix A,  *
* given by its diagonal diag and co-diagonal oben, (see the function   *
* fzytrdsy()), is decomposed into the product  R' * D * R  using the   *
* root free Cholesky method.                                           *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n       number of equations (n > 2)                                  *
* diag    [0..n-1] vector, the main diagonal                           *
* oben    [0..n-1] vector, with the co-diagonal in positions 0,..,n-2  *
*         and A[0][n-1] in oben[n-1].                                  *
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* diag  \  Cholesky decomposition of A in condensed form:              *
* oben   > The first co-diagonal of R is stored in oben [0,...,n-2],   *
* rechts/  its right most column is stored in rechts[0,...,n-3] (its   *
*          diagonal is made up of ones). the diagonal matrix D is      *
*          stored in diag[0..n-1].                                     *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: no error                                                        *
* = 2: matrix A numerically not strongly nonsingular                   *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, FOUR, MACH_EPS, FABS, ZERO, ONE, sqr                           *
***********************************************************************/

{
  REAL eps,      /* relative error bound                              */
       zeile,    /* one norm of a row of A                            */
       d,        /* reciprocal of zeile                               */
       unten,    /* lower co-diagonal element in row i                */
       summe;    /* aux variable to compute diag[n-1]                 */
  int  i;        /* aux variable: Number of a row                     */


  eps = FOUR * MACH_EPS;

  zeile = FABS(diag[0]) + FABS(oben[0]) +  /* check first row of A    */
          FABS(oben[n - 1]);               /* for strong              */
  if (zeile == ZERO)                       /* nonsingularity          */
    return 2;
  d = ONE / zeile;
  if (FABS(diag[0]) * d <= eps)
    return 2;


  unten     =  oben[0];                                  /* factorize */
  oben[0]   /= diag[0];
  rechts[0] =  oben[n - 1] / diag[0];

  for (i = 1; i < n - 1; i++)
  {
    zeile   =  FABS(diag[i]) + FABS(oben[i]) + FABS(unten);
    diag[i] -= unten * oben[i - 1];

    if (zeile == ZERO)              /* check strong nonsingularity    */
      return 2;
    d = ONE / zeile;
    if (FABS(diag[i]) * d <= eps)
      return 2;

    if (i < n - 2)
      rechts[i] =  -unten * rechts[i - 1] / diag[i],
      unten     =  oben[i],
      oben[i]   /= diag[i];
    else
      oben[i] = (oben[i] - unten * rechts[i - 1]) / diag[i];
  }

  zeile       =  FABS(oben[n - 1]) + FABS(diag[n - 1]) +
                 FABS(oben[n - 2]);
  diag[n - 1] -= diag[n - 2] * sqr(oben[n - 2]);
  for (summe = ZERO, i = 0; i < n - 2; i++)
    summe += diag[i] * sqr(rechts[i]);
  diag[n - 1] -= summe;

  if (zeile == ZERO)
    return 2;
  d = ONE / zeile;
  if (FABS(diag[n - 1]) * d <= eps)
    return 2;


  return 0;
}



/* ------------------------------------------------------------------ */

static int lzytrdsy
/*.IX{lzytrdsy}*/
                   (
                    int  n,
                    REAL diag[],
                    REAL oben[],
                    REAL rechts[],
                    REAL rs[]
                   )

/***********************************************************************
* Solve the linear system                                              *
*                        A * X = RS                                    *
* for a cyclic tridiagonal symmetric and strongly nonsingular matrix   *
* A  from its Cholesky factorization   A  =  R' * D * R  as given by   *
* zzytrdsy().                                                          *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n       number of equations  (n > 2)                                 *
* diag    [0..n-1] vector, the diagonal of D                           *
* oben    [0..n-2] vector, the first co-diagonal of R                  *
* rechts  [0..n-3] vector, the right most column of R                  *
* rs      [0..n-1] vector, the right hand side                         *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* rs    [0..n-1] vector, the solution                                  *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: no error                                                        *
* = 4: invalid decomposition, D is singular                            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO                                                           *
***********************************************************************/

{
  REAL rsi1,     /* auxiliary values during updating                  */
       summe,    /* Sum, needed at end of updating                    */
       rsn1;     /* aux values for back substitution                  */
  int  i;        /* Loop variable                                     */


  if (diag[0] == ZERO)
    return 4;
  rsi1  = rs[0];                                            /* update */
  rs[0] = rsi1 / diag[0];
  summe = rechts[0] * rsi1;

  for (i = 1; i < n - 1; i++)
  {
    if (diag[i] == ZERO)
      return 4;
    rsi1  = rs[i] - oben[i - 1] * rsi1;
    rs[i] = rsi1 / diag[i];
    if (i < n - 2)
      summe += rechts[i] * rsi1;
  }

  if (diag[n - 1] == ZERO)
    return 4;
  rsi1      = rs[n - 1] - oben[n - 2] * rsi1 - summe;
  rs[n - 1] = rsi1 / diag[n - 1];


  rsn1 = rs[n - 1];                               /* back substitute  */
  rsi1 = rs[n - 2] -= oben[n - 2] * rsn1;
  for (i = n - 3; i >= 0; i--)
    rsi1 = rs[i] -= oben[i] * rsi1 + rechts[i] * rsn1;


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int zytrdsy    /* cyclic tridiagonal symmetric linear system .........*/
/*.IX{zytrdsy}*/
           (
            int  modus,     /* Modus of call: 0, 1, 2 ................*/
            int  n,         /* size of matrix ........................*/
            REAL diag[],    /* main diagonal or D in R'*D*R ..........*/
            REAL oben[],    /* first codiagonal in R in R'*D*R .......*/
            REAL rechts[],  /* right most column of R in R'*D*R ......*/
            REAL rs[]       /* right hand side/solution ..............*/
           )                /* error code ............................*/

/***********************************************************************
* Solve a linear system                                                *
*                        A * X = RS                                    *
* for a cyclic tridiagonal symmetric and strongly nonsingular matrix   *
* A. The matrix A is given by the two vectors oben and diag:           *
*                                                                      *
* diag[0]  *x[0]   + oben[0]  *x[1]   + oben[n-1]*x[n-1]  =  rs[0]     *
* oben[i-1]*x[i-1] + diag[i]  *x[i]   + oben[i]  *x[i+1]  =  rs[i]     *
* oben[n-1]*x[0]   + oben[n-2]*x[n-2] + diag[n-1]*x[n-1]  =  rs[n-1]   *
*                                                                      *
* for i=1,..,n-2. The matrix thus has the form :                       *
*                                                                      *
*                                                                      *
* ( diag[0]  oben[0]     0   .    .   .    .   .    0     oben[n-1] )  *
* (                                                                 )  *
* ( oben[0]  diag[1]  oben[1]     0                          0      )  *
* (                                                          .      )  *
* (    0     oben[1]  diag[2]  oben[2]     0                 .      )  *
* (    .   .        .        .        .        .             .      )  *
* (    .        .        .        .        .        .        .      )  *
* (    .            .         .        .        .        .   .      )  *
* (    .                 .        .        .        .        0      )  *
* (    .                     .         .        .        .          )  *
* (    0                          .        .        .     oben[n-2] )  *
* (                                   .           .        .        )  *
* ( oben[n-1]   0   .    .   .    .   .    0    oben[n-2] diag[n-1] )  *
.BE*)
*                                                                      *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* modus  kind of call for the function:                                *
*        = 0: factor matrix and solve system                           *
*        = 1: factor matrix only                                       *
*        = 2: solve system only. In this case, the factorization must  *
*             be available in oben, diag and rechts from zzytrdsy().   *
*             Use when solving the same system for several right hand  *
*             sides.                                                   *
* n      number of equations    (n > 2)                                *
* diag   modus = 0,1: [0..n-1] vector, the main diagonal of the matrix *
*        modus = 2:   [0..n-1] vector, the diagonal of D from          *
                      zzytrdsy())                                      *
* oben   modus = 0,1: [0..n-1] vector, with the codiagonal for indices *
*                     0,...,n-2  and A[0][n-1] in oben[n-1]            *
*        modus = 2:   [0..n-2] vector, containing information of the   *
*                     Cholesky factor R (see zzytrdsy() )              *
* rechts modus = 0,1: unused                                           *
*        modus = 2:   [0..n-3] vector, containing information about R  *
*        The  Cholesky factorization of A can be obtained by calling   *
*        this function with modus = 0 or modus = 1.                    *
* rs     modus = 0,2: [0..n-1] vector of the right hand side           *
*        modus = 1:   unused                                           *
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* diag  \  the Cholesky factors in condensed form                      *
* oben   > (see function zzytrdsy() )                                  *
* rechts/                                                              *
* rs       modus = 0,2: [0..n-1] vector with the solution              *
*          modus = 1:   unused                                         *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: no error                                                        *
* = 1: improper input  (n <= 2  or  modus != 0,1,2)                    *
* = 2: A not numerically strongly nonsingular                          *
* = 4: unusable factorization and modus = 2:  D is singular            *
* The determinant of A is given by :                                   *
*         det(A) = diag[0] * diag[1] * ... * diag[n-1].                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, zzytrdsy, lzytrdsy                                             *
.BA*)
***********************************************************************/
/*.BE*/

{
  int fehler;                           /* error code from zzytrdsy() */


  if (n <= 2)                                         /* n too small? */
    return 1;

  switch (modus)
  {
    case 0:                                     /* factor and solve   */
      fehler = zzytrdsy(n, diag, oben, rechts);
      if (fehler)
        return fehler;
      return lzytrdsy(n, diag, oben, rechts, rs);

    case 1:                                     /* factor only?       */
      return zzytrdsy(n, diag, oben, rechts);

    case 2:                                     /* solve only ?       */
      return lzytrdsy(n, diag, oben, rechts, rs);

    default:                                    /* invalid modus?     */
      return 1;
  }
}

/* ------------------------- END fzytrdsy.c ------------------------- */
