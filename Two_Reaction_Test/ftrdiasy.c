#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE ftrdiasy.c ------------------------ */

/***********************************************************************
*                                                                      *
* Solve a linear system with a tridiagonal, symmetric strongly         *
* ------------------------------------------------------------         *
* nonsingular system matrix using a root free Cholesky factorization   *
* ------------------------------------------------------------------   *
*                                                                      *
* exported function:                                                   *
*   - trdiasy():  solution of a linear system with a strongly          *
*                 nonsingular tridiagonal symmetric matrix             *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 mit 80387                          *
* Author:               Gisela Engeln-Muellges (FORTRAN)               *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               FORTRAN code                                   *
* Date:                 11.27.1992                                     *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  REAL, MACH_EPS, FABS, ZERO, ONE, FOUR */
#include <ftrdiasy.h>  /*  for  trdiasy                               */



/* ------------------------------------------------------------------ */

static int ztrdiasy
/*.IX{ztrdiasy}*/
                   (
                    int  n,
                    REAL diag[],
                    REAL oben[]
                   )

/***********************************************************************
* The tridiagonal, symmetric, strongly nonsingular matrix A is given   *
* by the two vectors diag and oben, denoting the diagonal and          *
* co-diagonal of A. (see function  ftrdiasy() ) A is decomposed into   *
* the product  R' * D * R  using the Cholesky method.                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n       number of equations (n > 2)                                  *
* diag    [0..n-1] vector, the main diagonal of A                      *
* oben    [0..n-2] vector, the co-diagonal of A                        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* diag \  the factors of the Cholesky decomposition of A in condensed  *
* oben /  form. The co-diagonal of the unit upper triangular R appears *
*         in oben[0..n-2], the diagonal of D is stored in diag[0..n-1].*
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* Error code.                                                          *
* = 0: all ok                                                          *
* = 2: matrix A not numerically strongly nonsingular                   *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, FOUR, MACH_EPS, ZERO, FABS, ONE                                *
***********************************************************************/

{
  REAL eps,      /* error bound for the relative error                */
       zeile,    /* norm of one matrix row                            */
       d,        /* reciprocal of zeile                               */
       oben_alt, /* in i loop: old off-diagonal element               */
       oben_neu; /* in i loop: new off-diagonal element               */
  int  i;        /* loop variable, number of the row                  */


  eps = FOUR * MACH_EPS;

  zeile = FABS(diag[0]) + FABS(oben[0]);   /* check first row for     */
  if (zeile == ZERO)                       /* A not to be strongly    */
    return 2;                              /* nonsingular             */
  d = ONE / zeile;
  if (FABS(diag[0]) * d <= eps)            /* A numerically not       */
    return 2;                              /* strongly nonsingular    */


  oben_alt =  oben[0];                    /* decompose A              */
  oben_neu =  oben[1];
  oben[0]  /= diag[0];

  for (i = 1; i < n; i++)
  {
    zeile   =  FABS(diag[i]) + FABS(oben_neu) + FABS(oben_alt);
    diag[i] -= oben_alt * oben[i - 1];

    if (zeile == ZERO)          /* check stron nonsingularity         */
      return 2;
    d = ONE / zeile;
    if (FABS(diag[i]) * d <= eps)
      return 2;

    if (i < n - 1)              /* last matrix row?                   */
      oben_alt =  oben[i],      /* => there is at least one more      */
      oben[i]  /= diag[i];      /*    sub-diagonal element in the next*/
                                /*    row                             */

    if (i < n - 2)              /* not the last but one row?          */
      oben_neu = oben[i + 1];   /* => there is at least one more      */
                                /*    upper co-diagonal element in    */
                                /*    next row                        */
    else                        /* last or last but one row?          */
      oben_neu = ZERO;          /* => no upper co-diagonal element    */
  }                             /*    in next row                     */


  return 0;
}



/* ------------------------------------------------------------------ */

static int ltrdiasy
/*.IX{ltrdiasy}*/
                   (
                    int  n,
                    REAL diag[],
                    REAL oben[],
                    REAL rs[]
                   )

/***********************************************************************
* Solve the linear system                                              *
*                        A * X = RS                                    *
* for a tridiagonal, symmetric, strongly nonsingular matrix A whose    *
* root free Cholesky decomposition  A  =  R' * D * R  is known from    *
* ztrdiasy().                                                          *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n     number of equations (n > 2)                                    *
* diag  [0..n-1] vector, the diagonal of D                             *
* oben  [0..n-2] vector, the co-diagonal of the unit upper triangular  *
*       matrix R                                                       *
* rs    [0..n-1] vector of the right hand side of the system           *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* rs    [0..n-1] solution vector                                       *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* error code.                                                          *
* = 0: no error                                                        *
* = 4: D is singular (zero on its diagonal)                            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO                                                           *
***********************************************************************/

{
  REAL rsi1;      /* aux storage for the previously computed component*/
                  /* of the solution before dividing by the diagonal  */
                  /* element (this helps perform the two updating     */
                  /* steps in one loop)                               */
  int  i;         /* Loop variable                                    */


  if (diag[0] == ZERO)
    return 4;
  rsi1  = rs[0];                                         /* update    */
  rs[0] = rsi1 / diag[0];

  for (i = 1; i < n; i++)
  {
    if (diag[i] == ZERO)
      return 4;
    rsi1  = rs[i] - oben[i - 1] * rsi1;
    rs[i] = rsi1 / diag[i];
  }


  for (i = n - 2; i >= 0; i--)              /* back substitute       */
    rs[i] -= oben[i] * rs[i + 1];


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int trdiasy      /* tridiagonal, symmetric linear system .............*/
/*.IX{trdiasy}*/
           (
            int  modus,   /* kind of call: 0, 1, 2 ...................*/
            int  n,       /* size of the system ......................*/
            REAL diag[],  /* main diagonal of A, or D in A = R'DR ....*/
            REAL oben[],  /* co-diagonal of A, or co-diagonal of R ...*/
            REAL rs[]     /* right hand side, or solution ............*/
           )              /* error code ..............................*/

/***********************************************************************
* Solve a linear system                                                *
*                        A * X = RS                                    *
* with a tridiagonal, symmetric, strongly non singular matrix A.       *
.BE*)
* A is given by its diagonal diag and its co-diagonal oben:            *
*                                                                      *
*  diag[0] *x[0]   + oben[0] *x[1]                       =  rs[0]      *
* oben[i-1]*x[i-1] + diag[i] *x[i]   +  oben[i] *x[i+1]  =  rs[i]      *
*                   oben[n-2]*x[n-2] + diag[n-1]*x[n-1]  =  rs[n-1],   *
*                                                                      *
* where i=1, ..., n-2. Thus the symmetric system matrix A is given as: *
*                                                                      *
*                                                                      *
* ( diag[0]  oben[0]     0   .    .   .    .   .    .   .    0      )  *
* (                                                          .      )  *
* ( oben[0]  diag[1]  oben[1]     0                          .      )  *
* (                                                          .      )  *
* (    0     oben[1]  diag[2]  oben[2]     0                 .      )  *
* (    .   .        .        .        .        .             .      )  *
* (    .        .        .        .        .        .        .      )  *
* (    .            .         .        .        .        .   .      )  *
* (    .                 .        .        .        .        0      )  *
* (    .                     .         .        .        .          )  *
* (    .                          .        .        .     oben[n-2] )  *
* (    .                              .           .        .        )  *
* (    0   .    .   .    .   .    .   .    0    oben[n-2] diag[n-1] )  *
*                                                                      *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* modus  kind of call:                                                 *
*        = 0: factor matrix and solve the linear system                *
*        = 1: factor matrix only                                       *
*        = 2: solve system only. The factorization of A is available   *
*             in oben and diag bereitstehen. This kind of call is used *
*             when wanting to solve for several right hand sides for   *
*             one A, such as when inverting A.                         *
* n      number of equations (n > 2)                                   *
* diag   modus = 0,1: [0..n-1] diagonal vector of A                    *
*        modus = 2:   [0..n-2] diagonal vector of D in the Cholesky    *
*                     factorization of A from ztrdiasy().              *
*        To obtain the root free Cholesky factorization of A, call     *
*        this function with  modus=0 or modus=1.                       *
* oben   modus = 0,1: [0..n-2] upper co-diagonal of A                  *
*        modus = 2:   [0..n-1] upper co-diagonal of R from ztrdiasy()  *
* rs     modus = 0,2: [0..n-1] right hand side                         *
*        modus = 1:   unused                                           *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* diag \  the Cholesky factors of A                                    *
* oben /                                                               *
* rs      modus = 0,2: [0..n-1] solution vector                        *
*         modus = 1:   unused                                          *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* Error code.                                                          *
* = 0: no error                                                        *
* = 1: wrong input  (n <= 2  or  modus != 0,1,2)                       *
* = 2: Matrix A not strongly nonsingular                               *
* = 4: invalid decomposition and  modus = 2:  D is singular            *
* The determinant of A can be computed as:                             *
*         det(A) = diag[0] * diag[1] * ... * diag[n-1].                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ztrdiasy, ltrdiasy                                             *
.BA*)
***********************************************************************/
/*.BE*/

{
  int fehler;                         /* Return value from zzytrdsy() */


  if (n <= 2)                                        /* n too small? */
    return 1;

  switch (modus)
  {
    case 0:                                     /* factor and solve?  */
      fehler = ztrdiasy(n, diag, oben);
      if (fehler)
        return fehler;
      return ltrdiasy(n, diag, oben, rs);

    case 1:                                          /* factor only ? */
      return ztrdiasy(n, diag, oben);

    case 2:                                     /* solve only?        */
      return ltrdiasy(n, diag, oben, rs);

    default:                                    /* modus incorrect?   */
      return 1;
  }
}

/* -------------------------- END ftrdiasy.c ------------------------ */
