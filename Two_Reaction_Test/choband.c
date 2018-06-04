#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE choband.c ----------------------- */

/***********************************************************************
*                                                                      *
* Solving a linear system of equations with a symmetric, strongly      *
* ---------------------------------------------------------------      *
* nonsingular banded system matrix via the modified Cholesky method    *
* -----------------------------------------------------------------    *
*                                                                      *
* exported functions:                                                  *
*   - chobnd():  Solution of a linear system via Cholesky              *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Elmar Pohl (FORTRAN)                           *
* Realization:          Juergen Dietel, Computer center, RWTH Aachen   *
* Source:               FORTRAN program                                *
* Date:                 3.20.1992                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  REAL, max, min, MACH_EPS, FABS, ZERO, */
                       /*       FOUR                                  */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR                                */
#include <choband.h>   /*  for  chobnd                                */



/* ------------------------------------------------------------------ */

#ifdef DEBUG
/* print a REAL vector of length `n' (together with its name)         */
#define zeigrvek(vek, n)              \
  {                                   \
    int jjj;                          \
    printf("%-8s", #vek": ");         \
    for (jjj = 0; jjj < n; jjj++)     \
      printf("%8"LZP"g", vek[jjj]);   \
    printf("\n");                     \
  }

#else
#define zeigrvek(vek, n)
#endif



/* ------------------------------------------------------------------ */

static int chobdz               /* Factorize a condensed band matrix  */
/*.IX{chobdz}*/
                 (
                  int  n,       /* Size of the matrix ................*/
                  int  m,       /* half of its band width ............*/
                  REAL *ap[]    /* Input or factored matrices ........*/
                 )              /* Error code ........................*/

/***********************************************************************
* We factor a symmetric, strongly nonsingular banded matrix, given in  *
* condensed form, into the product  R' * D * R  using the modified     *
* Cholesky method. Here D is diagonal, and R is a unit upper banded    *
* matrix with the same band width as the given matrix.                 *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n   Size of the system matrix A                                      *
* m   half of the matix band width                                     *
* ap  [0..n-1,0..m] matrix, that contains the upper triangle           *
*     of A in condensed form (for details, see  chobnd())              *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ap  contains the Cholesky factorization  R' * D * R of A in condensed*
*     form, where R is unit upper triangular and banded, with band     *
*     width m, and D is diagonal.                                      *
*     The diagonal odf D is in the first column of ap. The subsequent  *
*     columns contain the m superdiagonals of R. R' is not stored.     *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: no error                                                        *
* = 2: memory insufficient                                             *
* = 3: The system matrix is numerically not strongly nonsingular       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, MACH_EPS, max, min, FABS, vminit, vmalloc, vmcomplete, vmfree, *
* VEKTOR, ZERO, FOUR                                                   *
***********************************************************************/

{
  REAL *z;             /* [0,...,n-1] aux vector for the row sums     */
  REAL eps;            /* relative error bound                        */
  REAL tmp;            /* aux variable for the matrix row sums, etc.  */
  int  i;              /* Loop variables i, j, k                      */
  int  j;
  int  k;
  void *vmblock;       /* List of dynamically allocated vectors       */


  vmblock = vminit();
  z = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 2;
  }

  eps = FOUR * MACH_EPS;


  for (i = 0; i < n; i++)  /* store row sums in z to check for        */
  {                        /* nonsingularity                          */
    tmp = ZERO;
    k   = min(n, i + m + 1);
    for (j = i; j < k; j++)
      tmp += FABS(ap[i][j - i]);
    for (j = max(0, i - m); j < i; j++)
      tmp += FABS(ap[j][i - j]);

    if (tmp == ZERO)                   /* Is the matrix near singular?*/
      return 3;
    z[i] = tmp;
  }


  for (j = 0; j < n; j++)            /* factor matrix a la Cholesky  */
  {
    for (i = max(0, j - m); i < j; i++)
    {
      tmp          = ap[i][j - i];
      ap[i][j - i] = tmp / ap[i][0];
      for (k = i + 1; k <= j; k++)
        ap[k][j - k] -= tmp * ap[i][k - i];
    }

    if (FABS(ap[j][0]) / z[j] <= eps)       /* Matrix numerically not */
      return 3;                             /* strongly nonsingular?  */
  }


  vmfree(vmblock);
  return 0;                       /* all ok, matrix positive definite */
}



/* ------------------------------------------------------------------ */

static int chobdl           /* Solve by updating and backsubstitution */
/*.IX{chobdl}*/
                 (
                  int  n,        /* Size of the matrix ...............*/
                  int  m,        /* Number of super diagonals ........*/
                  REAL *ap[],    /* condensed factorization ..........*/
                  REAL x[]       /* right hand side or solution ......*/
                 )               /* Error code .......................*/

/***********************************************************************
* Given the Cholesky factorization  A = R' * D * R from chobdz() of a  *
* symmetric, strongly nonsingular banded matrix A in condensed form,   *
* we find the solution vector X with                                   *
*                         A * X = RS                                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n   Size of the linear system                                        *
* m   half bandwidth of the system matrix                              *
* ap  [0..n-1,0..m] matrix, which contains the facors of A in          *
*     condensed form from chobdz()                                     *
* x   [0..n-1] vector, the right hand side of the linear system        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x  [0..n-1] vector, the solution                                     *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* Error code:                                                          *
* = 0: no error                                                        *
* = 5: unusable factorization: zero appears on the diagonal of D.      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, max, min                                                       *
***********************************************************************/

{
  int i;                                         /* Loop variables    */
  int j;
  int l;                                         /* end of loop value */


  for (j = 0; j < n; j++)                        /* updating :        */
    for (i = max(0, j - m); i < j; i++)          /* Solve R' * Z = RS */
      x[j] -= ap[i][j - i] * x[i];               /* Store Z in x      */
  zeigrvek(x, n);

  for (i = n; i-- != 0; )                        /* Backsubstitution: */
  {
    if (ap[i][0] == ZERO)                        /* invalid D ?       */
      return 5;
                                            /*                -1      */
    x[i] /= ap[i][0];                       /* Solve R * X = D   * Z  */
    l    =  min(n, i + m + 1);              /* Z is stored in x,      */
    for (j = i + 1; j < l; j++)             /* which is successively  */
      x[i] -= ap[i][j - i] * x[j];          /* filled with the        */
  }                                         /* solution vector        */


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int chobnd         /* Cholesky method for condensed band matrices ....*/
/*.IX{chobnd}*/
          (
           int  modus,    /* type of call: 0, 1, 2 ...................*/
           int  n,        /* size of the matrix ......................*/
           int  m,        /* half band width .........................*/
           REAL *ap[],    /* condensed matrix: Input or factorization */
           REAL rs[]      /* right hand side or solution .............*/
          )               /* Error code ..............................*/

/***********************************************************************
* Solve the linear system                                              *
*                         A * X = RS                                   *
* for a symmetric, strongly nonsingular definite banded matrix A in    *
* condensed form according to the modified Cholesky method             *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* modus  call type of function:                                        *
*        = 0: factor the matrix and solve the linear system            *
*        = 1: factor only                                              *
*        = 2: solve the system only. The factorization must be stored  *
*             in ap. This type of call is used when solving multiple   *
*             systems of equations for the same system matrix and      *
*             various right hand sides.                                *
* n      Size of the linear system                                     *
* m      half band width of A                                          *
* ap     [0..n-1,0..m] matrix, which contains for                      *
*        modus = 0: the upper codiagonals of A in condensed form       *
*          or  = 1: Condensing A is done by writing the codiagonals of *
*                   A as columns of ap. Specifically the main diagonal *
*                   A appears in the first column of ap.               *
*                   The following code condensed A into ap:            *
*                       for (i = 0; i < n; i++)                        *
*                           for (j = i; j < min(n, i + m + 1); j++)    *
*                               ap[i][j - i] = A[i][j];                *
*        modus = 2: the Cholesky factorization in condensed form       *
* rs     [0..n-1] vector with the right hand side for modus = 0, 2     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ap  the Cholesky factorization R' * D * R in condensed form, where   *
*     R is a unit upper band matrix with m codiagonals, and D is       *
*     diagonal. The diagonal of D occupies the first column of ap.     *
*     The other columns of ap contain the codiagonals of R. R' is not  *
*     stored explicitly.                                               *
* rs  [0..n-1] vector , the solution vector for  modus = 0, 2          *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: no error, matrix positive definite                              *
* = 1: wrong input parameter                                           *
* = 2: insufficient memory                                             *
* = 3: matrix numerically not strongly nonsingular                     *
* = 5: improper factorization in case  modus = 2:  D is singular.      *
* Upon a successful run, the determinant of A is available as :        *
*         det(A) = ap[0][0] * ap[1][0] * ... * ap[n-1][0].             *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* chobdz, chobdl, REAL                                                 *
.BA*)
***********************************************************************/
/*.BE*/

{
  int fehler;                             /* error code from chobdz() */


  if (n <= 0 || m < 0)                 /* improper values for n or m? */
    return 1;


  switch (modus)
  {
    case 0:
      fehler = chobdz(n, m, ap);  /* Compute factorization A = R'*D*R */

      if (fehler)                   /* Factorization impossible?      */
        return fehler;              /* send error message             */
      return chobdl(n, m, ap, rs);  /* solve linear system            */

    case 1:
      return chobdz(n, m, ap);      /* Factor  A = R' * D * R         */
                                    /* and send result                */

    case 2:
      return chobdl(n, m, ap, rs);  /* solve the linear system        */
                                    /* ans send solution              */

    default:                        /* improper value for modus?      */
      return 1;
  }
}

/* ------------------------- END choband.c -------------------------- */
