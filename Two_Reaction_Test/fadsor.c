#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE fadsor.c ------------------------ */

/***********************************************************************
*                                                                      *
* Iterative solution for a linear system of equations with the         *
* -------------------------------------------------------------        *
* adaptive SOR Method                                                  *
* --------------------                                                 *
*                                                                      *
* exported function :                                                  *
*   - adsor():  adaptive SOR method to solve linear equations          *
*               iteratively                                            *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 mit 80387                          *
* Author:               Gisela Engeln-Muellges (FORTRAN)               *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing FORTRAN code                          *
* Date:                 10.7.1992                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  REAL, ZERO, TWO, NULL, EIGHT, max,   */
                        /*       MACH_EPS, FABS, sqr, SQRT, ONE,      */
                        /*       seidel, TRUE, norm_max, copy_vector  */
#include <vmblock.h>    /*  for  vminit, vmalloc, VEKTOR, vmcomplete, */
                        /*       vmfree                               */
#include <fadsor.h>     /* adsor                                      */



/* ------------------------------------------------------------------ */

static void sub_vector
/*.IX{sub\unt vector}*/
                      (
                       REAL differenz[],
                       REAL minuend[],
                       REAL subtrahend[],
                       int  n
                      )

/***********************************************************************
* subtract vector subtrahend from vector minuend; store result in      *
* the vector differenz                                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
*                                                                      *
***********************************************************************/

{
  for (n--; n >= 0; n--)
    *differenz++ = *minuend++ - *subtrahend++;
}



/* ------------------------------------------------------------------ */

static void seidel               /* Gauss-Seidel method ..............*/
/*.IX{seidel}*/
                  (
                   int  n,       /* size of matrix ...................*/
                   REAL *mat[],  /* modified input matrix ............*/
                   REAL b[],     /* modified right hand side .........*/
                   REAL omega,   /* Relaxation coefficient ...........*/
                   REAL x[]      /* solution vector ..................*/
                  )

/***********************************************************************
* Perform one Gauss-Seidel Iteration step for a given relaxation       *
* coefficient                                                          *
*                                                                      *
* Input parameters:                                                    *
* ================                                                     *
* n      size of the linear system                                     *
* mat    [0..n-1,0..n-1] system matrix, modified to have diagonal      *
*        entries equal to 1                                            *
* b      [0..n-1] vector, the right hand side (modified)               *
* omega  Relaxation coefficient                                        *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* x      [0..n-1] vector, the solution                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
*                                                                      *
***********************************************************************/

{
  int  i,                              /* Loop indices                */
       j;                              /*                             */
  REAL tmp;                            /* aux variable for summing    */

  for (i = 0; i < n; i++)
  {
    for (tmp = b[i], j = 0; j < n; j++)
      tmp -= mat[i][j] * x[j];
    x[i] += omega * tmp;
  }
}


/* ------------------------------------------------------------------ */
/*.BA*/

int adsor                         /* adaptive SOR method .............*/
/*.IX{adsor}*/
         (
          int  crit,              /* Convergence criterion (0,1,2,3) .*/
          int  n,                 /* size of matrix ..................*/
          REAL *mat[],            /* matrix ..........................*/
          REAL b[],               /* right hand side .................*/
          REAL *omega,            /* Relaxation coefficient ..........*/
          REAL x[],               /* solution ........................*/
          REAL residu[],          /* Residuum ........................*/
          int  *iter,             /* number of iterations ............*/
          int  l,                 /* number of steps before adapting  */
                                  /* new coefficient .................*/
          REAL eps,               /* error bound  ....................*/
          int  maxit,             /* Maximal number of iterations ....*/
          int  methode            /* method :  (0,1,2) ...............*/
         )                        /* error code ......................*/

/***********************************************************************
* adsor solves the linear system  mat * x = b  iteratively.            *
* Here  mat is nonsingular n x n, b  is the right hand side and x is   *
* the solution.                                                        *
*                                                                      *
* adsor uses Gauss-Seidel iteration with relaxation. Here the          *
* relaxation coefficient is periodically adjusted. (adaptive SOR       *
* method). By choosing parameters appropriately the ordinary Gauss-    *
* Seidel method or a nonadaptive SOR method can be selected as well.   *
.BE*)
*                                                                      *
* Application:                                                         *
* ============                                                         *
* Solve linear systems that satisfy one of the criteria: row or column *
* sum criterion or criterion of Schmidt-von-Mises. In this case        *
* convergence is certain.                                              *
*                                                                      *
* Note :                                                               *
* ======                                                               *
* For the adaptive SOR method (methode = 0) we advise to take I = 4    *
* or l=5.                                                              *
* If the optimal relaxation coefficient  omega_opt  is known, we       *
* advise to select methode = 1 with  omega = omega_opt.                *
*                                                                      *
* Input parameters:                                                    *
* ================                                                     *
* crit     swich which selects the convergence test criterion :        *
*          = 1:   row sum criterion                                    *
*          = 2:   column sum criterion                                 *
*          = 3:   criterion of Schmidt-von-Mises                       *
*          other: no check                                             *
* n        size of system                                              *
* mat      [0..n-1,0..n-1] system matrix                               *
* b        [0..n-1] right hand side                                    *
* omega    if methode = 1: the optimal relaxation coefficient          *
*          (0 < omega < 2), otherwise not used at input                *
* x        [0..n-1] starting vector for iteration                      *
* l        number of iterations after which the relaxation coefficient *
*          is to be recomputed/adapted                                 *
* eps      error bound: Iterations are stopped if the relative error   *
*          in max norm does not exceed eps                             *
* maxit    Maximal number of iterations allowed                        *
* methode  Switch for selection of method :                            *
*           = 0: adaptive SOR method                                   *
*           = 1: SOR method with a given constant relaxation           *
*                coefficient                                           *
*           = 2: Gauss-Seidel method                                   *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        [0..n-1] vector, the solution                               *
* residu   [0..n-1] residual vector  b - mat * x for last iteration    *
* iter     number of iterations executed                               *
* mat      modified system matrix:                                     *
*          mat[i][j] /= mat[i][i]  (i,j=0,...,n-1)                     *
* b        modified right hand side:                                   *
*          b[i] /= a[i][i]         (i=0,...,n-1)                       *
* omega    methode = 0: adaptively computed relaxation coefficient     *
*          methode = 1: optimal relaxation coefficient                 *
*          methode = 2: the number 1                                   *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* =  0: solution has been found                                        *
* =  1: n < 1 or l < 1 or eps <= 0 or maxit < 1 or methode < 0 or      *
*       methode > 2 or (if methode = 1: omega <= 0 or omega >= 2)      *
* =  2: mat or b or x invalid                                          *
* =  3: one diagonal entry or one row of mat is zero                   *
* =  4: desired accuracy eps not reached after maxit iterations; check *
*       answer anyway.                                                 *
* =  5: lack of memory                                                 *
* = 11: column sum criterion not met                                   *
* = 12: row sum criterion not met                                      *
* = 13: criterion of Schmidt-von-Mises not met                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, TWO, NULL, vminit, vmalloc, VEKTOR, vmcomplete, vmfree,  *
* EIGHT, MACH_EPS, FABS, sqr, SQRT, ONE, seidel, TRUE, norm_max,       *
* copy_vector, max, sub_vector                                         *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i,             /* Loop indices                                 */
       j,             /*                                              */
       rc;            /* error code of this function                  */
  REAL q,             /* estimate for the spectral radius of the      */
                      /* iteration matrix                             */
       relerr,        /* relative error when checking matrix          */
       *diff0,        /* [0..n-1] vector of the difference between two*/
                      /* previous iterates                            */
       *diff1,        /* [0..n-1] vector of the difference between two*/
                      /* most recent iterates                         */
       *hilf,         /* [0..n-1] aux vector for q                    */
       tmp;           /* aux variable                                 */
  void *vmblock;      /* List of dynamically allocated vectors        */


  /* -------------------- check input ------------------------------- */

  if (n < 1 || l < 1 || eps <= ZERO || maxit < 1 ||
      methode < 0 || methode > 2)
    return 1;

  if (methode == 1 && (*omega <= ZERO || *omega >= TWO))
    return 1;

  if (mat == NULL)
    return 2;

  for (i = 0; i < n; i++)
    if (mat[i] == NULL)
      return 2;

  if (x == NULL || b == NULL || residu == NULL)
    return 2;


  vmblock = vminit();           /* allocate three dynamic aux vectors */
  diff0 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  diff1 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  hilf  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 5;
  }

  relerr = EIGHT * MACH_EPS;        /* check matrix for singularity   */
  for (i = 0; i < n; i++)           /* or zero diagonal entries       */
  {
    for (tmp = FABS(mat[0][0]), j = 1; j < n; j++)
      tmp += FABS(mat[i][j]);
    if (tmp == ZERO || FABS(mat[i][i]) / tmp < relerr)
      return 3;
  }

  switch (methode)         /* initialize  l and omega depending on    */
  {                        /* chosen method                           */
    case 0: *omega = ONE;
            break;
    case 1: l     = maxit;
            break;
    case 2: l     = maxit;
            *omega = ONE;
  }

  for (i = 0; i < n; i++)       /* modify matrix                      */
  {
    for (tmp = ONE / mat[i][i],
         j = 0; j < n; j++)
      mat[i][j] *= tmp;         /* modify b accordingly               */
    b[i] *= tmp;
  }

  switch (crit)                  /* check chosen convergence criterion*/
  {
    case 1:  for (i = 0; i < n; i++)       /* row sum criterion ?     */
             {
               for (tmp = ZERO, j = 0; j < n; j++)
                 tmp += FABS(mat[i][j]);
               if (tmp >= TWO)
                 return 11;
             }
             break;

    case 2:  for (j = 0; j < n; j++)       /* column sum criterion ?  */
             {
               for (tmp = ZERO, i = 0; i < n; i++)
                 tmp += FABS(mat[i][j]);
               if (tmp >= TWO)
                 return 12;
             }
             break;

    case 3:  for (tmp = ZERO, i = 0; i < n; i++)
               for (j = 0; j < n; j++)     /* criterion of            */
                 tmp += sqr(mat[i][j]);    /* Schmidt-von-Mises?      */
             tmp = SQRT(tmp - ONE);
             if (tmp >= ONE)
               return 13;
             break;

    default: break;                        /* no checking ?           */
  }


  /* ------------------- prepare iteration -------------------------- */

  copy_vector(residu,     /* residu containes the previous iterate    */
              x, n);      /* during the iteration. thus residu now    */
                          /* stores the given starting vector.        */
  seidel(n, mat, b, *omega, x);   /* compute next iterate x           */
  sub_vector(diff0,               /* find difference of two most      */
             x, residu, n);       /* iterates                         */


  /* -------------------------- iterate ----------------------------- */

  for (rc = 0, *iter = 1; TRUE; )
  {
    if (norm_max(diff0, n) <= eps * norm_max(x, n))   /* solution good*/
      break;                                          /* enough ?     */

    if (*iter >= maxit)                /* excessive iteration ?       */
    {
      rc = 4;                          /* report error                */
      break;
    }

    copy_vector(residu, x, n);    /* store previous iterate in residu */
    seidel(n, mat, b, *omega, x); /* compute new iterate  x           */
    sub_vector(diff1,             /* compute their difference         */
               x, residu, n);

    if (++*iter % l == 0)          /* l steps with constant relaxation*/
    {                              /* coefficient performed ?         */
      for (i = 0; i < n; i++)               /* estimate spectral      */
        hilf[i] = ((FABS(diff0[i]) <        /* radius q of iteration  */
                    TWO * MACH_EPS))        /* matrix and use for new */
                  ? ONE                     /* relaxation coefficient */
                  : (diff1[i] / diff0[i]);
      if ((q = norm_max(hilf, n)) <= ONE)
        q      = max(q, *omega - ONE),
        *omega = TWO / (ONE +
                SQRT(ONE - sqr((q + *omega - ONE) / *omega) / q));
    }

    copy_vector(diff0,       /* restore new difference vector as old  */
                diff1, n);   /* one for next step                     */
  }


  for (i = 0; i < n; i++)                 /* compute residual vector  */
 {
    for (tmp = b[i], j = 0; j < n; j++)
      tmp -= mat[i][j] * x[j];
    residu[i] = tmp;
  }


  vmfree(vmblock);
  return rc;
}

/* --------------------------- END fadsor.c ------------------------- */
