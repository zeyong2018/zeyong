#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE hrmsplin.c ------------------------ */

/***********************************************************************
*                                                                      *
* Computes interpolating Hermite polynomial splines                    *
* -----------------------------------------------------------------    *
*                                                                      *
* Programming language:  ANSI C                                        *
* Compiler:              Turbo C 2.0                                   *
* Computer:              IBM PS/2 70 with 80387                        *
* Source:                Equivalent Quick-basic code                   *
* Author:                Elmar Pohl (QuickBASIC), Dorothee Seesing (C) *
* Adaptation:            Juergen Dietel, Computer Center, RWTH Aachen  *
* Date:                  2.2.1993                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  sqr, MAXROOT, NULL, SQRT, FABS, REAL, */
                       /*       ZERO, ONE, TWO, THREE, HALF, FIVE,    */
                       /*       FOUR, TEN, EIGHT, NINE, SIX, TRUE,    */
                       /*       FALSE                                 */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR                                */
#include <u_proto.h>   /*  for  trdiag, tzdiag                        */
#include <hrmsplin.h>  /*  for  hermit, parmit                        */
/*.BA*/



/*.FE{C 10.2.2}
     {Computation of Non-Parametric Hermite Splines}
     {Computation of Non-Parametric Hermite Splines}*/

/*.BE*/
/* ------------------------------------------------------------------ */

static REAL fdext         /* aux function for hermit() ...............*/
/*.IX{fdext}*/
                 (
                  REAL a1,
                  REAL a2,
                  REAL a3,
                  REAL b1,
                  REAL b2,
                  REAL b3,
                  REAL rec1,
                  REAL rec2
                 )

/***********************************************************************
* The return value is one entry of the right hand side of the linear   *
* system in  hermit().                                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, TEN, FOUR                                                      *
***********************************************************************/

{
  REAL hilf,
       rec1h,
       rec2h;

  rec1h =  rec1 * rec1;
  rec2h =  rec2 * rec2;
  hilf  =  TEN * ((a3 - a2) * rec2 * rec2h - (a2 - a1) * rec1 * rec1h);
  hilf  += FOUR * (b1 * rec1h - (REAL)1.5 * (rec2h - rec1h) * b2 -
           b3 * rec2h);

  return hilf;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int hermit         /* non-parametric Hermite spline ..................*/
/*.IX{hermit}*/
          (
           int  m,            /* number of  nodes ....................*/
           REAL x[],          /* nodes: x-values .....................*/
           REAL y[],          /*        y-values .....................*/
           REAL y1[],         /* first derivative vector at nodes ....*/
           int  marg_cond,    /* type of boundary condition ..........*/
           REAL marg_0,       /* left boundary condition .............*/
           REAL marg_n,       /* right boundary condition ............*/
           int  save,         /* save dynamic aux arrays ? ...........*/
           REAL c[],          /* Spline coefficients of (x-x[i])^2 ...*/
           REAL d[],          /* Spline coefficients of (x-x[i])^3 ...*/
           REAL e[],          /* Spline coefficients of (x-x[i])^4 ...*/
           REAL f[]           /* Spline coefficients of (x-x[i])^5 ...*/
          )                   /* error code ..........................*/

/***********************************************************************
* We compute the coefficients of the Hermite interpolating spline for  *
* the given triples of nodes and first derivatives:                    *
*                    (x[i], y[i], y1[i]), i = 0,...,m-1.               *
* The boundary type can be specified in marg_cond . The x[i] must be   *
* strictly increasing. When calling hermit for the same set of x-values*
* but for different y- and/or y1-values one can save the time for      *
* forming and updating the same linear system repeatedly by setting    *
* the control parameter save to a value different from 0 which saves   *
* the LU decomposition for the next call of hermit.                    *
.BE*)
* IMPORTANT: In order to free the storage for the aux vectors and thus *
* --------   to avoid computing with the wrong LU factorization in sub-*
*            sequent calls, the last call of a sequence of related runs*
*            must be executed with save = 0 instead of save = 1.       *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* m:          number of nodes ( m > 2)                                 *
* x:          [0..m-1] vector of nodes: x-values                       *
* y:          [0..m-1] vector, y-values                                *
* y1:         [0..m-1] vector, first derivatives                       *
* marg_cond:  type of boundary condition:                              *
*             = 1: periodic spline                                     *
*             = 2: natural spline                                      *
*             = 3: 2nd derivatives at end points prescribed            *
*             = 4: curvature radii at the end points prescribed        *
*             = 5: 3rd derivative at end points prescribed             *
* marg_0: \   boundary conditions for                                  *
* marg_n: /   marg_cond = 3, 4, 5                                      *
* save:       Flag, which indicates whether the LU decomposition can be*
*             used for a subsequent call, provided the x-values stay   *
*             the same. Usually save = 0. If save = 1 for the same set *
*             of x-values one can save about 4*m flops by not recom-   *
*             puting the factorization for the tridiagonal system      *
*             matrix.                                                  *
*             The last call of a run, however, must be executed with   *
*             save = 0, otherwise certain storage is erroneaously kept *
*             and will foul up future runs of hermit.                  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* c: \   [0..m-2] or [0..m-1] vectors of the spline coefficients:      *
* d:  \  s(x)  =   a[i] + b[i] * (x - x[i]) +                          *
* e:  /          + c[i] * (x - x[i]) ^ 2 + d[i] * (x - x[i]) ^ 3 +     *
* f: /           + e[i] * (x - x[i]) ^ 4 + f[i] * (x - x[i]) ^ 5.      *
*        Here a denotes the vector y, and b denotes y1.                *
*        c and d have m entries (as do a and b).                       *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: no error                                                        *
* = 1: type of boundary condition out of range                         *
* = 2: m < 3                                                           *
* = 3: marg_cond = 1: not periodic                                     *
* = 4: x-values not monotone                                           *
* = 5: marg_cond = 4: one of the boundary values is  0.                *
* = 6: not enough storage for aux vectors                              *
* = 7: error in trdiag() or tzdiag()                                   *
* In case of error, the output is indeterminant and the aux storage is *
* freed.                                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* fdext, sqr, REAL, SQRT, vminit, vmalloc, vmcomplete, vmfree, VEKTOR, *
* trdiag, tzdiag, NULL, TRUE, FALSE, ZERO, ONE, TWO, THREE, HALF,      *
* FIVE, TEN, EIGHT, NINE, SIX                                          *
.BA*)
***********************************************************************/
/*.BE*/

{

#define ciao(fehler)          /* prepare to clear storage in case of */\
  {                           /* error before closing hermit()       */\
    vmfree(vmblock); /* free storage for aux vectors and report      */\
    vmblock = NULL;                                                    \
    return fehler;   /* pass error code on to user                   */\
  }

  static
    void *vmblock = NULL;    /* List of dynamically allocated vectors */
                             /* NULL indicates a first call of the    */
                             /* routine                               */

  static
    REAL *h,         /* [0..m-2] vector of length of node intervals   */
         *lower,     /* [0..m-2] lower subdiagonal of system/LU factor*/
                     /* after trdiag() or tzdiag()                    */
         *diag,      /* [0..m-2] main diagonal of system/LU factoriz. */
                     /* in trdiag() or tzdiag()                       */
         *upper,     /* [0..m-2] upper co-diagonal/LU factorization   */
         *lowrow,    /* [0..m-4] lowest row of system matrix/LU factor*/
         *ricol;     /* [0..m-4] right hand side vector/LU factor     */
                     /* update in tzdiag()                            */
    int  n,                  /* m - 1, Index of last node             */
         nm1,                /* n - 1                                 */
         nm2,                /* n - 2                                 */
         fehler,             /* error code of trdiag() or tzdiag()    */
         i,                  /* Loop variable                         */
         erster_aufruf;      /* Flag, keeping track of the first call */
                             /* of a run for the same x-values        */
    REAL alpha,              /* Parameter of the linear system that   */
         beta1 = ZERO,       /* depend on the boundary conditions:    */
         beta2 = ZERO,       /* (alpha affects the system matrix,     */
                             /* beta1, beta2 the right hand side)     */
         y21 = ZERO,         /* c[0] for given end point curvature    */
         y2n = ZERO,         /* c[n] for given end point curvature    */
         a1, a2, a3,         /* aux storage for y                     */
         b1, b2, b3,         /* aux storage for y1                    */
         hi,                 /* h[i]                                  */
         hsq,                /* h[0], h[nm2], y1[0], y1[n] ^ 2        */
         rec1, rec2;         /* aux variables, reciprocals of h       */

  n   = m   - 1;
  nm1 = n   - 1;
  nm2 = nm1 - 1;


  if (n < 2)                              /* less than 3 nodes ?      */
    ciao(2);

  if (marg_cond < 1 || marg_cond > 5)     /* improper end point       */
    ciao(1);                              /* condition ?              */

  if (marg_cond == 1)                         /* periodic Spline?     */
    if (y[0] != y[n] || y1[0] != y1[n])       /* non perodic boundary */
      ciao(3);                                /* values ?             */

  if (marg_cond == 4)                      /* curvature radii given ? */
    if (marg_0 == ZERO || marg_n == ZERO)  /* one radius equal to 0 ? */
      ciao(5);


  /* 1st call: allocate storage for aux vectors: 4 times [0..n-1]     */
  /* vectors (if periodic,  2 times [0..n-3] vectors additionally)    */

  if (vmblock == NULL)                       /* first call of a run ? */
  {
    erster_aufruf = TRUE;
    #define MYALLOC(l)  (REAL *)vmalloc(vmblock, VEKTOR, (l), 0)
    vmblock = vminit();                     /* allocate storage       */
    h     = MYALLOC(n);
    lower = MYALLOC(n);
    diag  = MYALLOC(n);
    upper = MYALLOC(n);
    if (marg_cond == 1)                   /* periodic spline with     */
      if (n > 2)                          /* enough nodes ?           */
        lowrow = MYALLOC(n - 2),          /* allocate extra vectors   */
        ricol  = MYALLOC(n - 2);
    #undef MYALLOC
  }
  else
    erster_aufruf = FALSE;
  if (! vmcomplete(vmblock))                      /* lack of memory ? */
    ciao(6);

  if (erster_aufruf)
  {
    for (i = 0; i <= nm1; i++)       /* compute interval length and   */
    {                                /* check for increasing x-values */
      hi = x[i + 1] - x[i];
      if (hi <= ZERO)
        ciao(4);
      h[i] = hi;
    }


    /* --------- form system matrix --------------------------------- */

    if (marg_cond == 5)
      alpha = EIGHT / NINE;
    else
      alpha = ONE;

    if (n == 2)
      diag[0] = THREE * alpha * (ONE / h[0] + ONE / h[1]);
    else
    {

      for (rec1 = alpha / h[0], i = 0; i < nm2; i++, rec1 = rec2)
        rec2     = ONE / h[i + 1],
        diag[i]  = THREE * (rec1 + rec2),
        upper[i] = lower[i + 1] = -rec2;
      diag[nm2] = THREE * (ONE / h[nm2] + alpha / h[nm1]);
    }

    if (marg_cond == 1)                       /* periodic spline?     */
    {                                         /* In this case         */
      rec1 = ONE / h[nm1];                    /* correct elements in  */
      rec2 = ONE / h[0];                      /* diag, upper, lower   */
      diag[nm1] = THREE * (rec1 + rec2);
      if (n == 2)
        upper[0] = lower[1] = -rec1 - rec2;   /* and initialize       */
      else                                    /* lowrow[0], ricol[0]  */
        lowrow[0] = ricol[0] = -rec2,
        lower[nm1] = upper[nm2] = -rec1;
    }
  }


  switch (marg_cond)             /* compute beta1 and beta2 from the  */
  {                              /* end point conditions              */
    case 1:                            /* periodic spline?            */
    case 2:                            /* natural spline ?            */
      beta1 = beta2 = ZERO;
      break;
    case 3:                            /* 2nd boundary derivatives ?  */
      beta1 = HALF * marg_0 / h[0];
      beta2 = HALF * marg_n / h[nm1];
      break;
    case 4:                            /* curvature radii givan ?     */
      hsq = sqr(y1[0]);
      y21 = HALF * SQRT(sqr((ONE + hsq)) * (ONE + hsq)) / marg_0;
      hsq = sqr(y1[n]);
      y2n = HALF * SQRT(sqr((ONE + hsq)) * (ONE + hsq)) / marg_n;
      beta1 = y21 / h[0];
      beta2 = y2n / h[nm1];
      break;
    case 5:                            /* 3rd boundary derivatives ?  */
      hsq = sqr(h[0]);
      beta1 = TEN * (y[1] - y[0]) / THREE / h[0] / hsq -
              TWO * (TWO * y1[1] + THREE * y1[0]) / THREE / hsq -
              marg_0 / (REAL)18.0;
      hsq = sqr(h[nm1]);
      beta2 = -TEN * (y[n] - y[nm1]) / THREE / h[nm1] / hsq +
              TWO * (THREE * y1[n] + TWO * y1[nm1]) / THREE / hsq +
              marg_n / (REAL)18.0;
      break;
  }


  rec1 = ONE / h[0];                  /* form right hand side         */
  a1   = y[0];
  a2   = y[1];
  b1   = y1[0];
  b2   = y1[1];
  for (i = 0; i <= nm2; i++, a1 = a2, b1 = b2, a2 = a3, b2 = b3,
                        rec1 = rec2)
    rec2 = ONE / h[i + 1],
    a3   =  y[i + 2],
    b3   = y1[i + 2],
    c[i] = fdext(a1, a2, a3, b1, b2, b3, rec1, rec2);
  c[0]   += beta1;
  c[nm2] += beta2;
  if (marg_cond == 1)                         /* periodic spline ?    */
    c[nm1] = fdext(a1, a2, y[1], b1, b2, y1[1], rec1, ONE / h[0]);


  switch (n)                                  /* solve linear system  */
  {
    case 2:
      if (marg_cond == 1)                     /* periodic spline?     */
        c[2] = (c[1] - lower[1] * c[0] / diag[0]) /
               (diag[1] - lower[1] * upper[0] / diag[0]),
        c[1] = (c[0] - upper[0] * c[2]) / diag[0];
      else
        c[1] = c[0] / diag[0];
      break;
    default:
      if (marg_cond == 1)                     /* periodic spline?     */
        fehler = tzdiag(n, lower, diag, upper, lowrow, ricol,
                        c, ! erster_aufruf);
      else
        fehler = trdiag(nm1, lower, diag, upper, c, ! erster_aufruf);
      if (fehler != 0)
        ciao(7);
      for (i = n; i != 0; i--)  /* shift solution vector one position */
        c[i] = c[i - 1];        /* down                               */
  }


  switch (marg_cond) /*depending on end point data, adjust c[0], c[n] */
  {
    case 1:
      c[0] = c[n];
      break;
    case 2:
      c[0] = c[n] = ZERO;
      break;
    case 3:
      c[0] = HALF * marg_0;
      c[n] = HALF * marg_n;
      break;
    case 4:
      c[0] = y21;
      c[n] = y2n;
      break;
    case 5:
      c[0] = beta1 * h[0]   + c[1]   / THREE;
      c[n] = beta2 * h[nm1] + c[nm1] / THREE;
  }


  /* -- compute remaining  d[i], e[i] and f[i] --------------------- */

  for (i = 0; i <= nm1; i++)
  {
    d[i] = TEN * (y[i + 1] - y[i]) / h[i] -
           TWO * (TWO * y1[i + 1] + THREE * y1[i]);
    d[i] = (d[i] / h[i] + c[i + 1] - THREE * c[i]) / h[i];
  }

  d[n] = d[nm1] - (TWO * (y1[n] - y1[nm1]) / h[nm1] -
                   TWO * (c[n] + c[nm1])) / h[nm1];

  for (i = 0; i <= nm1; i++)
    hi   = h[i],
    e[i] = (HALF * (y1[i + 1] - y1[i]) / hi - c[i]) / hi,
    e[i] = (e[i] - (REAL)0.25 * (d[i + 1] + FIVE * d[i])) / hi,
    f[i] = (((c[i + 1] - c[i]) / hi - THREE * d[i]) /
             hi - SIX * e[i]) / hi / TEN;


  if (!save)                               /* save aux vectors ?    */
    ciao(0);                               /* (last call of a run)? */

  return 0;
#undef ciao
}
/*.BA*/



/*.FE{C 10.2.3}
     {Computation of Parametric Hermite Splines}
     {Computation of Parametric Hermite Splines}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int parmit              /* parametric Hermite splines ................*/
/*.IX{parmit}*/
          (
           int  m,           /* number of nodes ......................*/
           REAL x[],         /* nodes : x-values .....................*/
           REAL y[],         /*         y-values .....................*/
           int  richt,       /* type of derivative ...................*/
           REAL xricht[],    /* Tangent or normal direction or only   */
           REAL yricht[],    /* dy/dx in yricht ......................*/
           int  marg,        /* type of end point condition ..........*/
           REAL corn_1[],    /* left hand end point condition ........*/
           REAL corn_n[],    /* right hand end point condition .......*/
           REAL cx[],        /* x spline coeffic. for (t-t[i])^2 .....*/
           REAL dx[],        /* x spline coeffic. for (t-t[i])^3 .....*/
           REAL ex[],        /* x spline coeffic. for (t-t[i])^4 .....*/
           REAL fx[],        /* x spline coeffic. for (t-t[i])^5 .....*/
           REAL cy[],        /* y spline coeffic. for (t-t[i])^2 .....*/
           REAL dy[],        /* y spline coeffic. for (t-t[i])^3 .....*/
           REAL ey[],        /* y spline coeffic. for (t-t[i])^4 .....*/
           REAL fy[],        /* y spline coeffic. for (t-t[i])^4 .....*/
           REAL t[],         /* Parameters nodes .....................*/
           REAL xt[],        /* normalized tangent vectors (x comp.) .*/
           REAL yt[]         /* normalized tangent vectors (y comp.) .*/
          )                  /* error code ...........................*/

/***********************************************************************
* Compute the coefficients for a  parametric Hermite interpolating     *
* spline.                                                              *
* The parameter marg designates the kind of end point condition present*
* The direction of the curve can be given by its tangent, normal or its*
* derivatives dy/dx.                                                   *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* m:       number of nodes (n > 2)                                     *
* x: \     [0..m-1] vectors of nodes                                   *
* y: /                                                                 *
* richt:   control parameter for derivative data:                      *
*          = 1: tangent vectors given in (xricht[i], yricht[i])        *
*          = 2: normal vectors given in (xricht[i], yricht[i])         *
*          = 3: derivatives  y'(x)  in yricht[i]                       *
*          The actual length of these vectors (xricht[i], yricht[i])   *
*          has no effect on the computations, as we normalize them.    *
*          On output, (xt, yt) contain the normalized data.            *
* xricht:  [0..m-1] x-components for tangents or normals (richt = 1,2) *
* yricht:  [0..m-1] y-components for tangents or normals (richt = 1,2) *
*          or  y'(x) (for richt = 3).                                  *
*          A derivative that is larger in magnitude than MAXROOT, i.e.,*
*          the square root of the largest representable number, is     *
*          interpreted as infinite, and the spline is given a vertical *
*          tangent there.                                              *
* marg:    kind of end point condition:                                *
*          = 1: periodic spline                                        *
*          = 2: natural spline                                         *
*          = 3: 2nd derivatives y''(x) given at end points, i.e.,      *
*               y''(x[0]) in corn_1[0] and y''(x[m-1]) in corn_n[0]    *
*          = 4: 2nd derivatives of the spline components wrt. curve    *
*               parameter are given at the end points:                 *
*                   ..                       ..                        *
*                   x(t[0]) in corn_1[0],    y(t[0]) in corn_n[1],     *
*                   ..                       ..                        *
*                   x(t[m-1]) in corn_n[0],  y(t[m-1]) in corn_n[1].   *
*          = 5: curvature radii r1 and rn at the end points, namely    *
*               r1 in corn_1[0] and rn in  corn_n[0]. The dirction of  *
*               the curvature is determined by the sign of r1 and rn.  *
*               For positive curvature radius, the curvature circle    *
*               lies on the left hand side of the curve as seen for    *
*               increasing papraters. (Curve concave from the left)    *
*               For a negative radius, the curve is concave from the   *
*               right. For our paramtrization, the paramters t[i]      *
*               increase with their indices i.                         *
*          = 6: 3rd  derivatives of the spline components wrt t at the *
*               end points are given :                                 *
*                   ...                      ...                       *
*                    x(t[0]) in corn_1[0],    y(t[0]) in corn_n[1],    *
*                   ...                      ...                       *
*                    x(t[m-1]) in corn_n[0],  y(t[m-1]) in corn_n[1].  *
* corn_1:  [0..1] left end point conditions                            *
* corn_n:  [0..1] right end point conditions                           *
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* cx: \   [0..m-2] vectors with spline coefficients for the x-component*
* dx:  \  sx(t) = ax[i] + bx[i] * (t - t[i]) + cx[i] * (t - t[i]) ^ 2  *
* ex:  /            + dx[i] * (t - t[i]) ^ 3 + ex[i] * (t - t[i]) ^ 4  *
* fx: /             + fx[i] * (t - t[i]) ^ 5                           *
*         ax is given by x, bx by xt.                                  *
*         cx and dx (as well as ax and bx) have one additional entry   *
*         cx[m-1] and dx[m-1].                                         *
* cy: \   [0..m-2] vectors for the spline coeeficients for y-component *
* dy:  \  sy(t) = ay[i] + by[i] * (t - t[i]) + cy[i] * (t - t[i]) ^ 2  *
* ey:  /            + dy[i] * (t - t[i]) ^ 3 + ey[i] * (t - t[i]) ^ 4  *
* fy: /             + fy[i] * (t - t[i]) ^ 5                           *
*         ay is given by  y, by  by yt.                                *
*         cy and dy (as well as ay and by) have one additional entry   *
*         cy[m-1] and dy[m-1].                                         *
* t:      [0..m-1] vector with values of curve paprmeter t             *
* xt: \   [0..m-1] normalized tangent vectors                          *
* yt: /                                                                *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* =  0: no error                                                       *
* =  1: invalid entry for marg                                         *
* =  2: invalid entry for kind of vectors in richt                     *
* =  3: m < 3                                                          *
* =  4: richt = 3: Tangent vectors cannot be reasonably assigned from  *
*       the derivatives, e.g. three successive points are collinear,   *
*       but the inner tangent is prescribed as normal to that line.    *
* =  5: two successive nodes are identical (double covered points are  *
*       otherwise allowed)                                             *
* =  6: One tangent or normal vector is zero.                          *
* =  7: y''(x) is prescribed at the endpoints (marg = 3), but the      *
*       curve has a vertical tangent there leading to mathematical     *
*       contradictions.                                                *
* =  8: one curvature radius is zero.                                  *
* =  9: A periodic spline is desired (marg = 1), but the data is not   *
* = 10: error in first call of  hermit()                               *
* = 11: error in second call of hermit()                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hermit, MAXROOT, sqr, REAL, FABS, SQRT, ZERO, ONE                    *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  n,                   /* m - 1, Index of last node              */
       marg_herm = 0,       /* type of end point condition            */
       fehler,              /* error code of hermit()                 */
       i;                   /* Loop variable                          */
  REAL deltx, delty,        /* aux variables for t[i]                 */
       delt,                /* inner product of corresponding tangent */
                            /* and chordal vectors                    */
       corn_1_herm = ZERO,  /* left end condition for hermit()        */
       corn_n_herm = ZERO;  /* right end condition for hermit()       */


  n = m - 1;


  if (marg < 1 || marg > 6)       /* invalid end point condition ?    */
    return 1;

  if (richt < 1 || richt > 3)     /* invalid directional vector type? */
    return 2;

  if (n < 2)                              /* too few nodes ?          */
    return 3;


  for (t[0] = ZERO, i = 1; i <= n; i++)   /* compute parameter        */
  {                                       /* vector t                 */
    deltx = x[i] - x[i - 1];
    delty = y[i] - y[i - 1];
    delt  = deltx * deltx +
            delty * delty;
    if (delt <= ZERO)              /* two equal nodes in succession ? */
      return 5;
    t[i] = t[i - 1] + SQRT(delt);
  }


  switch (richt)                 /* change  (xricht[i], yricht[i]) to */
  {                              /* tangent vector (xt[i], yt[i])     */
    case 1:                      /* Tangent vectors given ?           */
      for (i = 0; i <= n; i++)   /* copy to xt, yt                    */
        xt[i] = xricht[i],
        yt[i] = yricht[i];
      break;

    case 2:                            /* (xricht[i], yricht[i]) are  */
      for (i = 0; i <= n; i++)         /* normal vectors ?            */
        xt[i] = yricht[i],             /* rotate by 90 degrees and    */
        yt[i] = -xricht[i];            /* and store                   */
      break;

    case 3:                          /* derivative wrt x in  yricht ? */
      for (i = 0; i <= n; i++)       /* change to derivative wrt. t   */
                                     /*  in xt und yt                 */
        if (FABS(yricht[i]) >= MAXROOT)         /* derivative too     */
                                                /* large ?            */
          xt[i] = ZERO,                         /* record vertical    */
          yt[i] = ONE;                          /* tangent            */
        else
          xt[i] = ONE,
          yt[i] = yricht[i];
      break;

  }

  if (marg == 1)                              /* periodic spline ?    */
    if (y[0] != y[n] || y[0] != y[n] ||       /* Periodicity          */
        xt[0] != xt[n] || yt[0] != yt[n])     /* condition            */
      return 9;                               /* satisfied ?          */


  if (richt == 3)                       /* derivatives wrt x given ? */
  {

    /* If tangent and chordal vectors point in opposite directions,   */
    /* "we "turn' the tangents around.                                */
    /* If one tangent vector is perpendicular to both adjoining       */
    /* chordal vectors, the problem is not solvable: error return.    */

    if ((delt =                       /* angle between tangent and    */
         (x[1] - x[0]) * xt[0] +      /* next chordal vector > 90 deg?*/
         (y[1] - y[0]) * yt[0])
        < ZERO)
      xt[0] = -xt[0],                 /* turn tangent vector around   */
      yt[0] = -yt[0];
    else if (delt == ZERO)
      if (marg != 1)                  /* not a periodic spline ?      */
        return 4;
      else if ((delt =                       /* angle between tangent */
                (x[0] - x[n - 1]) * xt[0] +  /* and last chordal      */
                (y[0] - y[n - 1]) * yt[0])   /* vector > 90 degrees ? */
               < ZERO)
        xt[0] = -xt[0],                      /* turn around           */
        yt[0] = -yt[0];
      else if (delt == ZERO)
        return 4;

    for (i = 1; i < n; i++)         /* angle between tangent and next */
      if ((delt =                           /* chordal vector > 90 deg*/
           (x[i + 1] - x[i]) * xt[i] +      /* for all interior       */
           (y[i + 1] - y[i]) * yt[i])       /* nodes ?                */
          < ZERO)
        xt[i] = -xt[i],                     /* turn all around        */
        yt[i] = -yt[i];
      else if (delt == ZERO)
        if ((delt =                         /* angle > 90 deg again ? */
             (x[i] - x[i - 1]) * xt[i] +
             (y[i] - y[i - 1]) * yt[i])
            < ZERO)
          xt[i] = -xt[i],                   /* turn all around        */
          yt[i] = -yt[i];
        else if (delt == ZERO)
          return 4;

    if (marg == 1)         /* periodic spline ?                       */
                           /* final tangent vector must be precisely  */
      xt[n] = xt[0],       /* equal to the first tangent vector, or   */
      yt[n] = yt[0];       /* hermit() might hang up.                 */
    else if ((delt =
              (x[n] - x[n - 1]) * xt[n] +
              (y[n] - y[n - 1]) * yt[n])
             < ZERO)
      xt[n] = -xt[n],
      yt[n] = -yt[n];
    else if (delt == ZERO)
      return 4;
  }


  for (i = 0; i <= n; i++)                /* normalize tangent vectors*/
  {
    delt = SQRT(sqr(xt[i]) + sqr(yt[i]));
    if (delt <= ZERO)                     /* zero vector ?            */
      return 6;
    xt[i] /= delt;
    yt[i] /= delt;
  }


  switch (marg)           /* initialize for first call values for     */
                          /* hermit() :                               */
  {
    case 1:               /* type of boundary data in marg_herm,      */
      marg_herm = 1;      /* boundary data in in corn_1_herm, and     */
                          /*  corn_n_herm                             */
      corn_1_herm = ZERO; /* for  IBM C Set/2 1.0                     */
      corn_n_herm = ZERO;
      break;

    case 2:
      corn_1_herm = ZERO; /* for  IBM C Set/2 1.0                     */
      corn_n_herm = ZERO;
      marg_herm = 2;
      break;

    case 3:
      if (xt[0] == ZERO || xt[n] == ZERO)
        return 7;
      corn_1_herm = corn_n_herm = ONE;
      marg_herm = 3;
      break;

    case 4:
      corn_1_herm = corn_1[0];
      corn_n_herm = corn_n[0];
      marg_herm = 3;
      break;

    case 5:
      if (corn_1[0] == ZERO || corn_n[0] == ZERO)
        return 8;
      corn_1_herm = corn_n_herm = ONE;
      if (xt[0] == ZERO)
        corn_1_herm = -ONE / corn_1[0] / yt[0];
      if (xt[n] == ZERO)
        corn_n_herm = -ONE / corn_n[0] / yt[n];
      marg_herm = 3;
      break;

    case 6:
      corn_1_herm = corn_1[0];
      corn_n_herm = corn_n[0];
      marg_herm = 5;
  }


  /* --------------- compute spline component  SX ------------------- */

  fehler = hermit(n + 1, t, x, xt, marg_herm, corn_1_herm, corn_n_herm,
                  1, cx, dx, ex, fx);

  if (fehler != 0)                             /* error  in hermit()? */
     return 10;


  switch (marg)           /* initialize end point conditions for      */
  {                       /* second call of hermit()                  */
    case 3:
      corn_1_herm = (xt[0] * sqr(xt[0]) * corn_1[0] + yt[0]) / xt[0];
      corn_n_herm = (xt[n] * sqr(xt[n]) * corn_n[0] + yt[n]) / xt[n];
      break;

    case 4:
    case 6:
      corn_1_herm = corn_1[1];
      corn_n_herm = corn_n[1];
      break;

    case 5:
      corn_1_herm = corn_n_herm = ONE;
      if (xt[0] != ZERO)
        corn_1_herm = (ONE / corn_1[0] + yt[0]) / xt[0];
      if (xt[n] != ZERO)
        corn_n_herm = (ONE / corn_n[0] + yt[n]) / xt[n];

  }


  /* --------------- compute spline component  SY ------------------- */

  fehler = hermit(n + 1, t, y, yt, marg_herm, corn_1_herm, corn_n_herm,
                  0, cy, dy, ey, fy);

  if (fehler != 0)                             /* error in  hermit()? */
     return 11;


  return 0;
}

/* -------------------------- END hrmsplin.c ------------------------ */
