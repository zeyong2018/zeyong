#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE kubsplin.c ------------------------ */

/***********************************************************************
*                                                                      *
* Functions to compute cubic interpolating splines                     *
* -------------------------------------------------------------------  *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Source:               equivalent QuickBASIC program                  *
* Authors:              Elmar Pohl (QuickBASIC), Dorothee Seesing (C)  *
* Adaptation:           Juergen Dietel, Computer center, RWTH Aachen   *
* Date:                 2.2.1993                                       *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  sign, MAXROOT, sqr, NULL, SQRT, ONE, */
                        /*       FABS, PI, REAL, ACOS, ZERO, THREE,   */
                        /*       HALF, TWO                            */
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vmfree, vminit, */
                        /*       VEKTOR                               */
#include <u_proto.h>    /*  for  trdiag, tzdiag                       */
#include <kubsplin.h>   /*  for  spline, parspl, spltrans             */



/*.BA*/
/*.FE{C 10.1.2}
     {Computation of Non-Parametric Cubic Splines}
     {Computation of Non-Parametric Cubic Splines}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int spline      /* non-parametric cubic splines ......................*/
/*.IX{spline}*/
          (
           int  m,            /* number of nodes .....................*/
           REAL x[],          /* nodes : x-values ....................*/
           REAL y[],          /*         y-values ....................*/
           int  marg_cond,    /* type of end point condition .........*/
           REAL marg_0,       /* left end point condition ............*/
           REAL marg_n,       /* right end point condition ...........*/
           int  save,         /* save aux vectors ? ..................*/
           REAL b[],          /* Spline coefficients of (x-x[i]) .....*/
           REAL c[],          /* Spline coefficients of (x-x[i])^2 ...*/
           REAL d[]           /* Spline coefficients of (x-x[i])^3 ...*/
          )                   /* error code ..........................*/

/***********************************************************************
* We compute the coefficients of the non-parametric cubic              *
* interpolating spline for the given nodes :                           *
*                    (x[i], y[i]), i = 0,...,m-1.                      *
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
* x:          [0..m-1] vector of nodes: x-values (not used if last call*
*             was executed with save != 0)                             *
* y:          [0..m-1] vector, y-values of the nodes                   *
* marg_cond:  type of boundary condition:                              *
*             = 0: not-a-node condition (marg_0, mag_n not used)       *
*             = 1: first end point derivatives                         *
*             = 2: second end point derivatives (if marg_0 = marg_n = 0*
*                  this specifies a natural spline)                    *
*             = 3: 3nd derivatives at end points prescribed            *
*             = 4: periodic spline (marg_0, marg_n not used)           *
* marg_0:     boundary condition at x[0]                               *
* marg_n:     boundary condition at x[m-1]                             *
* save:       Flag, which indicates whether the LU decomposition can be*
*             used for a subsequent call, provided the x-values stay   *
*             the same. Usually save = 0. If save = 1 for the same set *
*             of x-values one can save about 4*m flops by not recom-   *
*             puting the factorization for the tridiagonal system      *
*             matrix.                                                  *
*             The last call of a run, however, must be executed with   *
*             save = 0, otherwise certain storage is erroneaously kept *
*             and will foul up future runs of spline.                  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* b: \   [0..m-2] or [0..m-1] vectors of the spline coefficients:      *
* c:  >  s(x)  =   a[i] + b[i] * (x - x[i]) +                          *
* d: /           + c[i] * (x - x[i]) ^ 2 + d[i] * (x - x[i]) ^ 3.      *
*        Here a denotes the vector y.                                  *
*        c has m entries (as does a).                                  *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: no error                                                        *
* =-1: monotonicity violated: x[i-1] >= x[i]                           *
* = 1: wrong value for marg_cond                                       *
* = 2: m < 3                                                           *
* = 3: lack of storage                                                 *
* = 4: marg_cond =4; input non periodic                                *
* > 4: error in trdiag() or tzdiag()                                   *
* In case of error, the output is indeterminant and the aux storage is *
* freed.                                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, vminit, vmalloc, vmcomplete, vmfree, VEKTOR, trdiag, tzdiag,   *
* NULL, ZERO, THREE, HALF, TWO                                         *
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
         fehler,             /* error code of trdiag() or tzdiag()    */
         i,                  /* Loop variable                         */
         erster_aufruf;      /* Flag, keeping track of the first call */
                             /* of a run for the same x-values        */

  n = m - 1;


  if (n < 2)                                           /* n too small */
    ciao(2);

  if (marg_cond < 0 || marg_cond > 4)  /* wrong  marg_cond value      */
    ciao(1);

  if (marg_cond == 4)                   /* periodic spline?           */
    if (y[n] != y[0])                   /* check periodicity          */
      ciao(4);


  /* 1st call: allocate storage: four [0..n-1] vectors                */
  /* (for the periodic case an additional two [0..n-3] vectors)       */

  if (vmblock == NULL)                  /* first call of one run ?    */
  {
    erster_aufruf = 1;
    #define MYALLOC(l)  (REAL *)vmalloc(vmblock, VEKTOR, (l), 0)
    vmblock = vminit();                           /* allocate storage */
    h     = MYALLOC(n);
    lower = MYALLOC(n);
    diag  = MYALLOC(n);
    upper = MYALLOC(n);
    if (marg_cond == 4)                   /* periodic spline for      */
      if (n > 2)                          /* sufficiently many nodes? */
        lowrow = MYALLOC(n - 2),          /* allocate additional      */
        ricol  = MYALLOC(n - 2);          /* storage                  */
    #undef MYALLOC
  }
  else
    erster_aufruf = 0;
  if (! vmcomplete(vmblock))                      /* lack of memory ? */
    ciao(3);


  if (erster_aufruf)
    for (i = 0; i < n; i++)  /* compute step sizes and check          */
    {                        /*  monotonicity                         */
      h[i] = x[i + 1] - x[i];
      if (h[i] <= ZERO)
        ciao(-(i + 1));
    }


  for (i = 0; i < n - 1; i++)                   /* form linear system */
  {
    c[i] = THREE * ((y[i + 2] - y[i + 1]) / h[i + 1]  /* right side   */
                  - (y[i + 1] - y[i])     / h[i]);
    if (erster_aufruf)
      diag[i] = TWO * (h[i] + h[i + 1]),        /* form main diagonal */
      lower[i + 1] = upper[i] = h[i + 1];       /* and co-diagonals   */
  }

  switch (marg_cond)    /* correct coefficients according to end point*/
  {                     /* conditiones                                */
    case 0:                     /* not-a-node-condition ?             */
      if (n == 2)                     /* only three nodes ?           */
      {                    /* In which case the system is under-      */
                           /* determined and we only compute a        */
                           /* quadratic polynomial.                   */
        c[0] /= THREE;                             /* right side      */
        if (erster_aufruf)
          diag[0] *= HALF;                         /* adjust matrix   */
      }
      else                            /* more than three nodes ?      */
      {
        c[0]     *= h[1]     / (h[0]     + h[1]);           /* right  */
        c[n - 2] *= h[n - 2] / (h[n - 1] + h[n - 2]);       /* side   */
        if (erster_aufruf)
          diag[0]      -= h[0],                           /* adjust   */
          diag[n - 2]  -= h[n - 1],                       /* Matrix   */
          upper[0]     -= h[0],
          lower[n - 2] -= h[n - 1];
      }
      break;

    case 1:                     /* first derivative at end points ?  */
      c[0]     -= (REAL)1.5 * ((y[1] - y[0]) / h[0] - marg_0);
      c[n - 2] -= (REAL)1.5 * (marg_n - (y[n] - y[n - 1]) / h[n - 1]);
      if (erster_aufruf)
        diag[0]     -= HALF * h[0],                /* adjust matrix   */
        diag[n - 2] -= HALF * h[n - 1];
      break;

    case 2:                     /* second derivative at end points ?  */
      c[0]     -= h[0]     * HALF * marg_0;
      c[n - 2] -= h[n - 1] * HALF * marg_n;
      break;

    case 3:                     /* third derivative at end point ?    */
      c[0]        += HALF * marg_0 * h[0]     * h[0];
      c[n - 2]    -= HALF * marg_n * h[n - 1] * h[n - 1];
      if (erster_aufruf)
        diag[0]     += h[0],                       /* initialize      */
        diag[n - 2] += h[n - 1];                   /* matrix          */
      break;

    case 4:                     /* periodic spline?                   */
      c[n - 1] = THREE * ((y[1] - y[0])     / h[0] -
                          (y[n] - y[n - 1]) / h[n - 1]);
      if (erster_aufruf)
        if (n > 2)
          diag[n - 1]  = TWO * (h[0] + h[n - 1]),
          ricol[0] = lowrow[0] = h[0];
  }


  switch (n)                /* solve linear system for the spline     */
  {                         /* coefficients c[i]                      */
    case 2:                           /* only threee nodes =>         */
                                      /* => solve directly            */
      if (marg_cond == 4)             /* periodic spline ?            */
        c[1] = THREE * (y[0] - y[1]) / (x[2] - x[1]) / (x[1] - x[0]),
        c[2] = - c[1];
      else
        c[1] = c[0] / diag[0];
      break;

    default:                          /* more than three nodes ?      */
      if (marg_cond == 4)             /* periodic spline ?            */
        fehler = tzdiag(n, lower, diag, upper, lowrow,
                        ricol, c, !erster_aufruf);
      else                               /* non periodic spline ?     */
        fehler = trdiag(n - 1, lower, diag, upper, c, !erster_aufruf);

      if (fehler != 0)        /* error in  tzdiag() or in trdiag() ?  */
        ciao(fehler + 4);
      for (i = n; i != 0; i--)  /* shift solution vector down         */
        c[i] = c[i - 1];
  }


  switch (marg_cond)    /* correct first or last entry of c depending */
  {                     /* on end point condition                     */
    case 0:                     /* not-a-node condition ?             */
      if (n == 2)                     /* only three nodes ?           */
        c[0] = c[2] = c[1];
      else                            /* more than three nodes ?      */
        c[0] = c[1] + h[0] * (c[1] - c[2]) / h[1],
        c[n] = c[n - 1] + h[n - 1] *
               (c[n - 1] - c[n - 2]) / h[n - 2];
      break;

    case 1:                     /* first end point derivative given ? */
      c[0] =  (REAL)1.5 * ((y[1] - y[0]) / h[0] - marg_0);
      c[0] = (c[0] - c[1] * h[0] * HALF) / h[0];
      c[n] = (REAL)-1.5 * ((y[n] - y[n - 1]) / h[n - 1] - marg_n);
      c[n] = (c[n] - c[n - 1] * h[n - 1] * HALF) / h[n - 1];
      break;

    case 2:                     /* second derivative ?                */
      c[0] = marg_0 * HALF;
      c[n] = marg_n * HALF;
      break;

    case 3:                     /* third derivative ?                 */
      c[0] = c[1]     - marg_0 * HALF * h[0];
      c[n] = c[n - 1] + marg_n * HALF * h[n - 1];
      break;

    case 4:                     /* periodic spline?                   */
      c[0] = c[n];

  }


  for (i = 0; i < n; i++)                      /* compute b[i], c[i]  */
    b[i] = (y[i + 1] - y[i]) / h[i] - h[i] *
           (c[i + 1] + TWO * c[i]) / THREE,
    d[i] = (c[i + 1] - c[i]) / (THREE * h[i]);


  if (!save)                         /* save aux LU factors ?         */
    ciao(0);                         /* (last call of a run) ?        */

  return 0;
#undef ciao
}
/*.BA*/



/*.FE{C 10.1.3}
     {Computing Parametric Cubic Splines}
     {Computing Parametric Cubic Splines}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int parspl           /* parametric cubic splines .....................*/
/*.IX{parspl}*/
          (
           int  m,            /* number of nodes .....................*/
           REAL x[],          /* nodes : x-values ....................*/
           REAL y[],          /*         y-values ....................*/
           int  marg_cond,    /* type of end point condition .........*/
           REAL marg_0[],     /* left end point condition ............*/
           REAL marg_n[],     /* right end point condition ...........*/
           int  cond_t,       /* Parameter nodes given ? .............*/
           REAL t[],          /* Parameter nodes .....................*/
           REAL bx[],         /* x spline coeffic. for (t-t[i]) ......*/
           REAL cx[],         /* x spline coeffic. for (t-t[i])^2 ....*/
           REAL dx[],         /* x spline coeffic. for (t-t[i])^3 ....*/
           REAL by[],         /* y spline coeffic. for (t-t[i]) ......*/
           REAL cy[],         /* y spline coeffic. for (t-t[i])^2 ....*/
           REAL dy[]          /* y spline coeffic. for (t-t[i])^3 ....*/
          )                   /* error code ..........................*/

/***********************************************************************
* Compute the coefficients of an interpolating parametric cubic spline *
* for the given nodes                                                  *
*                      (x[i], y[i]), i = 0,...,m-1.                    *
* The type of end point condition can be prescribed via marg_cond.     *
* The parameter nodes can either be given or they can be computed      *
* internally.                                                          *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* m:         number of nodes (m > 2)                                   *
* x: \       [0..m-1] vectors with node data                           *
* y: /                                                                 *
* marg_cond: type of end point conditions :                            *
*            = 0: not-a-node condition                                 *
*            = 1: first derivative wrt. t given :                      *
*                 .                      .                             *
*                 sx(t[0]) in marg_0[0], sy(t[0]) in marg_0[1]),       *
*                 .                        .                           *
*                 sx(t[m-1]) in marg_n[0], sy(t[m-1]) in marg_n[1])    *
*            = 2: second derivative prescribed at ends (natural spline *
*                 for  marg_cond = 2  with marg_0 = marg_n =0) :       *
*                 ..                     ..                            *
*                 sx(t[0]) in marg_0[0], sy(t[0]) in marg_0[1]),       *
*                 ..                       ..                          *
*                 sx(t[m-1]) in marg_n[0], sy(t[m-1]) in marg_n[1])    *
*            = 3: periodic spline                                      *
*            = 4: end point derivatives dy/dx given :                  *
*                 y'(x[0]) in marg_0[0] and y'(x[m-1]) in marg_n[0].   *
*                 (marg_0[1] and marg_n[1] are not used here.)         *
* marg_0:    [0..1] vector of end point conditions for t[0] if         *
*            marg_cond = 1, 2, 4 (not used if marg_cond = 0 , 3)       *
* marg_n:    [0..1] vector of end point conditions for t[m-1]          *
* cond_t:    Control for parameter values t[i]:                        *
*            =  0: Parameter values t[i] are computed internally.      *
*            != 0: parameter values are prescribed by user.            *
* t          [0..m-1] vector of monotonically increasing parameter     *
*            values if cond_t != 0.                                    *
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* t:     [0..m-1] vector of computed parameter values t[i] for         *
*        cond_t = 0                                                    *
* bx: \  [0..m-2] vectors with spline coefficients for x:              *
* cx:  >     sx(t)  =  ax[i] + bx[i] * (t - t[i]) +  cx[i]             *
* dx: /                * (t - t[i]) ^ 2  + dx[i] * (t - t[i]) ^ 3.     *
*        ax denotes x.                                                 *
*        cx as well as ax have one additional entry in position [m-1]. *
* by: \  [0..m-2] vectors with the spline coefficients for y :         *
* cy:  >     sy(t)  =  ay[i] + by[i] * (t - t[i]) +  cy[i]             *
* dy: /                * (t - t[i]) ^ 2  + dy[i] * (t - t[i]) ^ 3.     *
*        ay denotes y.                                                 *
*        cy as well as ay have one additional entry at [m-1].          *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* =  0: no error                                                       *
* = -i: two successive nodes are identical :                           *
*       (x[i-1], y[i-1]) = (x[i], y[i]).                               *
* =  1: m < 3                                                          *
* =  2: wrong type in  marg_cond                                       *
* =  3: periodic spline and x[0] != x[m-1]                             *
* =  4: periodic spline and y[0] != y[m-1]                             *
* =  5: the t[i] are not monotonic.                                    *
* >  5: error in  spline()                                             *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* spline, REAL, sign, MAXROOT, sqr, FABS, SQRT, ONE, ZERO              *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i,                /* Loop variable                             */
       n,                /* m - 1, Index of final node                */
       mess = 0,         /* type of end point condition for spline()  */
       fehler;           /* error code from  spline() or  perspl()    */
  REAL deltx,            /* aux variables for t[i]                    */
       delty,
       delt,
       alfx = ZERO,      /* left end condition of sx for spline()     */
       betx = ZERO,      /* right end condition of sx for spline()    */
       alfy = ZERO,      /* left end condition of sy for spline()     */
       bety = ZERO;      /* right end condition of sy for spline()    */


  n = m - 1;


  if (n < 2)                        /* less than three nodes ?        */
    return 1;
  if (marg_cond < 0 || marg_cond > 4)          /* wrong marg_cond ?   */
    return 2;


  /* If t is not prescribed, we compute their values now :            */
  /* We shall use the chordal parametrization, which will approximate */
  /* the arc length parametrization roughly.                          */

  if (cond_t == 0)            /* Parameter values not prescribed ?    */
    for (t[0] = ZERO, i = 1; i <= n; i++)
    {
      deltx = x[i] - x[i - 1];
      delty = y[i] - y[i - 1];
      delt  = deltx * deltx + delty * delty;
      if (delt <= ZERO)
        return -i;
      t[i] = t[i - 1] + SQRT(delt);
    }


  switch (marg_cond)   /* initialize mess, alfx, betx, alfy, bety for */
  {                    /* spline() depending on end point condition   */
    case 0:                                     /* not-a-node spline? */
      mess = 0;
      break;
    case 1:                  /* first end derivatives wrt t given ?   */
    case 2:                  /* second end derivatives ?              */
      mess = marg_cond;
      alfx = marg_0[0];
      alfy = marg_0[1];
      betx = marg_n[0];
      bety = marg_n[1];
      break;
    case 3:                                   /* periodic spline ?    */
      mess = 4;
      if (x[n] != x[0])              /* non periodic data for a       */
        return 3;
      if (y[n] != y[0])              /* periodic spline ?             */
        return 4;
      alfx = betx = ZERO;    /* as a precaution for  IBM C Set/2 1.0  */
      alfy = bety = ZERO;
      break;
    case 4:                   /* first derivative dy/dx given ?       */
      mess = 1;         /* for spline(): first derivatives prescribed */
      if (FABS(marg_0[0]) >= MAXROOT)      /* derivatives excessive ? */
        alfx = ZERO,                       /* use vertical tangent    */
        alfy = sign(ONE, y[1] - y[0]);
      else
        alfx = sign(SQRT(ONE / (ONE + sqr(marg_0[0]))), x[1] - x[0]),
        alfy = alfx * marg_0[0];
      if (FABS(marg_n[0]) >= MAXROOT)
        betx = ZERO,
        bety = sign(ONE, y[n] - y[n - 1]);
      else
        betx = sign(SQRT(ONE / (ONE + sqr(marg_n[0]))),x[n] - x[n - 1]),
        bety = betx * marg_n[0];
  }


  fehler = spline(n + 1, t, x, mess, alfx, betx, 1, bx, cx, dx);
  if (fehler < 0)
    return 5;
  else if (fehler > 0)
    return fehler + 5;

  fehler = spline(n + 1, t, y, mess, alfy, bety, 0, by, cy, dy);
  if (fehler != 0)
    return fehler + 20;


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int spltrans     /* transformed parametric cubic splines .............*/
/*.IX{spltrans}*/
            (
             int  m,         /* number of nodes ......................*/
             REAL x[],       /* nodes : x-values .....................*/
             REAL y[],       /*         y-values .....................*/
             int  mv,        /* type of transformation of origin .....*/
             REAL px[],      /* coordinats of the transformation      */
             REAL py[],      /* vector  P ............................*/
             REAL a[],       /* Spline coeff. of (phi-phin[i])^0 .....*/
             REAL b[],       /* Spline coeff. of (phi-phin[i]) .......*/
             REAL c[],       /* Spline coeff. of (phi-phin[i])^2 .....*/
             REAL d[],       /* Spline coeff. of (phi-phin[i])^3 .....*/
             REAL phin[],    /* angular coordinates of nodes .........*/
             REAL *phid      /* angle of rotation of coordinate system*/
            )                /* error code ...........................*/

/***********************************************************************
*                                                                      *
* Compute the coefficients of a transformed parametric interpolating   *
* cubic spline for a closed smooth curve.                              *
*                                                                      *
* A transformed parametric cubic spline is a periodic cubic spline     *
* as in  spline(), but in polar representation.                        *
* This method often enables us to interpolate data, whose nodes are not*
* monotonic without  having to compute proper parametric splines as    *
* in parspl().                                                         *
.BE*)
* First we transform the given nodes to polar coordinates              *
* (phin[i],a[i]), where  phin[i] denotes the angle and  a[i] denotes   *
* the length of the node vector (x[i],y[i]). This must be done in such *
* a way that the angles phin[i] are monotonically increasing.          *
* Otherwise this method does not work.                                 *
* Then a periodic cubic interpolating spline is computed for the nodes *
* (phin[i], a[i]).                                                     *
* In order to achieve monotonicity of the  phin[i] it may be necessary *
* to translate the origin to P = (px, py) and to rotate the coordinate *
* system by an angle  phid . (px, py) must lie inside the convex hull  *
* of (x[i], y[i]), so that every polar ray from P intersects its edge  *
* presicely once.                                                      *
* P may be given by the user or generated internally.                  *
* The computed value for P is, however, only an approximation, which   *
* may not always satisfy the above conditions.                         *
* Moreover the  (x[i],y[i]) must be ordered so that they lie in        *
* counter clockwise direction for increasing indices. And x[m-1] = x[0]*
* and  y[m-1] = y[0]  must also hold.                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* m:    number of nodes  (m > 2)                                       *
* x:    [0..m-1] node vectors                                          *
* y:                                                                   *
* mv:   Control for origin shifts :                                    *
*       mv > 0: user supplies coordinates for px, py                   *
*       mv = 0: no transformation  (i.e., px = py = 0)                 *
*       mv < 0: px and py shall be computed.                           *
*               We set :                                               *
*               px = (xmax + xmin) / 2                                 *
*               py = (ymax + ymin) / 2                                 *
*               where xmax = max(x[i]), xmin = min(x[i]),              *
*                     ymax = max(y[i]), ymin = min(y[i]), i=0,...,m-1. *
*               NOTE: This is a rough guess only.                      *
*                     For return value = -3, the user ought to give a  *
*                     better guess from a plot of the data, maybe.     *
*                                                                      *
* px: \ for mv > 0 :  coordinates of P                                 *
* py: /                                                                *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* a: \  [0..m-1] vectors with the spline coefficients :                *
* b:  \     S(phi)  =  a[i] + b[i] * (phi - phin[i])                   *
* c:  /                     + c[i] * (phi - phin[i]) ^ 2               *
* d: /                      + d[i] * (phi - phin[i]) ^ 3               *
*       for  phin[i] <= phi <= phin[i+1],   i=0,...,m-2.               *
*       the a[i] are the length of the vectors (x[i],y[i]).            *
* phin: [0..m-1] angle vector for the nodes (x[i],y[i])                *
*       we have  phin[0]   = 0,                                        *
*                phin[i]   = arctan((y[i] - py) / (x[i] - px)) - phid, *
*                            i=1, ..., m-2                             *
*                phin[m-1] = 2 * Pi                                    *
*                                                                      *
* px: \ Coordinates of transformation vector  P                        *
* py: /                                                                *
* phid: angle of rotation for coordinate system;                       *
*              phid = arctan(y[0] / x[0]).                             *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
*  0: no error                                                         *
* -1: m < 3                                                            *
* -3: phin[i] not monotonically increasing.                            *
* -4: x[m-1] != x[0]   or   y[m-1] != y[0]                             *
* >0: error in  spline()                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* spline, REAL, PI, sqr, SQRT, ACOS, ZERO, TWO                         *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL xmin,                 /* Minimum of x[i]                       */
       xmax,                 /* Maximum of x[i]                       */
       ymin,                 /* Minimum of y[i]                       */
       ymax,                 /* Maximum of y[i]                       */
       sa,                   /* sin(-phid)                            */
       ca;                   /* cos(-phid)                            */
  int  n,                    /* m - 1, Index of final node            */
       i;                    /* Loop variable                         */


  n = m - 1;


  /* ---------------- check input ----------------------------------- */
  if (n < 2)
    return -1;
  if (x[0] != x[n] || y[0] != y[n])
    return -4;


  /* ---------------- tranform -------------------------------------- */
  if (mv == 0)                                 /* no transformation ? */
  {
    *px = *py = ZERO;
    for (i = 0; i <= n; i++)
      b[i] = x[i],
      c[i] = y[i];
  }
  else                               /* transform with P = (px, py) ? */
  {
    if (mv < 0)                               /* Compute  py and py ? */
    {
      xmax = x[0];
      xmin = x[0];
      ymax = y[0];
      ymin = y[0];
      for (i = 1; i <= n; i++)
      {
        if (x[i] > xmax)
          xmax = x[i];
        if (x[i] < xmin)
          xmin = x[i];
        if (y[i] > ymax)
          ymax = y[i];
        if (y[i] < ymin)
          ymin = y[i];
      }
      *px = (xmax + xmin) / TWO;
      *py = (ymax + ymin) / TWO;
    }

    for (i = 0; i <= n; i++)     /* store transformed coordinates     */
      b[i] = x[i] - *px,         /* in (b[i],c[i])                    */
      c[i] = y[i] - *py;
  }


  /* ---- compute transformed nodes :                 --------------- */
  /* ---- first compute norms a[i]. Stop if  a[i] = 0 --------------- */

  for (i = 0; i <= n; i++)
  {
    a[i] = SQRT(sqr(b[i]) + sqr(c[i]));
    if (a[i] == ZERO)
      return -3;
  }


  /*------------------------------------------------------------------*/
  /* Secondly, rotate the coordinates by alpha                        */
  /*                                                                  */
  /*  (X1)   ( cos(alpha)   -sin(alpha) ) (X)                         */
  /*  (  ) = (                          ) ( )                         */
  /*  (Y1)   ( sin(alpha)    cos(alpha) ) (Y)                         */
  /*                                                                  */
  /*  for alpha = -phid                                               */
  /*------------------------------------------------------------------*/

  *phid = ACOS(b[0] / a[0]);
  if (c[0] < ZERO)
    *phid = TWO * PI - *phid;
  ca = b[0] / a[0];
  sa = -c[0] / a[0];
  for (i = 0; i <= n; i++)               /* store rotated coordinates */
    d[i] = b[i] * ca - c[i] * sa,        /*   in  (d[i],c[i]          */
    c[i] = b[i] * sa + c[i] * ca;


  /* ------ compute angles  phin[i]. Stop if angles do not ---------- */
  /* ------ increase monotonically                         ---------- */

  phin[0] = ZERO;
  for (i = 1; i < n; i++)
  {
    phin[i] = ACOS(d[i] / a[i]);
    if (c[i] < ZERO)
      phin[i] = TWO * PI - phin[i];
    if (phin[i] <= phin[i - 1])
      return -3;
  }
  phin[n] = TWO * PI;


  /* --------------- compute spline coefficients -------------------- */
  return spline(n + 1, phin, a, 4, ZERO, ZERO, 0, b, c, d);
}

/* ------------------------- END kubsplin.c ------------------------- */
