#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE subsplin.c ----------------------- */

/***********************************************************************
*                                                                      *
* Module to compute  Akima and Renner Subsplines                       *
* ----------------------------------------------                       *
*                                                                      *
* exported functions:                                                  *
*   - akima():     compute coefficients of an Akima subsline           *
*   - renner():    compute coefficients of a Renner subspline          *
*   - rennwert():  compute functions value and derivatives of a        *
*                  Renner-Subsplines at a certain place                *
*   - renntab():   make a table with function values of a              *
*                  Renner subspline                                    *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 with 80387                         *
* Source:               equivalenten TP-Unit, with partially new codes *
* Author:               Elmar Pohl (QuickBASIC)                        *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Date:                 5.25.1994                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  REAL, ZERO, ONE, FABS, MACH_EPS, min, */
                       /*       sqr, TWO, THREE, copy_vector, SQRT,   */
                       /*       norm_max, SIX, intervall, skalprod    */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR, MATRIX                        */
#include <subsplin.h>  /*  for  akima, renner, rennwert, renntab      */



/*.BA*/
/*.FE{C 13.1}{Akima Subsplines}{Akima Subsplines}*/

/*.BE*/
/* ------------------------------------------------------------------ */

static int a_eckrund      /* round corners for Akima subsplines ......*/
/*.IX{a\unt eckrund}*/
        (
         int  *n,                  /* Number of final node ...........*/
         int  nmax,                /* upper index limit of  x and y ..*/
         REAL x[],                 /* nodes: x-values ................*/
         REAL y[],                 /*        y-values ................*/
         REAL beta,                /* rounding parameter ............ */
         int  perio                /* periodic interpolation? ........*/
        )                          /* error code .....................*/

/***********************************************************************
* This aux function for akima() rounds corners for Akima subsplines    *
* by replacing each corner with two slightly transformed nodes.        *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n     Number of final node. Numbering starts at zero  (n > 3)        *
* nmax  upper index limit on  x, y (for checking whether rounding can  *
*       proceed in the allotted space)                                 *
* x     [0..n] node vectors, monotonic for x-values.                   *
* y     [0..n] node vector for y-values.                               *
*       The vectors  x and x must be large enough to accommodate       *
*       maximally [(n+1)/2] extra points from rounding.                *
* beta  determines amount of rounding. Each corner point is replaced   *
*       by two points obtained by moving the corner point towards its  *
*       left and right neighbor a short distance that depends on beta. *
*       If beta lies outside the interval (0,1) no rounding takes      *
*       place.                                                         *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* n     new number of nodes, likely larger than its input value        *
* x,y   node vectors as on input, except for their number              *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error (or beta not in (0,1))                                   *
* 2: x[i] not monotonically increasing                                 *
* 5: not enough storage for rounding all corners                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, ONE, FABS, MACH_EPS, min                                 *
***********************************************************************/

{
  REAL qimin2,         /* Secant slope at P[i-2]                      */
       qimin1,         /* Secant slope at P[i-1]                      */
       qi,             /* Secant slope at P[i]                        */
       qiplus1,        /* Secant slope at P[i+1]                      */
       L, R, B,        /* aux variable for weights                    */
       lambda,         /* weights for shifting backwards              */
       my;             /* weights for shifting forwards               */
  int  i, j,           /* Loop variables                              */
       anf,            /* point index where check for corners begins  */
       end,            /* point index where check for corners ends    */
       anfm2,          /* `(anf - 2) modulo n' (Periodizitaet!)       */
       ip2;            /* `(i + 2) modulo n' (Periodizitaet!)         */


  if (beta <= ZERO || beta >= ONE)          /* rounding not desired ? */
    return 0;

  if (perio)
    anf   = 1,
    end   = *n - 1,
    anfm2 = anf - 2 + *n;
  else
    anf   = 2,
    end   = *n - 2,
    anfm2 = anf - 2;

  if (x[anf - 1] == x[anfm2]   ||         /* x-values not monotone?   */
      x[anf]     == x[anf - 1] ||
      x[anf + 1] == x[anf])
    return 2;


  qimin2 = (y[anf - 1] - y[anfm2])   / (x[anf - 1] - x[anfm2]);
  qimin1 = (y[anf]     - y[anf - 1]) / (x[anf]     - x[anf - 1]);
  qi     = (y[anf + 1] - y[anf])     / (x[anf + 1] - x[anf]);


  for (i = anf; i <= end; i++)
  {
    ip2 = (i == *n - 1) ? 1 : (i + 2);
    if (x[ip2] == x[i + 1])               /* not strictly monotone?   */
      return 2;
    qiplus1 = (y[ip2] - y[i + 1]) / (x[ip2] - x[i + 1]);

    if (FABS(qiplus1 - qi) +                            /* found a    */
        FABS(qimin1  - qimin2) <  MACH_EPS &&           /* corner?    */
        FABS(qi      - qimin1) >= MACH_EPS)
    {

      if (*n == nmax)                     /* no room for new points?  */
        return 5;

      for (j = *n; j >= i; j--)           /* move labels for points  */
        x[j + 1] = x[j],                  /* P[i] ... P[n] up by one */
        y[j + 1] = y[j];

      /* ---------- move P[i] and P[i+1] a little bit -------------- */

      L        =  (x[i]     - x[i - 1]);
      R        =  (x[i + 2] - x[i + 1]);
      B        =  beta * min(L, R);
      lambda   =  B / L;
      my       =  B / R;
      x[i]     -= lambda * (x[i]     - x[i - 1]);
      y[i]     -= lambda * (y[i]     - y[i - 1]);
      x[i + 1] += my     * (x[i + 2] - x[i + 1]);
      y[i + 1] += my     * (y[i + 2] - y[i + 1]);

      (*n)++;                          /* there is one more node now  */
      end++;

      qimin2 = (y[i]     - y[i - 1]) / (x[i]     - x[i - 1]);
      qimin1 = (y[i + 1] - y[i])     / (x[i + 1] - x[i]);
      qi     = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
    }

    else                                            /* not a corner?  */
      qimin2 = qimin1,
      qimin1 = qi,
      qi     = qiplus1;
  }


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int akima                 /* compute coeffic. of an Akima subspline ..*/
/*.IX{akima}*/
        (
         int  *n,                  /* Number of final node ...........*/
         int  nmax,                /* upper index limit for nodes x, y*/
         REAL x[],                 /* nodes: x-values ................*/
         REAL y[],                 /*        y-values ............... */
         int  perio,               /* periodic interpolation? ........*/
         REAL beta,                /* rounding parameter .............*/
         REAL b[],                 /* Spline coefficients ............*/
         REAL c[],
         REAL d[]
        )                          /* error code                      */

/***********************************************************************
* Compute the coefficients for an interpolating Akima subspline.       *
.BE*)
*                                                                      *
* The Akima subspline funkcion can be evaluates analogously to cubic   *
* splines via  spwert() from the module  spliwert. It can be tabulated *
* with  sptab() from  splintab as well.                                *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n      Number of final node. Numbering starts at zero. (n > 3)       *
* nmax   upper index limit on  x, y, b, c, d (for checking whether     *
*        rounding can proceed in the allotted space)                   *
* x      [0..n] node vectors, monotonic for x-values.                  *
* y      [0..n] node vector for y-values.                              *
*        For periodic interpolation y[0] = y[n].                       *
*        The vectors  x and x must be large enough to accommodate      *
*        maximally [(n+1)/2] extra points from rounding.               *
* perio  Flag for periodic or non periodic Interpolation:              *
*        perio FALSE: non periodic interpolation                       *
*        perio TRUE : periodic interpolation                           *
* beta   determines amount of rounding. Each corner point is replaced  *
*        by two points obtained by moving the corner point towards its *
*        left and right neighbor a short distance that depends on beta.*
*        If beta lies outside the interval (0,1) no rounding takes     *
*        place.                                                        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* b \   [0..n-1] spline coefficient vectors forming :                  *
* c  >     s(x)  =  a[i] + b[i] * (x - x[i]) + c[i] * (x - x[i]) ^ 2   *
* d /                    + d[i] * (x - x[i]) ^ 3                       *
*       a corresponds to  y and has one additional entry  a[n].        *
*                                                                      *
* And only for beta in (0,1):                                          *
* n     new value of n with injected rounding nodes                    *
* x,y   extended node list                                             *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: n < 4                                                             *
* 2: x[i] not monotonoc                                                *
* 3: while perio = TRUE:  y[0] != y[n]                                 *
* 4: lack of memory                                                    *
* 5: no more space for rounded corners                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, a_eckrund, MACH_EPS, sqr, FABS, vminit, vmalloc, vmcomplete,   *
* vmfree, VEKTOR, ZERO, TWO, THREE                                     *
.BA*)
***********************************************************************/
/*.BE*/

{
  void *vmblock;   /* List of dynamic allocations                     */
  REAL *tL,        /* [0..n] vector of left side slopes               */
                   /* (tR is not needed here because its values can   */
                   /* be stored in b.)                                */
       *m,         /* [-2..n+1] vector with secant slopes             */
       hi,         /* length of current node interval                 */
       L,          /* |m[i-1] - m[i-2]|                               */
       R,          /* |m[i+1] - m[i]|                                 */
       nenner,     /* denominator for alpha                           */
       alpha;      /* Factor for computing left and right sided slopes*/
  int  i,          /* Loop variable                                   */
       fehler;     /* error code from a_eckrund()                     */


  if (*n < 4)
    return 1;

  if (perio)                         /* periodic Akima interpolation? */
    if (y[*n] != y[0])               /* invalid y values?             */
      return 3;                      /* report error                  */


  fehler = a_eckrund(n, nmax, x, y, beta, perio);  /* round corners   */
  if (fehler != 0)                                 /* if desired      */
    return fehler;


  /* ----------- allocate buffers for aux vectors :         --------- */
  /* ----------- one [0..n+3] vector and one [0..n] vector  --------- */

  vmblock = vminit();                          /* initialize storage  */
  tL = (REAL *)vmalloc(vmblock, VEKTOR, *n + 1, 0);
  m  = (REAL *)vmalloc(vmblock, VEKTOR, *n + 4, 0);
  if (! vmcomplete(vmblock))                      /* lack of memory ? */
    return 4;

  m += 2;                         /* allows for negative indices in m */


  for (i = 0; i < *n; i++)          /* compute slopes  m[0]..m[n-1]   */
  {                                 /* while testing for monotonicity */
    hi = x[i + 1] - x[i];
    if (hi <= ZERO)                 /* not montonic ?                 */
    {
      vmfree(vmblock);              /* free buffers, report error     */
      return 2;
    }
    m[i] = (y[i + 1] - y[i]) / hi;
  }


  if (perio)                  /* periodic Akima interpolation ?       */
    m[-2]     = m[*n - 2],    /* find end point slopes by periodicity */
    m[-1]     = m[*n - 1],
    m[*n]     = m[0],
    m[*n + 1] = m[1];

  else                              /* standard Akima interpolation ? */
                                    /* compute end point slopes       */
    m[-2]     = THREE * m[0]      - TWO * m[1],
    m[-1]     = TWO   * m[0]      -       m[1],
    m[*n]     = TWO   * m[*n - 1] -       m[*n - 2],
    m[*n + 1] = THREE * m[*n - 1] - TWO * m[*n - 2];


  for (i = 0; i <= *n; i++)    /* compute left and right sided slopes */
  {                            /* tL[i] and tR[i], i=0,...,n; tR[i]   */
                               /* is stored in b[i] so that tR can be */
                               /* saved.                              */
    REAL him1, hi;
    if (i == 0)
    {
      if (perio) him1 = x[*n] - x[*n - 1];
      else       him1 = x[2]  - x[1];
      hi = x[i + 1] - x[i];
    }
    else if (i == *n)
    {
      him1 = x[i] - x[i - 1];
      if (perio) hi = x[1]      - x[0];
      else       hi = x[*n - 1] - x[*n - 2];
    }
    else
      him1 = x[i]     - x[i - 1],
      hi   = x[i + 1] - x[i];
    L = FABS(m[i - 1] - m[i - 2]) * him1;
    R = FABS(m[i + 1] - m[i])     * hi;
    nenner = L + R;
    if (nenner >= MACH_EPS)
    {
      alpha = L / nenner;
      tL[i] = m[i - 1] + alpha * (m[i] - m[i - 1]);
      if (i < *n)
        b[i] = tL[i];
    }
    else
    {
      tL[i] = m[i - 1];
      if (i < *n)
        b[i] = m[i];
    }
  }


  for (i = 0; i < *n; i++)            /* compute spline coefficients  */
    hi   = x[i + 1] - x[i],
    c[i] = (THREE * m[i] - TWO * b[i] - tL[i + 1]) / hi,
    d[i] = (b[i] + tL[i + 1] - TWO * m[i]) / sqr(hi);


  vmfree(vmblock);                        /* free buffers, report run */
  return 0;
}
/*.BA*/



/*.FE{C 13.2}{Renner Subsplines}{Renner Subsplines}*/

/*.BE*/
/* ------------------------------------------------------------------ */

static REAL norm22        /* (Euclidean norm of an vector)^2 .........*/
/*.IX{norm22}*/
        (
         REAL v[],                 /* input vector ...................*/
         int  n                    /* number of vector elements ......*/
        )                          /* square of norm .................*/

/***********************************************************************
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, sqr                                                      *
***********************************************************************/

{
  REAL norm;                              /* suare of 2 norm of v     */

  for (norm = ZERO; n-- != 0; )
    norm += sqr(*v++);

  return norm;
}



/* ------------------------------------------------------------------ */

static int ecke_runden    /* round corner P[i] for Renner subsplines  */
/*.IX{ecke\unt runden}*/
        (
         int      *n,              /* number of final node ...........*/
         int      nmax,            /* upper index limit for P ........*/
         int      dim,             /* space dimension (2 or 3) .......*/
         int      i,               /* index of corner to be rounded ..*/
         REAL     *P[],            /* nodes ..........................*/
         REAL     *s[],            /* secant vectors .................*/
         REAL     beta             /* rounding parameter .............*/
        )                          /* error code .....................*/

/***********************************************************************
* For renner() this function rounds corner P[i] by replacing it by two *
* additional nodes that are slightly moved.                            *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n      Number of final node. Numbering starts at zero. (n > 3)       *
* nmax   upper index limit on P (for checking whether rounding can     *
*        proceed in the allotted space)                                *
* i      index of the node which represents the corner to be rounded   *
* P      [0..n,0..dim-1] node vector.                                  *
*        Adjacent points may not coincide.                             *
* s      [-2..n+1,0..dim-1] secant vector of the curve                 *
*        The vectors P and s must be large enough to accommodate       *
*        maximally [(n+1)/2] extra points from rounding.               *
* beta   determines amount of rounding. Each corner point is replaced  *
*        by two points obtained by moving the corner point towards its *
*        left and right neighbor a short distance that depends on beta.*
*        If beta lies outside the interval (0,1) no rounding takes     *
*        place.                                                        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* n      new, possibly larger value for n                              *
* P,s    vectors as on input except for additional points              *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* 0: no error                                                          *
* 4: nmax not large enough for all corner roundings                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, copy_vector, SQRT, norm22, min                                 *
***********************************************************************/

{
  REAL L, R, B,         /* aux variable for weights                   */
       lambda,          /* weight for moving to back                  */
       my,              /* weight for moving to front                 */
       hilf[3];         /* difference of two nodes                    */
  int  j,               /* loop variable                              */
       k;               /* loop variable                              */


  if (*n == nmax)            /* no more room for additional rounding? */
    return 4;

  copy_vector(s[*n + 2], s[*n + 1], dim);  /* shift secant vectors    */
  for (j = *n; j >= i; j--)                /* s[i] ... s[n+1] and     */
  {                                        /* points P[i] ... P[n]    */
    copy_vector(P[j + 1], P[j], dim);      /* one index up            */
    copy_vector(s[j + 1], s[j], dim);
  }

  /* ------------ shift P[i] and P[i+1] a little bit ---------------- */

  for (k = dim; k-- != 0; )
    hilf[k] = P[i][k] - P[i - 1][k];
  L = SQRT(norm22(hilf, dim));

  for (k = dim; k-- != 0; )
    hilf[k] = P[i + 2][k] - P[i + 1][k];
  R = SQRT(norm22(hilf, dim));

  B      = beta * min(L, R);
  lambda = B / L;
  my     = B / R;

  for (k = dim; k-- != 0; )
    P[i    ][k] -= lambda * (P[i    ][k] - P[i - 1][k]),
    P[i + 1][k] += my     * (P[i + 2][k] - P[i + 1][k]);

  for (j = i - 1; j <= i + 1; j++)        /* compute three new        */
    for (k = dim; k-- != 0; )             /* secant vectors           */
      s[j][k]  = P[j + 1][k] - P[j][k];   /* s[i-1], s[i], s[i+1]     */

  (*n)++;



  return 0;
}



/* ------------------------------------------------------------------ */

static void normieren     /* norm vector s into s0 ...................*/
/*.IX{normieren}*/
        (
         REAL s[],                 /* input vector ...................*/
         REAL s0[],                /* corresponding unit vector ......*/
         int  n                    /* vector length ..................*/
        )

/***********************************************************************
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, sqr, SQRT, MACH_EPS                                      *
***********************************************************************/

{
  REAL norm;                            /* Euclidean norm of s        */
  int  i;                               /* loop variable              */

  for (norm = ZERO, i = n; i-- != 0; )
    s0[i]  = s[i],
    norm  += sqr(s[i]);
  norm = SQRT(norm);
  if (norm >= MACH_EPS)                 /* norm not vanishing?        */
    for (i = n; i-- != 0; )             /* divide vector elements     */
      s0[i] /= norm;                    /* by norm                    */
}



/* ------------------------------------------------------------------ */

static REAL flaeche       /* compute area between s and t ............*/
/*.IX{flaeche}*/
        (
         REAL s[],                 /* 1. normalized secant vector ....*/
         REAL t[],                 /* 2. normalized secant vector ....*/
         int  n                    /* vector length...................*/
        )                          /* area ...........................*/

/***********************************************************************
* This auxiliary function for renner() computes the area of the        *
* parallelogram spanned by the two normalized secant vectors s and t.  *
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL, FABS, SQRT, ONE, sqr                                           *
***********************************************************************/

{
  if (n == 2)                                   /* two-dimensional?   */
    return FABS(s[0] * t[1] - s[1] * t[0]);
  else                                          /* three-dimensional? */
    return SQRT(ONE - sqr(s[0] * t[0] +
                          s[1] * t[1] +
                          s[2] * t[2]));
}



/* ------------------------------------------------------------------ */
/*.BA*/

int renner                /* compute renner subspline coefficients ...*/
/*.IX{renner}*/
        (
         int  *n,                  /* number of last node ............*/
         int  nmax,                /* upper index limit for P ........*/
         int  dim,                 /* space dimension (2 or 3) .......*/
         REAL *P[],                /* nodes ..........................*/
         REAL beta,                /* rounding parameter .............*/
         REAL T[],                 /* lengths of parameter intervals  */
         REAL *b[],                /* coeffizients for t^1 ...........*/
         REAL *c[],                /* coeffizients for t^2 ...........*/
         REAL *d[]                 /* Koeffizients for t^3 ...........*/
        )                          /* error code .....................*/

/***********************************************************************
* Compute the coefficients of an interpolating  Renner subspline.      *
.BE*)
*                                                                      *
* The interpolating function can be evaluated using rennwert(). A table*
* of values can be obtained via renntab(). These functions, however,   *
* require monotonically increasing parameter values while the Renner   *
* parameters start over at 0 at each node. Hence one must add up the   *
* T[i] for a monotonic parametrization Tsum[i]. See the following      *
* piece of code for this :                                             *
*                                                                      *
*             for (Tsum[0] = 0, i = 1; i <= n; i++)                    *
*               Tsum[i] = Tsum[i-1] + T[i-1];                          *
*                                                                      *
* Then use Tsum instead of T in rennwert() and renntab() .             *
* Due to the special properties of Renner subsplines, Tsum[i] is in    *
* fact a good approximation for the arc length up to node i.           *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n     Number of last node. Numbering starts at zero.  (n > 3)        *
* nmax  upper index limit for the vectors P, T, b, c, d (to insure     *
*       that all corners can be rounded and to dimension aux vectors)  *
* P     [0..n,0..dim-1] node vector with the points to be interpolated.*
*       Any two consecutive points must be different.                  *
* beta  rounding parameter:  Each corner point is replaced             *
*       by two points obtained by moving the corner point towards its  *
*       left and right neighbor a short distance that depends on beta. *
*       If beta lies outside the interval (0,1) no rounding takes place*
*       The vectors  x, y, t, ax..dx and ay..dy must be large enough   *
*       to accomodate all potential new nodes from rounded corners.    *
*       There may be  [(n+1)/2] additional new nodes.                  *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* T      [0..n-1] vector with parameter interval lengths for the curve *
*        segments from P[i] to P[i+1], i=0,...,n-1                     *
* b  \   [0..n-1,0..dim-1] coefficient vectors of subspline function,  *
* c   >  which is given for i = 0,...,n-1 by                           *
* d  /                                                                 *
*              (a[i][0])     (b[i][0])       (c[i][0])       (d[i][0]) *
*        S(t)= (       ) + t*(       ) + t^2*(       ) + t^3*(       ) *
*              (a[i][1])     (b[i][1])       (c[i][1])       (d[i][1]) *
*                                                                      *
*        for dim=2 and 0 <= t <= T[i] resp.                            *
*                                                                      *
*              (a[i][0])     (b[i][0])       (c[i][0])       (d[i][0]) *
*        S(t)= (a[i][1]) + t*(b[i][1]) + t^2*(c[i][1]) + t^3*(d[i][1]) *
*              (a[i][2])     (b[i][2])       (c[i][2])       (d[i][2]) *
*                                                                      *
*        for dim=3 and 0 <= t <= T[i].                                 *
*        a corresponds to P and therefore has an additional entry      *
*        a[n].                                                         *
*                                                                      *
* And only for  beta in (0,1):                                         *
* n      new number of nodes after rounding corners                    *
* P      node coordinates, augmented number from rounding corners      *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* 0: no error                                                          *
* 1: n < 4                                                             *
* 2: two consecutive points coincide                                   *
* 3: lack of memory                                                    *
* 4: lack of room for additional nodes while rounding                  *
* 5: dim < 2  or  dim > 3                                              *
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, ONE, vminit, vmalloc, MATRIX, vmcomplete, vmfree,        *
* norm_max, MACH_EPS, copy_vector, THREE, TWO, normieren, flaeche,     *
* ecke_runden, norm22, SIX, skalprod, SQRT, sqr                        *
.BA*)
***********************************************************************/
/*.BE*/

{
  void     *vmblock;   /* List of dynamic allocations                 */
  int      m,          /* n after before corners are rounded          */
           fehler,     /* error code from ecke_runden()               */
           rundung,    /* flag: round corners?                        */
           anf = 0,    /* point index where check for corners begins  */
           end = 0,    /* point index where check for corners ends    */
           i,          /* loop variable                               */
           j;          /* loop variable                               */
  REAL     **tL,       /* [0..nmax] vector with left sided unit       */
                       /* tangent vectors                             */
                       /* (tR is not needed her because it can be     */
                       /* replaced by b.)                             */
           **s,        /* [-2..nmax+1] vector with secant vectors     */
           s0[4][3],   /* storage for the four normalized secant      */
                       /* vectors that are needed when computing      */
                       /* tL[i] and tR[i] und rounding corner P[i]    */
           *s0im2,     /* pointer to the unit secant vector s0[i-2]   */
           *s0im1,     /* pointer to the unit secant vector s0[i-1]   */
           *s0i,       /* pointer to the unit secant vector s0[i]     */
           *s0ip1,     /* pointer to the unit secant vector s0[i+1]   */
           Fim2,       /* area F(s0[i-2],s0[i-1])                     */
           Fim1,       /* area F(s0[i-1],s0[i])                       */
           Fi,         /* area F(s0[i],  s0[i+1])                     */
           hilf[3],    /* first difference P[0] - P[n], then not yet  */
                       /* normalized tangent vector tL[i], finally    */
                       /* sum of tR[i] and tL[i+1]                    */
           alpha,      /* weight factor for computing left and right  */
                       /* sided unit tangent vectors from adjacent    */
                       /* secant vectors                              */
           nenner,     /* denominator used for alpha                  */
           A, B, C,    /* aux variables for computing the length of   */
                       /* parameter intervals T[i]                    */
           TiInv,      /* reciprocal of T[i]                          */
           TiInvSq;    /* square of above                             */

/* cyclically shift the unit secant vectors and parallelogram areas   */
/* needed in the environment of P[i]                                  */
#define ZYKLISCH_TAUSCHEN_s0i_Fi \
{                \
  REAL *tmp;     \
  tmp   = s0im2; \
  s0im2 = s0im1; \
  s0im1 = s0i;   \
  s0i   = s0ip1; \
  s0ip1 = tmp;   \
  Fim2  = Fim1;  \
  Fim1  = Fi;    \
}


  if (*n < 4)                           /* too few nodes?             */
    return 1;

  if (dim < 2 || dim > 3)               /* no curve in R2 or R3?      */
    return 5;

  m       = *n;
  rundung = (beta > ZERO && beta < ONE);


  /* ---------- allocate buffers for aux vectors:            -------- */
  /* ---------- 1 [0..nmax] vector and 1 [-2..nmax+1] vector -------- */

  vmblock = vminit();                                   /* initialize */

  tL = (REAL **)vmalloc(vmblock, MATRIX, nmax + 1, dim);
  s  = (REAL **)vmalloc(vmblock, MATRIX, nmax + 4, dim);

  if (! vmcomplete(vmblock))                      /* lack of memory?  */
  {
    vmfree(vmblock);
    return 3;
  }

  s  += 2;              /* make sure both s can be used with negative */
                        /* indices as required by our algorithm       */


  for (i = 0; i < m; i++)               /* compute interior secants   */
  {                                     /* s[i], i=0,...,m-1          */
    for (j = dim; j-- != 0; )
      s[i][j] = P[i + 1][j] - P[i][j];
    if (norm_max(s[i], dim) < MACH_EPS)   /* if two consecutive nodes */
    {                                     /* coincide  =>  error      */
      vmfree(vmblock);
      return 2;
    }
  }


  /* ------- compute additional secant vectors                 ------ */
  /* ------- s[-2], s[-1], s[m], s[m+1]                        ------ */

  for (j = dim; j-- != 0; )
    hilf[j] = P[0][j] - P[m][j];

  if (norm_max(hilf, dim) < MACH_EPS)      /* closed curve?           */
  {
    copy_vector(s[-2],    s[m - 2], dim);  /* form secants at end     */
    copy_vector(s[-1],    s[m - 1], dim);  /* points cyclically       */
    copy_vector(s[m],     s[0],     dim);
    copy_vector(s[m + 1], s[1],     dim);
  }

  else                                     /* open curve?             */
                                           /* compute secant vectors  */
                                           /* at end points from      */
                                           /* adjacent data           */
    for (j = dim; j-- != 0; )
      s[-2][j]    = THREE * s[0][j]     - TWO * s[1][j],
      s[-1][j]    = TWO   * s[0][j]     -       s[1][j],
      s[m][j]     = TWO   * s[m - 1][j] -       s[m - 2][j],
      s[m + 1][j] = THREE * s[m - 1][j] - TWO * s[m - 2][j];


  s0im2 = s0[0];                           /* initialize the pointers */
  s0im1 = s0[1];                           /* to the four             */
  s0i   = s0[2];                           /* unit secant vectors     */
  s0ip1 = s0[3];

  normieren(s[-2], s0im2, dim);            /* compute unit secant     */
  normieren(s[-1], s0im1, dim);            /* vectors                 */
  normieren(s[0],  s0i,   dim);            /* s0[-2], s0[-1], s0[0]   */
  Fim2 = flaeche(s0im2, s0im1, dim);       /* and the corresponding   */
  Fim1 = flaeche(s0im1, s0i, dim);         /* areas                   */

  if (rundung)                                 /* round corners?      */

    if (norm_max(hilf, dim) < MACH_EPS)        /* closed curve?       */
    {
      normieren(s[1], s0ip1, dim);
      Fi = flaeche(s0i, s0ip1, dim);
      if (Fim2 + Fi <  MACH_EPS &&             /* Is P[0] a corner?   */
          Fim1      >= MACH_EPS)
      {
        for (i = -2; i < m + 1; i++)           /* shift points        */
        {                                      /* cyclically so that  */
          if (i >= 0 && i < m)                 /* P[n-1] becomes the  */
            copy_vector(P[i], P[i + 1], dim);  /* corner              */
          copy_vector(s[i], s[i + 1], dim);    /* => no special       */
        }                                      /* treatment for       */
        copy_vector(P[m], P[0], dim);          /* corner P[0] needed  */
        copy_vector(s[m + 1], s[1], dim);
        ZYKLISCH_TAUSCHEN_s0i_Fi;
      }
      anf = 1;
      end = m - 1;
    }
    else                                       /* open curve?         */
      anf = 2,
      end = m - 2;


  for (i = 0; i <= *n; i++)    /* compute left and right sided unit   */
  {                            /* tangent vectors and round corners   */
                               /* if desired and possible             */
    normieren(s[i + 1], s0ip1, dim);
    Fi     = flaeche(s0i, s0ip1, dim);
    nenner = Fim2 + Fi;
    if (nenner < MACH_EPS)                /* vanishing denominator?   */
    {
      if (rundung              &&         /* rounding desired?        */
          i >= anf && i <= end &&         /* allowed index range?     */
          Fim1 >= MACH_EPS                /* P[i-1], P[i], P[i+1]     */
         )                                /* not collinear?           */
      {
        fehler = ecke_runden(n, nmax, dim, i, P, /* round corner P[i] */
                             s, beta);
        if (fehler)                              /* failed?           */
        {
          vmfree(vmblock);
          return fehler;
        }
        normieren(s[i - 2], s0im2, dim);    /* recompute unit secant  */
        normieren(s[i - 1], s0im1, dim);    /* vectors and            */
        normieren(s[i], s0i, dim);          /* parallelogram areas    */
        Fim2 = flaeche(s0im2, s0im1, dim);  /* affected by rounding   */
        Fim1 = flaeche(s0im1, s0i, dim);


        end++;                        /* one more point to be checked */
        i--;                          /* P[i], P[i+1] new =>          */
      }                               /* inspect again                */

      else                                /* corner P[i] not rounded? */
      {
        copy_vector(tL[i], s0im1, dim);
        if (i < *n)
          copy_vector(b[i], s0i, dim);
        ZYKLISCH_TAUSCHEN_s0i_Fi;
      }
    }

    else                                  /* P[i] is no corner?       */
    {
      alpha  =  Fim2 / nenner;
      for (j = dim; j-- != 0; )
        hilf[j] =  s[i - 1][j] + alpha * (s[i][j] - s[i - 1][j]);
      normieren(hilf, tL[i], dim);
      if (i < *n)
        copy_vector(b[i], tL[i], dim);
      ZYKLISCH_TAUSCHEN_s0i_Fi;
    }
  }


  for (i = 0; i < *n; i++)  /* compute lengths of parameter intervals */
  {                         /* and spline coefficients                */
                            /* (tR is already stored in b and thus    */
                            /* could be omitted.)                     */

    for (j = dim; j-- != 0; )               /* hilf = tR[i] + tL[i+1] */
      hilf[j] = b[i][j] + tL[i + 1][j];
    A    = (REAL)16.0 - norm22(hilf, dim);
    B    = SIX        * skalprod(s[i], hilf, dim);
    C    = (REAL)36.0 * norm22(s[i], dim);
    T[i] = (-B + SQRT(sqr(B) + A * C)) / A;

    TiInv   = 1 / T[i];                     /* compute the remaining  */
    TiInvSq = sqr(TiInv);                   /* subspline coeffizients */
    for (j = dim; j-- != 0; )
      c[i][j] = (THREE * TiInv * s[i][j] -
                 TWO * b[i][j] - tL[i + 1][j]) * TiInv,
      d[i][j] = (hilf[j] - TWO * TiInv * s[i][j]) * TiInvSq;
  }


  vmfree(vmblock);
  return 0;
}



/* ------------------------------------------------------------------ */

static REAL polwert       /* compute  a + b*x + c*x^2 + d*x^3  .......*/
        (
         REAL x,                   /* place of evaluation ............*/
         REAL a,                   /* coeffizient of x^0 .............*/
         REAL b,                   /* coeffizient of x^1 .............*/
         REAL c,                   /* coeffizient of x^2 .............*/
         REAL d,                   /* coeffizient of x^3 .............*/
         REAL dp[]                 /* 1st,2nd,3rd derivative in x ....*/
        )                          /* polynomial value in x ..........*/

/***********************************************************************
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL, THREE, TWO                                                     *
***********************************************************************/

{
  REAL hilf1,             /* temporary storage for often used         */
       hilf2,             /* expressions in polynomial evaluation     */
       hilf3;

  hilf1 = THREE * d;                        /* evaluate polynom via   */
  hilf2 = TWO   * c;                        /* Horner scheme          */
  hilf3 = TWO   * hilf1;
  dp[0] = (hilf1 * x + hilf2) * x + b;
  dp[1] = hilf3 * x + hilf2;
  dp[2] = hilf3;

  return ((d * x + c) * x + b) * x + a;
}



/* ------------------------------------------------------------------ */
/*.BA*/

void rennwert             /* evaluation of a Renner subspline ........*/
/*.IX{rennwert}*/
        (
         int  n,                   /* number of spline pieces ........*/
         int  dim,                 /* space dimension (2 or 3) .......*/
         REAL twert,               /* point of evaluation ............*/
         REAL t[],                 /* nodes ..........................*/
         REAL *a[],                /* spline coeff. of (t-t[i])^0 ....*/
         REAL *b[],                /* spline coeff. of (t-t[i])^1 ....*/
         REAL *c[],                /* spline coeff. of (t-t[i])^2 ....*/
         REAL *d[],                /* spline coeff. of (t-t[i])^3 ....*/
         REAL s[],                 /* point on the spline curve ......*/
         REAL ds[4][3]             /* 0th - 3rd derivative of spline  */
        )

/***********************************************************************
* compute functional value and derivatives of a parametric cubic       *
* polynomial spline function in R2 or R3                               *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n        index of the last parameter value in t                      *
* dim      dimension of space where the spline curve lies:             *
*          planar curve for dim=2, spatial curve for dim=3             *
* twert    parameter value whose functional value derivaties must be   *
*          computed                                                    *
* t        [0..n] vector with the parameter values of the nodes P[i],  *
*          on which the spline coefficients were based                 *
* a,b,c,d  [0..n-1,0..dim-1] vectors with the spline coefficients      *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* s        functional value of the spline for t = twert                *
* ds       [0..3,0..2] matrix mit functional values and derivatives.   *
*          ds[0] = s, and ds[i] is the ith derivative of the spline    *
*          at t = twert, i = 1(1)3.                                    *
*          (Alle further derivatives are zero.)                        *
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL, intervall, polwert                                             *
.BA*)
***********************************************************************/
/*.BE*/

{
  static int i;       /* number of node interval for twert            */
  REAL   abl[3];      /* 1st, 2nd, 3rd derivative of spline at twert  */
  int    j,           /* loop variable                                */
         k;           /* loop variable                                */


  /* -- For repeated calls of rennwert() i must be determined in   -- */
  /* -- intervall() only when xwert has moved across a node.       -- */

  if (twert < t[i] || twert >= t[i + 1])
    i = intervall(n, twert, t);

  for (j = dim; j-- != 0; )
  {
    s[j]     = polwert(twert - t[i], a[i][j], b[i][j], c[i][j], d[i][j],
                       abl);
    ds[0][j] = s[j];
    for (k = 1; k < 4; k++)
      ds[k][j] = abl[k - 1];
  }
}



/* ------------------------------------------------------------------ */

static void hrenntab      /* auxiliary function for renntab() ........*/
/*.IX{hrenntab}*/
        (
         REAL ta,                  /* left border ....................*/
         REAL te,                  /* right border ...................*/
         REAL dt,                  /* step size ......................*/
         REAL ti,                  /* only node in [ta,te] ...........*/
         int  dim,                 /* space dimension (2 or 3) .......*/
         REAL a[],                 /* spline coeffic. of (t-t[i])^0 ..*/
         REAL b[],                 /* spline coeffic. of (t-t[i])^1 ..*/
         REAL c[],                 /* spline coeffic. of (t-t[i])^2 ..*/
         REAL d[],                 /* spline coeffic. of (t-t[i])^3 ..*/
         REAL *tab[],              /* table with curve points ........*/
         int  anzahl,              /* upper index limit of tab .......*/
         int  *lt                  /* highest tab index used .........*/
        )

/***********************************************************************
* This function tabulates a parametric cubic spline function inside    *
* one parameter interval from ta to te with step size dt. As the       *
* the computation is restricted to just one parameter interval, we     *
* must not distinguish between several polynomials.                    *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* ta       left end point of tabulation interval                       *
* te       right end point of tabulation interval                      *
* dt       step size                                                   *
* ti       left end node                                               *
* dim      dimension of space where the spline curve lies:             *
*          planar curve for dim=2, spatial curve for dim=3             *
* a,b,c,d  coefficients for spline in the node interval                *
* tab      [0..anzahl] vector with old table values                    *
* anzahl   upper index limit for tab                                   *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* tab      [0..anzahl] vector with old and new table values            *
* lt       index of last entry used in table                           *
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  REAL t0,
       t1;
  int  j;

  tab += *lt;
  for (t0 = ta; t0 < te; t0 += dt)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    tab++;
    t1 = t0 - ti;
    for (j = dim; j-- != 0; )
      (*tab)[j] = ((d[j] * t1 + c[j]) * t1 + b[j]) * t1 + a[j];
  }
}


/* ------------------------------------------------------------------ */
/*.BA*/

int renntab               /* table of values of a Renner subspline ...*/
/*.IX{renntab}*/
        (
         int  n,                   /* number of spline pieces ........*/
         REAL tanf,                /* left end point of interval .....*/
         REAL tend,                /* right end point of interval ....*/
         REAL delt,                /* step size ......................*/
         int  anzahl,              /* maximal size of table ..........*/
         REAL t[],                 /* nodes ..........................*/
         int  dim,                 /* space dimension (2 or 3) .......*/
         REAL *a[],                /* spline coeffic. of (t-t[i])^0 ..*/
         REAL *b[],                /* spline coeffic. of (t-t[i])^1 ..*/
         REAL *c[],                /* spline coeffic. of (t-t[i])^2 ..*/
         REAL *d[],                /* spline coeffic. of (t-t[i])^3 ..*/
         REAL *sptab[],            /* table of spline values .........*/
         int  *lentab              /* actual table length ............*/
        )                          /* error code .....................*/

/***********************************************************************
* Make a table of values of a cubic spline of a parametric cubic       *
* polynomial spline function in R2 or R3.                              *
.BE*)
* Due to rounding, we might tabulate in duplicate near nodes and at    *
* tend.                                                                *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n        Index of final node                                         *
* tanf \   define the parameter interval where the spline function is  *
* tend /   to be tabulated (tabulation interval)                       *
* delt     step size for table. Table values are computed for          *
*          t = tanf(delt)tend   in such a way that all t[i] from       *
*          [tanf,tend] show up as well.                                *
* t        [0..n] vector with the parameter values of the spline       *
* dim      dimension of space where the spline curve lies:             *
*          planar curve for dim=2, spatial curve for dim=3             *
* a,b,c,d  [0..n-1] vectors with spline coefficients.                  *
*          a[n] may also be used.                                      *
* anzahl   upper limit index for sptab                                 *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* sptab    [0..anzahl] vector with table values. Dimension `anzahl'    *
*          must be chosen large enough to offer storage for all        *
*          table values.                                               *
* lentab   last index used in table                                    *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* 0: no error                                                          *
* 1: tanf > tend                                                       *
* 2: delt <= 0                                                         *
*                                                                      *
* Global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, intervall, hrenntab, TWO                                 *
.BA*)
***********************************************************************/
/*.BE*/

{
  int anf,    /* Index of spline segment for tanf (0..n-1)            */
      end,    /* Index of spline segment for tend (0..n-1)            */
      anf2,   /* index of nearest node <= tanf (-1..n)                */
      end2,   /* index of nearest node <= tend (-1..n)                */
      i,      /* Laufvariable                                         */
      j;      /* Laufvariable                                         */

  if (tanf > tend)                      /* invalid input parameters?  */
    return 1;
  if (delt <= ZERO)
    return 2;

  anf2 = anf = intervall(n, tanf, t);   /* find spline segments for   */
  end2 = end = intervall(n, tend, t);   /* `tanf' and `tend'          */
  if (tanf < t[0])                      /* `tanf' to the left of the  */
    anf2--;                             /* interpolation interval?    */
  if (tend > t[n])                      /* `tend' to the right of the */
    end2++;                             /* interpolation interval?    */

  *lentab = -1;

  if (anf2 < end2)                      /* at least one node in the   */
  {                                     /* tabulation interval?       */
    hrenntab(tanf, t[anf2 + 1], delt,   /* tabulate from `tanf' to    */
             t[anf], dim, a[anf],       /* the first node lying in    */
             b[anf], c[anf], d[anf],    /* the tabulation interval    */
             sptab, anzahl, lentab);

    for (i = anf2 + 1; i < end2; i++)   /* tabulate in node intervals */
       hrenntab(t[i], t[i + 1], delt,   /* lying completely in the    */
                t[i], dim, a[i], b[i],  /* tabulation interval        */
                c[i], d[i], sptab,
                anzahl, lentab);

    tanf = t[end2];                     /* continue the table at      */
                                        /* `tanf'                     */
    if (end2 == n)                      /* make the spline value at   */
      if (*lentab < anzahl)             /* t[n] be entered exactly if */
      {                                 /* it lies in the             */
        (*lentab)++;                    /* tabulation interval        */
        for (j = dim; j-- != 0; )
          sptab[*lentab][j] = a[n][j];
        tanf += delt;
      }

    hrenntab(tanf, tend, delt, t[end],  /* tabulate from the last     */
             dim, a[end], b[end],       /* node lying in the          */
             c[end], d[end], sptab,     /* tabulation interval        */
             anzahl, lentab);           /* to `tend'                  */
  }

  else                                  /* no node at all lying in    */
                                        /* the tabulation interval?   */
    hrenntab(tanf, tend, delt, t[anf],
             dim, a[anf], b[anf],
             c[anf], d[anf], sptab,
             anzahl, lentab);

  hrenntab(tend, tend + delt / TWO,     /* special treatment for the  */
           delt, t[end], dim, a[end],   /* right border `tend' as it  */
           b[end], c[end], d[end],      /* normally is not considered */
           sptab, anzahl, lentab);      /* in the loop above          */

  return 0;
}

/* -------------------------- END subsplin.c ------------------------ */
