#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE spliwert.c ------------------------ */

/***********************************************************************
*                                                                      *
* Evaluate  spline functions                                           *
* ---------------------------                                          *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Source:               equivalent TP unit and QuickBASIC module       *
* Author:               Elmar Pohl (QuickBASIC)                        *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Date:                 8.13.1991                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>        /*  for  PI, sqr, intervall, FABS, POW,     */
                          /*       COS, SIN, MACH_EPS, REAL, TWO,     */
                          /*       THREE, ZERO, FIVE, FOUR, SIX       */
#include <spliwert.h>     /*  for  spwert, pspwert, hmtwert, pmtwert, */
                          /*       strwert                            */



/* ------------------------------------------------------------------ */
/*.BA*/

REAL spwert       /* evaluate a cubic spline .........................*/
/*.IX{spwert}*/
           (
            int  n,        /* number of spline pieces ................*/
            REAL xwert,    /* x-value ................................*/
            REAL a[],      /* Spline coefficients of (x-x[i])^0 ......*/
            REAL b[],      /* Spline coefficients of (x-x[i])^1 ......*/
            REAL c[],      /* Spline coefficients of (x-x[i])^2 ......*/
            REAL d[],      /* Spline coefficients of (x-x[i])^3 ......*/
            REAL x[],      /* nodes (x-values) .......................*/
            REAL ausg[]    /* 1st, 2nd and 3rd derivatives of spline .*/
           )               /* Functional value for spline ............*/

/***********************************************************************
* Compute functional and derivatives of a cubic spline.                *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:       Index of final node                                         *
* xwert:   place of evaluation                                         *
* a,b,c,d: [0..n-1] vectors with spline coefficients                   *
* x:       [0..n] vector of nodes (x-values)                           *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* ausg: [0..2] vector :                                                *
*       ausg[0] = first derivative at xwert,                           *
*       ausg[1] = second derivative,                                   *
*       ausg[2] = third derivative at xwert.                           *
*       all other derivatives are zero.                                *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* functional value at  xwert                                           *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, intervall, TWO, THREE                                          *
.BA*)
***********************************************************************/
/*.BE*/

{
  static int i = 0;   /* Number of node interval for  xwert           */
  REAL       hilf1,   /* aux variables                                */
             hilf2,
             hilf3;

  /* -- For repeated calls spwert() i must be determined in        -- */
  /* -- intervall() only when xwert moves across a node            -- */

  if (xwert < x[i] || xwert >= x[i + 1])
    i = intervall(n, xwert, x);

  /* ------- evaluate spline at xwert via Horner schems ------------- */

  xwert -= x[i];
  hilf1 = THREE * d[i];
  hilf2 = TWO   * c[i];
  hilf3 = TWO   * hilf1;
  ausg[0] = (hilf1 * xwert + hilf2) * xwert + b[i];
  ausg[1] = hilf3 * xwert + hilf2;
  ausg[2] = hilf3;

  return ((d[i] * xwert + c[i]) * xwert + b[i]) * xwert + a[i];
}



/* ------------------------------------------------------------------ */
/*.BA*/

void pspwert  /* Evaluate a parametric cubic spline ..................*/
/*.IX{pspwert}*/
            (
             int      n,        /* number of nodes ...................*/
             REAL     twert,    /* place for evaluation ..............*/
             REAL     t[],      /* nodes .............................*/
             REAL     ax[],     /* x spline coeff. of (t-t[i])^0 .....*/
             REAL     bx[],     /* x spline coeff. of (t-t[i])^1 .....*/
             REAL     cx[],     /* x spline coeff. of (t-t[i])^2 .....*/
             REAL     dx[],     /* x spline coeff. of (t-t[i])^3 .....*/
             REAL     ay[],     /* y spline coeff. of (t-t[i])^0 .....*/
             REAL     by[],     /* y spline coeff. of (t-t[i])^1 .....*/
             REAL     cy[],     /* y spline coeff. of (t-t[i])^2 .....*/
             REAL     dy[],     /* y spline coeff. of (t-t[i])^3 .....*/
             REAL     *sx,      /* x-coordinate, .....................*/
             REAL     *sy,      /* y-coordinate of spline value ......*/
             abl_mat1 ausp      /* 0 to third derivatives of spline ..*/
            )

/***********************************************************************
* Evaluate function and derivatives of a parametric cubic spline.      *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:            Index of last parameter node in t                      *
* twert:        Parameter value for evaluation                         *
* t:            [0..n] vector of parameter nodes for                   *
*               (X[i],Y[i]), i = 0,...,n                               *
* ax,bx,cx,dx:\                                                        *
* ay,by,cy,dy:/ [0..n-1] vectors of spline coefficients                *
*                                                                      *
* Output parameters :                                                  *
* ===================                                                  *
* sx:   value of spline component SX at t = twert                      *
* sy:   value of spline component SY at t = twert                      *
* ausp: [0..3,0..1] array with function and derivative values.         *
*       We have ausp[0,0] = sx,   ausp[0,1] = sy                       *
*       and ausp[i,0] denotes the  ith derivative of SX at t = twert,  *
*       i = 1,2,3.                                                     *
*       (All other derivatives are zero.)                              *
*       ausp[i,1] contains the derivatives of SY.                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* spwert, REAL                                                         *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL ausg[4];
  int  i;

  *sx = spwert(n, twert, ax, bx, cx, dx, t, ausg);
  ausp[0][0] = *sx;
  for (i = 1; i < 4; i++)
    ausp[i][0] = ausg[i - 1];
  *sy = spwert(n, twert, ay, by, cy, dy, t, ausg);
  ausp[0][1] = *sy;
  for (i = 1; i < 4; i++)
    ausp[i][1] = ausg[i - 1];
}



/* ------------------------------------------------------------------ */
/*.BA*/

REAL hmtwert        /* Evaluate a Hermite spline .....................*/
/*.IX{hmtwert}*/
            (
             int  n,        /* number of nodes - 1 ...................*/
             REAL x0,       /* place of evaluation ...................*/
             REAL a[],      /* Spline coefficient of (x-x[i])^0 ......*/
             REAL b[],      /* Spline coefficient of (x-x[i])^1 ......*/
             REAL c[],      /* Spline coefficient of (x-x[i])^2 ......*/
             REAL d[],      /* Spline coefficient of (x-x[i])^3 ......*/
             REAL e[],      /* Spline coefficient of (x-x[i])^4 ......*/
             REAL f[],      /* Spline coefficient of (x-x[i])^5 ......*/
             REAL x[],      /* n+1 nodes .............................*/
             REAL ausg[]    /* 1st to 5th derivatives of spline ......*/
            )               /* Function value of spline ..............*/

/***********************************************************************
* Evaluate  Hermite pline and its derivatives.                         *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:           Index of last node                                      *
* x0:          place of evaluation                                     *
* a,b,c,d,e,f: [0..n-1] vectors of spline coefficients                 *
* x:           [0..n] vector of nodes                                  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* ausg: [0..4] vector with values of the derivatives at x0 :           *
*       ausg[0] = first derivative at x0,                              *
*       ausg[1] = second derivative at x0,                             *
*       ausg[2] etc.                                                   *
*       ausg[3] ditto                                                  *
*       ausg[4] = fifth derivative                                     *
*       all other derivatives vanish.                                  *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* function value at  x0                                                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, intervall, TWO, THREE, FIVE, FOUR, SIX                         *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i;                 /* Number of interval with x0               */
  REAL B, C, D, E, F;     /* Spline coefficients for the interval i   */


  i = intervall(n, x0, x);
  x0 -= x[i];
  B = b[i];
  C = c[i];
  D = d[i];
  E = e[i];
  F = f[i];
  ausg[0] = (((FIVE * F * x0 + FOUR * E) * x0 + THREE * D) * x0 +
             TWO * C) * x0 + B;
  ausg[1] = (((REAL)20.0 * F * x0 + (REAL)12.0 * E) * x0 + SIX * D) *
            x0 + TWO * C;
  ausg[2] = ((REAL)60.0 * F * x0 + (REAL)24.0 * E) * x0 + SIX * D;
  ausg[3] = (REAL)120.0 * F * x0 + (REAL)24.0 * E;
  ausg[4] = (REAL)120.0 * F;


  return ((((F * x0 + E) * x0 + D) * x0 + C) * x0 + B) * x0 + a[i];
}



/* ------------------------------------------------------------------ */
/*.BA*/

void pmtwert  /* Evaluate a  parametric Hermite spline ...............*/
/*.IX{pmtwert}*/
            (
             int      n,        /* number of nodes - 1 ...............*/
             REAL     twert,    /* place of evaluation ...............*/
             REAL     t[],      /* nodes (x-values) ..................*/
             REAL     ax[],     /* x spline coeff. of (t-t[i])^0 .....*/
             REAL     bx[],     /* x spline coeff. of (t-t[i])^1 .....*/
             REAL     cx[],     /* x spline coeff. of (t-t[i])^2 .....*/
             REAL     dx[],     /* x spline coeff. of (t-t[i])^3 .....*/
             REAL     ex[],     /* x spline coeff. of (t-t[i])^4 .....*/
             REAL     fx[],     /* x spline coeff. of (t-t[i])^5 .....*/
             REAL     ay[],     /* y spline coeff. of (t-t[i])^0 .....*/
             REAL     by[],     /* y spline coeff. of (t-t[i])^1 .....*/
             REAL     cy[],     /* y spline coeff. of (t-t[i])^2 .....*/
             REAL     dy[],     /* y spline coeff. of (t-t[i])^3 .....*/
             REAL     ey[],     /* y spline coeff. of (t-t[i])^4 .....*/
             REAL     fy[],     /* y spline coeff. of (t-t[i])^5 .....*/
             REAL     *sx,      /* x-coordinate, .....................*/
             REAL     *sy,      /* y-coordinate of spline value ......*/
             abl_mat2 ausp      /* 0th to fifth derivatives of spline */
            )

/***********************************************************************
* Evaluate a parametric Hermite spline of degree 5 and its derivatives.*
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:                  Index of the last parameter value in t           *
* twert:              Parameter value for which to evaluate            *
* t:                  [0..n] vector of t parameter values for the nodes*
*                     (X[i],Y[i]), i = 0, ..., n                       *
* ax,bx,cx,dx,ex,fx:\ [0..n-1] vectors with spline coefficients        *
* ay,by,cy,dy,ey,fy:/                                                  *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* sx:   Function value of SX at t = twert                              *
* sy:   Function value of the spline component  SY at t = twert        *
* ausp: [0..5,0..1] array with function and derivative values:         *
*       ausp[0,0] = sx,   ausp[0,1] = sy,                              *
*       ausp[i,0] = ith derivative of spline component SX at twert     *
*       ausp[i,1] = ditto for SY    (i = 1, ..., 5)                    *
*       All other derivatives vanish.                                  *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hmtwert, REAL                                                        *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL ausg[6];
  int  i;

  *sx = hmtwert(n, twert, ax, bx, cx, dx, ex, fx, t, ausg);
  ausp[0][0] = *sx;
  for (i = 1; i < 6; i++)
    ausp[i][0] = ausg[i - 1];
  *sy = hmtwert(n, twert, ay, by, cy, dy, ey, fy, t, ausg);
  ausp[0][1] = *sy;
  for (i = 1; i < 6; i++)
    ausp[i][1] = ausg[i - 1];
}



/* ------------------------------------------------------------------ */
/*.BA*/

int strwert  /* Evaluate a transformed parametric cubic spline .......*/
/*.IX{strwert}*/
           (
            REAL phi,        /* place of evaluation  .................*/
            int  n,          /* number of nodes - 1 ..................*/
            REAL phin[],     /* angular nodes ........................*/
            REAL a[],        /* Spline coeff. of (phi-phin[i])^0 .....*/
            REAL b[],        /* Spline coeff. of (phi-phin[i])^1 .....*/
            REAL c[],        /* Spline coeff. of (phi-phin[i])^2 .....*/
            REAL d[],        /* Spline coeff. of (phi-phin[i])^3 .....*/
            REAL phid,       /* angle of plane rotation ..............*/
            REAL px,         /* coordinates of translation vector P ..*/
            REAL py,
            REAL ablei[],    /* 0th to third derivatives wrt. x ......*/
            REAL *xk,        /* x-coordinate,                         */
            REAL *yk,        /* y-coordinate of spline at phi ........*/
            REAL *c1,        /* 1st derivative (dr/dphi) at phi ......*/
            REAL *ckr        /* curvature of spline at phi ...........*/
           )                 /* error code ...........................*/

/***********************************************************************
* Evaluate a transformed parametric cubic spline of the form           *
*   s(phi) = a[i] + b[i](phi-phin[i]) + c[i](phi-phin[i])^2 +          *
*                                     + d[i](phi-phin[i])^3            *
* for phi in [phin[i],phin[i+1]], i=0, ..., n-1, and its first three   *
* derivatives.                                                         *
* Besides we also compute the cartesian coordinates (xk, yk) from the  *
* polar coordinates (phi, s(phi)), as well as the 1st derivative and   *
* the curvature of the curve at phi.                                   *
.BE*)
* NOTE: This evaluation function is not well suited form a table of    *
*       values  s(phi) or one of points (xk, yk) on the curve.         *
*       For this we recommend the function spwert().                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* phi:     place of evaluation (angle in radians)                      *
* n:       number of nodes in  phin                                    *
* phin:    [0..n] vector of angular nodes in radian measure            *
* a,b,c,d: [0..n-1] vectors with spline coefficients                   *
* phid:\   data of the transformation                                  *
* px:   >                                                              *
* py:  /                                                               *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ablei: [0..3] vector with function and derivative values:            *
*        ablei[0] = s(phi),                                       (S)  *
*        ablei[1] = first derivative s'(phi),                     (S1) *
*        ablei[2] = second derivative s''(phi),                   (S2) *
*        ablei[3] = third derivative s'''(phi),                        *
*        all other derivatives vanish.                                 *
* xk:\   cartesian coordinates for the curve at phi                    *
* yk:/                                                                 *
* c1:    1st derivative of the curve at  phi :                         *
*                 c1 = (S1 * sin(rho) + S * cos(rho)) /                *
*                      (S1 * cos(rho) - S * sin(rho))                  *
*        for rho = phi + phid                                          *
* ckr:   curvature of curve at  phi :                                  *
*        ckr = (2 * S1^2 - S * S2 + S^2) / ((S1^2 + S^2) ^ 1.5)        *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* fehler1 + 3 * fehler2.                                               *
* fehler1 concerns  c1 while fehler2 is from  ckr.                     *
* fehler1 can have the following values :                              *
* 0: no error                                                          *
* 1: numerator for c1 is zero.                                         *
* 2: the magnitude of the denominator expression for c1 differs from   *
*    zero but does not exceed 4 * MACHEPS; further computations contain*
*    large errors.                                                     *
* fehler2 returns the same values, but this time for  ckr.             *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* spwert, REAL, PI, MACH_EPS, sqr, FABS, POW, COS, SIN, TWO, ZERO,     *
* FOUR                                                                 *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  fehler1,
       fehler2,
       l;
  REAL fmasch,
       phix,
       rho,
       cosa,
       sina,
       hz,
       hn;


  fehler1 = fehler2 = 0;
  fmasch  = FOUR * MACH_EPS;

  if (phi < ZERO)                 /* adjust phi in phix so that lies  */
    l = (int)FABS(phi / TWO / PI) /* in the interval [0, 2 * pi]      */
        + 1,
    phix = l * TWO * PI - phi;
  else if (phi > TWO * PI)
    l = (int)(phi / TWO * PI),
    phix = phi - l * TWO * PI;
  else
    phix = phi;

  /* --- compute function and derivative values at phix ------------- */
  ablei[0] = spwert(n, phix, a, b, c, d, phin, ablei + 1);

  rho  = phix + phid;            /* compute (xk, yk), the derivative  */
  cosa = COS(rho);               /* and curvature                     */
  sina = SIN(rho);
  *xk  = ablei[0] * cosa + px;
  *yk  = ablei[0] * sina + py;
  hz   = ablei[1] * sina + ablei[0] * cosa;
  hn   = ablei[1] * cosa - ablei[0] * sina;
  if (hn == ZERO)
    fehler1 = 1;
  else
  {
    if (FABS(hn) <= fmasch)
      fehler1 = 2;
    *c1 = hz / hn;
  }
  hz = TWO * sqr(ablei[1]) - ablei[0] * ablei[2] + sqr(ablei[0]);
  hn = POW((sqr(ablei[1]) + sqr(ablei[0])), (REAL)1.5);
  if (hn == ZERO)
    fehler2 = 1;
  else
  {
    if (FABS(hn) <= fmasch)
      fehler2 = 2;
    *ckr = hz / hn;
  }


  return fehler1 + 3 * fehler2;
}

/* -------------------------- END spliwert.c ------------------------ */
