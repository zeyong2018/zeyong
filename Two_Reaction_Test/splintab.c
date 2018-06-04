#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE splintab.c ------------------------ */

/***********************************************************************
*                                                                      *
* Make table of values for splines                                     *
* --------------------------------                                     *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Source:               equivalent TP unit and QuickBASIC module       *
* Author:               Elmar Pohl (QuickBASIC)                        *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Date:                 7.8.1991                                       *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  intervall, COS, SIN, REAL, ZERO, TWO */
#include <splintab.h>   /*  for  sptab, partab, hmtab, pmtab, strtab  */

/* ------------------------------------------------------------------ */

static void sptabh      /* aux function for sptab() ..................*/
/*.IX{sptabh}*/
                  (
                   REAL xa,
                   REAL xe,
                   REAL dx,
                   REAL xi,
                   REAL a,
                   REAL b,
                   REAL c,
                   REAL d,
                   REAL xt[],
                   REAL yt[],
                   int  anzahl,
                   int  *lt
                  )

/***********************************************************************
* This function tabulates a cubic spline inside one subinterval from   *
* xa to xe with step size  dx.                                         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* xa:      left end point of tabulation interval                       *
* xe:      right end point of tabulation interval                      *
* dx:      step size                                                   *
* xi:      left end node                                               *
* a,b,c,d: coefficients for spline in the support interval             *
* xt:      [0..anzahl] vector with old x-values                        *
* yt:      [0..anzahl] vector with old y-values                        *
* anzahl:  upper index limit for xtab and ytab                         *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xt: [0..anzahl] vector with old and new x-values of the table        *
* yt: [0..anzahl] vector, ditto for y-values                           *
* lt: Index of last entry in table                                     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL
***********************************************************************/

{
  REAL x0,
       x1;

  xt += *lt, yt += *lt;
  for (x0 = xa; x0 < xe; x0 += dx)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    x1 = x0 - xi;
    *++xt = x0;
    *++yt = ((d * x1 + c) * x1 + b) * x1 + a;
  }
}


/* ------------------------------------------------------------------ */
/*.BA*/

int sptab       /* Table of values of a cubic spline .................*/
/*.IX{sptab}*/
         (
          int  n,         /* number of spline pieces ( = # nodes - 1) */
          REAL xanf,      /* left end point of interval ..............*/
          REAL xend,      /* right end point of interval .............*/
          REAL deltx,     /* step size ...............................*/
          int  anzahl,    /* maximal size of table ...................*/
          REAL x[],       /* nodes ...................................*/
          REAL a[],       /* Spline coefficients for (x-x[i])^0 ......*/
          REAL b[],       /* Spline coefficients for (x-x[i])^1 ......*/
          REAL c[],       /* Spline coefficients for (x-x[i])^2 ......*/
          REAL d[],       /* Spline coefficients for (x-x[i])^3 ......*/
          REAL xtab[],    /* x-coordinates of the table ..............*/
          REAL ytab[],    /* y-coordinates of the table ..............*/
          int  *lentab    /* actual table length .....................*/
         )                /* error code ..............................*/

/***********************************************************************
* Make a table of values of a cubic spline.                            *
* Due to rounding, we might tabulate in duplicate near nodes and at    *
* xend.                                                                *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:       Index of final node                                         *
* xanf:\   end points of the x interval for the table                  *
* xend:/                                                               *
* deltx:   step size for table. The table is evaluated for x = xanf,   *
*          xanf + deltx, ..., xend .                                   *
* x:       [0..n] vector of nodes                                      *
* a,b,c,d: [0..n-1] vectors with spline coefficients (a[n] may also    *
*          be used)                                                    *
* anzahl:  upper index limit for xtab and ytab                         *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xtab:   [0..anzahl] vector of x-values in table (Choose anzahl large *
*         enough !)                                                    *
* ytab:   [0..anzahl] vector of y-values in table with ytab[i] de-     *
*         noting the spline value at xtab[i], i=0, ..., lentab         *
* lentab: last index used in table                                     *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: xanf > xend                                                       *
* 2: deltx <= 0                                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* sptabh, REAL, intervall, ZERO, TWO                                   *
.BA*)
***********************************************************************/
/*.BE*/

{
  int anf,    /* Index of spline segment for xanf (0..n-1)            */
      end,    /* Index of spline segment for xend (0..n-1)            */
      anf2,   /* index of nearest node <= xanf (-1..n)                */
      end2,   /* index of nearest node <= xend (-1..n)                */
      i;      /* Loop variable                                        */

  if (xanf > xend)
    return 1;
  if (deltx <= ZERO)
    return 2;
  anf2 = anf = intervall(n, xanf, x);
  end2 = end = intervall(n, xend, x);
  if (xanf < x[0]) /* xanf to the left of the interpolation interval? */
    anf2--;
  if (xend > x[n]) /* xend to the right of the interpolation interval?*/
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    sptabh(xanf, x[anf2 + 1], deltx, x[anf], a[anf], b[anf],
           c[anf], d[anf], xtab, ytab, anzahl, lentab);
    for (i = anf2 + 1; i < end2; i++)
      sptabh(x[i], x[i + 1], deltx, x[i], a[i], b[i], c[i], d[i],
             xtab, ytab, anzahl, lentab);
    xanf = x[end2];                        /* continue table at xanf  */
    if (end2 == n)                         /* record data at x[n]     */
      if (*lentab < anzahl)                /* if it belongs to the    */
        xtab[++(*lentab)] = x[n],          /* interval of the table   */
        ytab[*lentab]     = a[n],
        xanf += deltx;
    sptabh(xanf, xend, deltx, x[end], a[end], b[end],
           c[end], d[end], xtab, ytab, anzahl, lentab);
  }
  else
    sptabh(xanf, xend, deltx, x[anf], a[anf], b[anf], c[anf],
           d[anf], xtab, ytab, anzahl, lentab);

  /* ----- special care with right end point xend ------------------- */

  sptabh(xend, xend + deltx / TWO, deltx, x[end], a[end], b[end],
         c[end], d[end], xtab, ytab, anzahl, lentab);

  return 0;
}

/* ------------------------------------------------------------------ */

static void partabh    /* aux function for  partab() ................*/
/*.IX{partabh}*/
                   (
                    REAL ta,
                    REAL te,
                    REAL dt,
                    REAL ti,
                    REAL ax,
                    REAL bx,
                    REAL cx,
                    REAL dx,
                    REAL ay,
                    REAL by,
                    REAL cy,
                    REAL dy,
                    REAL xt[],
                    REAL yt[],
                    int  anzahl,
                    int  *lt
                   )

/***********************************************************************
* This function tabulates a parametric cubic spline inside the interval*
* from ta to te with the step size  dt.                                *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* ta:           left end point for table                               *
* te:           right endpoint for table                               *
* dt:           step size                                              *
* ti:           left end point of support interval                     *
* ax,bx,cx,dx:\ coefficients for spline inside support interval        *
* ay,by,cy,dy:/                                                        *
* xt:           [0..anzahl] vector with x-values for table             *
* yt:           [0..anzahl] vector, ditto for y                        *
* anzahl:       upper index limit for xtab and ytab                    *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xt: [0..anzahl] vector with old and new x-values of the table        *
* yt: [0..anzahl] vector, ditto for y                                  *
* lt: Index of final entry in table                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL
***********************************************************************/

{
  REAL t0,
       t1;

  xt += *lt, yt += *lt;
  for (t0 = ta; t0 < te; t0 += dt)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    t1 = t0 - ti;
    *++xt = ((dx * t1 + cx) * t1 + bx) * t1 + ax;
    *++yt = ((dy * t1 + cy) * t1 + by) * t1 + ay;
  }
}


/* ------------------------------------------------------------------ */
/*.BA*/

int partab    /* Table of values for a parametric cubic spline .......*/
/*.IX{partab}*/
          (
           int  n,         /* number of spline pieces ................*/
           REAL tanf,      /* left end point of interval .............*/
           REAL tend,      /* right end point of interval ............*/
           REAL delt,      /* step size ..............................*/
           int  anzahl,    /* maximal size of table ..................*/
           REAL t[],       /* parameter nodes ........................*/
           REAL ax[],      /* x spline coefficients for (t-t[i])^0 ...*/
           REAL bx[],      /* x spline coefficients for (t-t[i])^1 ...*/
           REAL cx[],      /* x spline coefficients for (t-t[i])^2 ...*/
           REAL dx[],      /* x spline coefficients for (t-t[i])^3 ...*/
           REAL ay[],      /* y spline coefficients for (t-t[i])^0 ...*/
           REAL by[],      /* y spline coefficients for (t-t[i])^1 ...*/
           REAL cy[],      /* y spline coefficients for (t-t[i])^2 ...*/
           REAL dy[],      /* y spline coefficients for (t-t[i])^3 ...*/
           REAL xtab[],    /* x-coordinates of table .................*/
           REAL ytab[],    /* y-coordinates of table .................*/
           int  *lentab    /* actual size of table ...................*/
          )                /* error code .............................*/

/***********************************************************************
* Make a table of values for a parametric cubic spline.                *
* Due to rounding we might tabulate in duplicate near nodes or xend.   *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:            final parameter index for t                            *
* tanf:\        Parameter end points for table                         *
* tend:/        Funktion tabelliert werden soll.                       *
* delt:         step size of table.                                    *
* t:            [0..n] vector of parameter nodes                       *
* ax,bx,cx,dx:\ [0..n-1] vectors with spline coefficients              *
* ax,by,cy,dy:/ ax[n] and ay[n] may also be used                       *
* anzahl:       upper index limit for xtab and ytab                    *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xtab:   [0..anzahl] vector of x-values of the table                  *
* ytab:   [0..anzahl] vector, ditto for y                              *
* lentab: final index in table                                         *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: tanf > tend                                                       *
* 2: delt <= 0                                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* partabh, REAL, intervall, ZERO, TWO                                  *
.BA*)
***********************************************************************/
/*.BE*/

{
  int anf,    /* Index for spline for tanf (0..n-1)                   */
      end,    /* Index foe spline at tend  (0..n-1)                   */
      anf2,   /* Number of nearest node <= tanf (-1..n)               */
      end2,   /* Number of nearest node <= tend (-1..n)               */
      i;      /* Loop variable                                        */

  if (tanf > tend)
    return 1;
  if (delt <= ZERO)
    return 2;
  anf2 = anf = intervall(n, tanf, t);
  end2 = end = intervall(n, tend, t);
  if (tanf < t[0])   /* tanf to the left of interpolation interval ?  */
    anf2--;
  if (tend > t[n])   /* tend to the right of interpolation interval ? */
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    partabh(tanf, t[anf2 + 1], delt, t[anf], ax[anf], bx[anf],
            cx[anf], dx[anf], ay[anf], by[anf], cy[anf],
            dy[anf], xtab, ytab, anzahl, lentab);
    for (i = anf2 + 1; i < end2; i++)
       partabh(t[i], t[i + 1], delt, t[i], ax[i], bx[i], cx[i],
               dx[i], ay[i], by[i], cy[i], dy[i], xtab, ytab,
               anzahl, lentab);
    tanf = t[end2];                        /* continue table at  tanf */
    if (end2 == n)                         /* record function value at*/
      if (*lentab < anzahl)                /* t[n], if it belongs to  */
        xtab[++(*lentab)] = ax[n],         /* interval                */
        ytab[*lentab]     = ay[n],
        tanf += delt;
    partabh(tanf, tend, delt, t[end], ax[end], bx[end],
            cx[end], dx[end], ay[end], by[end], cy[end],
            dy[end], xtab, ytab, anzahl, lentab);
  }
  else
    partabh(tanf, tend, delt, t[anf], ax[anf], bx[anf],
            cx[anf], dx[anf], ay[anf], by[anf], cy[anf],
            dy[anf], xtab, ytab, anzahl, lentab);

  /* ----- special care at right end point  xend -------------------- */

  partabh(tend, tend + delt / TWO, delt, t[end], ax[end], bx[end],
          cx[end], dx[end], ay[end], by[end], cy[end], dy[end],
          xtab, ytab, anzahl, lentab);

  return 0;
}

/* ------------------------------------------------------------------ */

static void hmtabh      /* aux function for  hmtab() .................*/
/*.IX{hmtabh}*/
                  (
                   REAL xa,
                   REAL xe,
                   REAL dx,
                   REAL xi,
                   REAL a,
                   REAL b,
                   REAL c,
                   REAL d,
                   REAL e,
                   REAL f,
                   REAL xt[],
                   REAL yt[],
                   int  anzahl,
                   int  *lt
                  )

/***********************************************************************
* Tabulate a Hermite spline inside the interval from xa to xe with step*
* size  dx.                                                            *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* xa:          left end point of interval for table                    *
* xe:          right end point of same                                 *
* dx:          step size                                               *
* xi:          left end point of support interval                      *
* a,b,c,d,e,f: coefficienets of spline in support interval             *
* xt:          [0..anzahl] vector with old x-values                    *
* yt:          [0..anzahl] vector, ditto for y                         *
* anzahl:      upper index limit for xtab and ytab                     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xt: [0..anzahl] vector with old and new x-values in table            *
* yt: [0..anzahl] vector, ditto for y                                  *
* lt: Index of final entry in table                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL
***********************************************************************/

{
  REAL x0,
       x1;

  xt += *lt, yt += *lt;
  for (x0 = xa; x0 < xe; x0 += dx)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    x1 = x0 - xi;
    *++xt = x0;
    *++yt = ((((f * x1 + e) * x1 + d) * x1 + c) * x1 + b) * x1 + a;
  }
}


/* ------------------------------------------------------------------ */
/*.BA*/

int hmtab         /* Table of values for a Hermite spline ............*/
/*.IX{hmtab}*/
         (
          int  n,         /* number of spline pieces .................*/
          REAL xanf,      /* left end point for tabulating interval ..*/
          REAL xend,      /* right end point .........................*/
          REAL deltx,     /* step size ...............................*/
          int  anzahl,    /* maximal length of table .................*/
          REAL x[],       /* nodes ...................................*/
          REAL a[],       /* Spline coefficients for (x-x[i])^0 ......*/
          REAL b[],       /* Spline coefficients for (x-x[i])^1 ......*/
          REAL c[],       /* Spline coefficients for (x-x[i])^2 ......*/
          REAL d[],       /* Spline coefficients for (x-x[i])^3 ......*/
          REAL e[],       /* Spline coefficients for (x-x[i])^4 ......*/
          REAL f[],       /* Spline coefficients for (x-x[i])^5 ......*/
          REAL xtab[],    /* x-coordinates of table ..................*/
          REAL ytab[],    /* y-coordinates of table ..................*/
          int  *lentab    /* actual size of table ....................*/
         )                /* error code ..............................*/

/***********************************************************************
* Make a table of values for a  Hermite spline.                        *
* Due to rounding, we may duplicate data near nodes or xend.           *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:           Index of last node                                      *
* xanf:\       end points of interval of table                         *
* xend:/                                                               *
* deltx:       step size for table                                     *
* x:           [0..n] node vector                                      *
* a,b,c,d,e,f: [0..n-1] vectors of spline coefficients.                *
*              a[n] may also be used.                                  *
* anzahl:      upper index limit for xtab and ytab                     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xtab:   [0..anzahl] vector of x-values in table                      *
* ytab:   [0..anzahl] vector, ditto for y                              *
* lentab: Index of final entry in table                                *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: xanf > xend                                                       *
* 2: deltx <= 0                                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hmtabh, REAL, intervall, ZERO, TWO                                   *
.BA*)
***********************************************************************/
/*.BE*/

{
  int anf,    /* Index of spline piece for xanf (0..n-1)              */
      end,    /* Index of spline piece for xend (0..n-1)              */
      anf2,   /* Number of nearest node  <= xanf (-1..n)              */
      end2,   /* Number of nearest node  <= xend (-1..n)              */
      i;      /* Loop variable                                        */

  if (xanf > xend)
    return 1;
  if (deltx <= ZERO)
    return 2;
  anf2 = anf = intervall(n, xanf, x);
  end2 = end = intervall(n, xend, x);
  if (xanf < x[0])  /* xanf to the left of interpolation interval ?  */
    anf2--;
  if (xend > x[n])  /* xend to the right of interpolation interval ? */
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    hmtabh(xanf, x[anf2 + 1], deltx, x[anf], a[anf], b[anf], c[anf],
           d[anf], e[anf], f[anf], xtab, ytab, anzahl, lentab);
    for (i = anf2 + 1; i < end2; i++)
      hmtabh(x[i], x[i + 1], deltx, x[i], a[i], b[i], c[i], d[i], e[i],
             f[i], xtab, ytab, anzahl, lentab);
    xanf = x[end2];                        /* continue table at xanf  */
    if (end2 == n)                         /* record function value at*/
      if (*lentab < anzahl)                /* x[n] if inside interval */
        xtab[++(*lentab)] = x[n],
        ytab[*lentab]     = a[n],
        xanf += deltx;
    hmtabh(xanf, xend, deltx, x[end], a[end], b[end], c[end], d[end],
           e[end], f[end], xtab, ytab, anzahl, lentab);
  }
  else
    hmtabh(xanf, xend, deltx, x[anf], a[anf], b[anf], c[anf], d[anf],
           e[anf], f[anf], xtab, ytab, anzahl, lentab);

  /* ----- special care for right end point xend -------------------- */

  hmtabh(xend, xend + deltx / TWO, deltx, x[end], a[end], b[end],
         c[end], d[end], e[end], f[end], xtab, ytab, anzahl, lentab);

  return 0;
}

/* ------------------------------------------------------------------ */

static void pmtabh      /* aux function for  pmtab() .................*/
/*.IX{pmtabh}*/
                  (
                   REAL ta,
                   REAL te,
                   REAL dt,
                   REAL ti,
                   REAL ax,
                   REAL bx,
                   REAL cx,
                   REAL dx,
                   REAL ex,
                   REAL fx,
                   REAL ay,
                   REAL by,
                   REAL cy,
                   REAL dy,
                   REAL ey,
                   REAL fy,
                   REAL xt[],
                   REAL yt[],
                   int  anzahl,
                   int  *lt
                  )

/***********************************************************************
* Tabulate a parametric Hermite spline inside the parameter interval   *
* from ta to te with step size  dt.                                    *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* ta:                 left end point of interval of table              *
* te:                 right end point                                  *
* dt:                 step size                                        *
* ti:                 left end point of support interval               *
* ax,bx,cx,dx,ex,fx:\ coefficients of spline in support interval       *
* ay,by,cy,dy,ey,fy:/                                                  *
* xt:                 [0..anzahl] vector with old x-values in table    *
* yt:                 [0..anzahl] vector, ditto for y                  *
* anzahl:             upper index limit for xtab and ytab              *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xt: [0..anzahl] vector of x-values in table                          *
* yt: [0..anzahl] vector, y-values                                     *
* lt: Index of final entry in table                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/
{
  REAL t0,
       t;

  xt += *lt, yt += *lt;
  for (t0 = ta; t0 < te; t0 += dt)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    t = t0 - ti;
    *++xt = ((((fx * t + ex) * t + dx) * t + cx) * t + bx) * t + ax;
    *++yt = ((((fy * t + ey) * t + dy) * t + cy) * t + by) * t + ay;
  }
}


/* ------------------------------------------------------------------ */
/*.BA*/

int pmtab   /* Table of values for a parametric Hermite spline .......*/
/*.IX{pmtab}*/
         (
          int  n,         /* number of spline pieces .................*/
          REAL tanf,      /* left end point of interval ..............*/
          REAL tend,      /* right end point .........................*/
          REAL delt,      /* step size ...............................*/
          int  anzahl,    /* maximal size of table ...................*/
          REAL t[],       /* nodes ...................................*/
          REAL ax[],      /* x spline coefficients for (t-t[i])^0 ....*/
          REAL bx[],      /* x spline coefficients for (t-t[i])^1 ....*/
          REAL cx[],      /* x spline coefficients for (t-t[i])^2 ....*/
          REAL dx[],      /* x spline coefficients for (t-t[i])^3 ....*/
          REAL ex[],      /* x spline coefficients for (t-t[i])^4 ....*/
          REAL fx[],      /* x spline coefficients for (t-t[i])^5 ....*/
          REAL ay[],      /* y spline coefficients for (t-t[i])^0 ....*/
          REAL by[],      /* y spline coefficients for (t-t[i])^1 ....*/
          REAL cy[],      /* y spline coefficients for (t-t[i])^2 ....*/
          REAL dy[],      /* y spline coefficients for (t-t[i])^3 ....*/
          REAL ey[],      /* y spline coefficients for (t-t[i])^4 ....*/
          REAL fy[],      /* y spline coefficients for (t-t[i])^5 ....*/
          REAL xtab[],    /* x-coordinates of spline in table ........*/
          REAL ytab[],    /* y-coordinates of spline .................*/
          int  *lentab    /* actual size of table ....................*/
         )                /* error code ..............................*/

/***********************************************************************
* Create a table of values for a parametric Hermite spline.            *
* Due to rounding, we might create duplication near nodes and tend.    *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:                  Index of final parameter node in t               *
* tanf:\              Parameter interval for table                     *
* tend:/                                                               *
* delt:               step size for table.                             *
* t:                  [0..n] vector of parameter nodes                 *
* ax,bx,cx,dx,ex,fx:\ [0..n-1] vectors with splinecoefficients         *
* ax,by,cy,dy,ey,fy:/ ax[n] and ay[n] may also be used                 *
* anzahl:             upper index limit for xtab and ytab              *
*                                                                      *
* Outpur parameters:                                                   *
* ==================                                                   *
* xtab:   [0..anzahl] vector with x-values in table                    *
* ytab:   [0..anzahl] vector, ditto for y                              *
* lentab: Index of final entry in table                                *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: tanf > tend                                                       *
* 2: delt <= 0                                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* pmtabh, REAL, intervall, ZERO, TWO                                   *
.BA*)
***********************************************************************/
/*.BE*/

{
  int anf,    /* Index of spline piece for tanf  (0..n-1)             */
      end,    /* Index of spline piece for tend  (0..n-1)             */
      anf2,   /* Number of nearest node  <= tanf (-1..n)              */
      end2,   /* Number of nearest node  <= tend (-1..n)              */
      i;      /* Loop variable                                        */

  if (tanf > tend)
    return 1;
  if (delt <= ZERO)
    return 2;
  anf2 = anf = intervall(n, tanf, t);
  end2 = end = intervall(n, tend, t);
  if (tanf < t[0])   /* tanf to the left of interpolation interval ?  */
    anf2--;
  if (tend > t[n])   /* tend to the right of interpolation interval ? */
    end2++;
  *lentab = -1;
  if (anf2 < end2)
  {
    pmtabh(tanf, t[anf2 + 1], delt, t[anf], ax[anf], bx[anf],
           cx[anf], dx[anf], ex[anf], fx[anf], ay[anf], by[anf],
           cy[anf], dy[anf], ey[anf], fy[anf], xtab, ytab, anzahl,
           lentab);
    for (i = anf2 + 1; i < end2; i++)
      pmtabh(t[i], t[i + 1], delt, t[i], ax[i], bx[i], cx[i], dx[i],
             ex[i], fx[i], ay[i], by[i], cy[i], dy[i], ey[i], fy[i],
             xtab, ytab, anzahl, lentab);
    tanf = t[end2];                        /* continue table at tanf  */
    if (end2 == n)                         /* record spline value at  */
      if (*lentab < anzahl)                /* t[n] if inside interval */
        xtab[++(*lentab)] = ax[n],
        ytab[*lentab]     = ay[n],
        tanf += delt;
    pmtabh(tanf, tend, delt, t[end], ax[end], bx[end], cx[end],
           dx[end], ex[end], fx[end], ay[end], by[end], cy[end],
           dy[end], ey[end], fy[end], xtab, ytab, anzahl, lentab);
  }
  else
    pmtabh(tanf, tend, delt, t[anf], ax[anf], bx[anf], cx[anf],
           dx[anf], ex[anf], fx[anf], ay[anf], by[anf], cy[anf],
           dy[anf], ey[anf], fy[anf], xtab, ytab, anzahl, lentab);

  /* ----- special care at right end point xend --------------------- */

  pmtabh(tend, tend + delt / TWO, delt, t[end], ax[end], bx[end],
         cx[end], dx[end], ex[end], fx[end], ay[end], by[end],
         cy[end], dy[end], ey[end], fy[end], xtab, ytab, anzahl,
         lentab);

  return 0;
}

/* ------------------------------------------------------------------ */

static void strtabh    /* aux function for  strtab() .................*/
/*.IX{strtabh}*/
                   (
                    REAL pa,
                    REAL pe,
                    REAL dp,
                    REAL pi,
                    REAL a,
                    REAL b,
                    REAL c,
                    REAL d,
                    REAL phid,
                    REAL px,
                    REAL py,
                    REAL xt[],
                    REAL yt[],
                    int  anzahl,
                    int  *lt
                   )

/***********************************************************************
* Tabulate a transformed parametric cubic spline inside the support    *
* interval from  pa to pe with step size  dp.                          *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* pa:      left end point of table interval                            *
* pe:      right end point of same                                     *
* dp:      step size                                                   *
* pi:      left end point of support interval                          *
* a,b,c,d: coefficientsm of spline in support interval                 *
* phid:    angle of rotation of coordinates                            *
* py,py:   translation vector                                          *
* xt:      [0..anzahl] vector with old x-values                        *
* yt:      [0..anzahl] vector with old y-values                        *
* anzahl:  upper index limit for xtab and ytab                         *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* xt: [0..anzahl] vector of old and new x-values                       *
* yt: [0..anzahl] vector, ditto for y                                  *
* lt: Index of final entry in table                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, SIN, COS
***********************************************************************/

{
  REAL p0,
       p1,
       s,
       rho;

  xt += *lt, yt += *lt;
  for (p0 = pa; p0 < pe; p0 += dp)
  {
    if (*lt >= anzahl)
      return;
    (*lt)++;
    p1    = p0 - pi;
    s     = ((d * p1 + c) * p1 + b) * p1 + a;
    rho   = p0 + phid;
    *++xt = s * COS(rho) + px;
    *++yt = s * SIN(rho) + py;
  }
}


/* ------------------------------------------------------------------ */
/*.BA*/

int strtab /* Table of values for transformed parametric cubic spline */
/*.IX{strtab}*/
          (
           int  n,         /* number of spline pieces ................*/
           REAL panf,      /* starting angle for table ...............*/
           REAL pend,      /* final angle of table ...................*/
           REAL phin[],    /* angular nodes ..........................*/
           REAL a[],       /* Spline coeff. for (phi-phin[i])^0 ......*/
           REAL b[],       /* Spline coeff. for (phi-phin[i])^1 ......*/
           REAL c[],       /* Spline coeff. for (phi-phin[i])^2 ......*/
           REAL d[],       /* Spline coeff. for (phi-phin[i])^3 ......*/
           REAL phid,      /* angle of rotation of coordinates .......*/
           REAL px,        /* x-coordinate,                           */
           REAL py,        /* y-coordinate of translation vector .....*/
           REAL x[],       /* nodes: x-values ........................*/
           REAL y[],       /*        y-values ........................*/
           int  nl,        /* maximal length of table ................*/
           int  *nt,       /* actual length of table .................*/
           REAL xtab[],    /* x-coordinates in table .................*/
           REAL ytab[]     /* y-coordinates in table .................*/
          )                /* error code .............................*/

/***********************************************************************
* Make a table of values for a transformed parametric cubic spline :   *
*   s(phi) = a[i] + b[i](phi-phin[i]) + c[i](phi-phin[i])^2 +          *
*                                     + d[i](phi-phin[i])^3            *
* for phi in  [phin[i], phin[i+1]], i=0, ..., n-1.                     *
* The table contains the values                                        *
*        xtab = xtab(phi) = s(phi) * cos(phi + phid) + px,             *
*        ytab = ytab(phi) = s(phi) * sin(phi + phid) + py,             *
* with phi in  [panf, pend].                                           *
.BE*)
* INTERPRETATIONS:                                                     *
*   - If panf < phin[0], the end point polynomial p[0] is evaluated    *
*     for all values of  xtab(phi) with  phi < phin[0].                *
*   - If pend > phin[n], the end point polynomial p[n-1] is evaluated  *
*     for all  xtab(phi) with  phi > phin[n].                          *
*   - the interval end points panf and pend and all nodes in between   *
*     are entered into the table.                                      *
*   - In each sub-interval  [phin[i],phin[i+1]] the table is formed    *
*     for equidistant steps of length h, wherei h depends on the       *
*     individual length of the interval and the maximal table length nl*
*   - The parameter nl denotes the approximate table size; its actual  *
*     size is nt + 1. (nt is the final entry index)  And in general :  *
*                     0 < nt < nl + n + 1.                             *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n:       Index of final node                                         *
* panf:\   Parameter interval for table                                *
* pend:/                                                               *
* phin:    [0..n] vector of parameter nodes                            *
* a,b,c,d: [0..n-1] vectors with spline coefficients                   *
* phid:    angle of rotation                                           *
* py,py:   translation vector                                          *
* nl:      upper index limit on xtab and ytab                          *
* x,y:     [0..n] vectors with the original nodes, which are used to   *
*          find the table value at the parameter nodes phin[i].        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* nt:   Index of final entry in table                                  *
* xtab: [0..nl] vector of x-values in table                            *
* ytab: [0..nl] vector, ditto for y                                    *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: panf >= pend                                                      *
* 2: n < 1                                                             *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* strtabh, REAL, intervall, TWO                                        *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  anf,   /* Index of spline piece for panf (0..n-1)              */
       end,   /* Index of spline piece for pend (0..n-1)              */
       anf2,  /* Number of nearest node  <= panf (-1..n)              */
       end2,  /* Number of nearest node  <= pend (-1..n)              */
       i;     /* Loop variable                                        */
  REAL h;     /* step size                                            */

  if (pend <= panf)
    return 1;
  if (n < 1)
    return 2;

  anf2 = anf = intervall(n, panf, phin);
  end2 = end = intervall(n, pend, phin);
  if (panf < phin[0])/* panf to the left of interpolation interval ?  */
    anf2--;
  if (pend > phin[n])/* pend to the right of interpolation interval ? */
    end2++;
                                      /* choose step size so that all */
  h = (pend - panf) / (nl - n - 1);   /* intermediate nodes fit into  */
  *nt = -1;                           /* the table                    */
  if (anf2 < end2)                    /* nodes inside table interval? */
  {
    if (panf == phin[anf])
      if (*nt < nl)                                /* record node     */
        xtab[++(*nt)] = x[anf],                    /* (x[anf],y[anf]  */
        ytab[*nt]     = y[anf],
        panf += h;
    strtabh(panf, phin[anf2 + 1], h,         /* compute table values  */
            phin[anf], a[anf],               /* from panf to          */
            b[anf], c[anf], d[anf],          /* phin[anf2+1]          */
            phid, px, py, xtab, ytab,
            nl, nt);

    for (i = anf2 + 1; i < end2; i++)        /* compute table entries */
    {                                        /* from phin[anf2+1] to  */
      panf = phin[i];                        /* phin[end2]            */
      if (*nt < nl)
        xtab[++(*nt)] = x[i],                    /* record the node   */
        ytab[*nt]     = y[i],                    /* (x[i],y[i])       */
        panf += h;
      strtabh(panf, phin[i + 1], h,
              phin[i], a[i], b[i],
              c[i], d[i], phid, px,
              py, xtab, ytab, nl, nt);
    }

    panf = phin[end2];
    if (*nt < nl)
      xtab[++(*nt)] = x[end2],                   /* record the node   */
      ytab[*nt]     = y[end2],                   /* (x[end2],y[end2]) */
      panf += h;
    strtabh(panf, pend, h, phin[end],        /* compute the values    */
            a[end], b[end], c[end],          /* from phin[end2] to    */
            d[end], phid, px, py,            /* pend                  */
            xtab, ytab, nl, nt);
  }
  else                        /* no node inside tabulating interval ? */
    strtabh(panf, pend, h, phin[anf], a[anf], b[anf], c[anf],
            d[anf], phid, px, py, xtab, ytab, nl, nt);

  /* ----- special care for right end point xend -------------------- */

  strtabh(pend, pend + h / TWO, h, phin[end], a[end], b[end],
          c[end], d[end], phid, px, py, xtab, ytab, nl, nt);

  return 0;
}

/* -------------------------- END splintab.c ------------------------ */
