#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE zeigkrv2.c ------------------------ */

/***********************************************************************
*                                                                      *
* Plotting curves and surfaces                                         *
* ----------------------------                                         *
*                                                                      *
* exported functions:                                                  *
*   - zeigkrv2():     graphical representation of the tabulated values *
*                     of a planar curve and of the nodes that were     *
*                     used to construct the curve                      *
*   - zeigflaeche():  graphical representation of the tabulated values *
*                     of a surface in three-dimensional space and of   *
*                     the nodes that were used to construct the        *
*                     surface                                          *
*                                                                      *
* Programming language: Borland C++ 2.0, QuickC 2.0                    *
* Compiler:             Borland C++ 2.0, PureC 1.0, QuickC 2.0         *
* Author:               Juergen Dietel, Computer Center, RWTH Aachen   *
* Date:                 12.22.1993 - 2.26.1997                         *
*                                                                      *
***********************************************************************/

#include <cnumgraf.h>     /*  for  MITGRAFIK, BGIGRAFMOEGLICH,        */
                          /*       MCGRAFMOEGLICH, GRAFIK_H,          */
                          /*       BIOS_H, INIT_FONTS                 */

#ifdef MITGRAFIK          /* graphics supported or wanted?            */

#if defined(BGIGRAFMOEGLICH)

#include GRAFIK_H         /*  for  initgraph, closegraph, grOk,       */
                          /*       DETECT, outtextxy, setcolor,       */
                          /*       graphresult, getmaxcolor,          */
                          /*       SOLID_LINE, THICK_WIDTH, YELLOW,   */
                          /*       WHITE, moveto, LIGHTGREEN,         */
                          /*       LIGHTMAGENTA, lineto,              */
                          /*       setlinestyle, outtext, LIGHTRED,   */
                          /*       NORM_WIDTH, setviewport,           */
                          /*       graphdefaults, line, rectangle,    */
                          /*       textheight, getmaxx, getmaxy,      */
                          /*       LIGHTCYAN, DASHED_LINE, drawpoly   */
#include <stdlib.h>       /*  for  getenv, itoa                       */
#include <limits.h>       /*  for  INT_MAX                            */
#include BIOS_H           /*  for  bioskey                            */

#elif defined(MCGRAFMOEGLICH)

#include <graph.h>        /*  for  _GBORDER, _rectangle, videoconfig, */
                          /*       _VRES16COLOR, _setvideomode,       */
                          /*       _TEXTMONO, _getvideoconfig,        */
                          /*       _setcolor, _outtext,               */
                          /*       _settextposition, _setvieworg,     */
                          /*       _moveto, _lineto, _bios_keybrd,    */
                          /*       _KEYBRD_READ, _DEFAULTMODE         */
#include <stdlib.h>       /*  for  itoa                               */
#include <bios.h>         /*  for  _bios_keybrd, _KEYBRD_READ         */
#endif                               /* #elif defined(MCGRAFMOEGLICH) */

#include <basis.h>        /*  for  REAL, ONE, TWO, max, THREE         */
#include <vmblock.h>      /*  for  vminit, vmalloc, IMATRIX,          */
                          /*       vmcomplete, vmfree, VVEKTOR        */
#include <zeigkrv2.h>     /*  for  zeigkrv2, zeigflaeche              */



#ifdef MCGRAFMOEGLICH
/*--------------------------------------------------------------------*/

#define LIGHTGREEN    10
#define LIGHTCYAN     11
#define LIGHTRED      12
#define LIGHTMAGENTA  13
#define YELLOW        14
#define WHITE         15



#endif                                       /* #ifdef MCGRAFMOEGLICH */
#ifdef GRAFIKMOEGLICH
/*--------------------------------------------------------------------*/

static void getminmax
        (
         int  nk,
         REAL *kurve[],
         int  ns,
         REAL *stuetz[],
         REAL *xmin,
         REAL *xmax,
         REAL *ymin,
         REAL *ymax
        )

/***********************************************************************
* Find the minimal rectangle [xmin,xmax] x [ymin,ymax] that contains   *
* the curve points (kurve[i][0],kurve[i][1]), i=0,1,...,nk-1, and      *
* the nodes (stuetz[i][0],stuetz[i][1]), i=0,1,...,ns-1.               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  REAL mix, max, miy, may;     /* aux variables for minima and maxima */


  for (mix = max = kurve[0][0], miy = may = kurve[0][1], nk--, kurve++;
       nk != 0;
       nk--, kurve++
      )
  {
    if (kurve[0][0] < mix)
      mix = kurve[0][0];
    else if (kurve[0][0] > max)
      max = kurve[0][0];
    if (kurve[0][1] < miy)
      miy = kurve[0][1];
    else if (kurve[0][1] > may)
      may = kurve[0][1];
  }

  for ( ; ns != 0; ns--, stuetz++)
  {
    if (stuetz[0][0] < mix)
      mix = stuetz[0][0];
    else if (stuetz[0][0] > max)
      max = stuetz[0][0];
    if (stuetz[0][1] < miy)
      miy = stuetz[0][1];
    else if (stuetz[0][1] > may)
      may = stuetz[0][1];
  }

  *xmin = mix;
  *xmax = max;
  *ymin = miy;
  *ymax = may;
}



/*--------------------------------------------------------------------*/

static void transformieren
        (
         int  n,
         REAL *punkt[],
         REAL xmin,
         REAL xmax,
         REAL ymin,
         REAL ymax,
         int  grafmaxx,
         int  grafmaxy,
         int  polypoints[]
        )

/***********************************************************************
* Transform the points  (punkt[i][0],punkt[i][1]), i=0,...,n-1 linearly*
* from the rectangle  [xmin,xmax] x [ymin,ymax]  to the plot window    *
* [0,grafmaxx] x [0,grafmaxy] so that xmin maps to 0, xmax to grafmaxx,*
* ymin to grafmaxy and ymax to  0.                                     *
* The resulting truncated integer values (in point scale for plotting) *
* are stored in the vector polypoints which contains an alternating    *
* sequence of x- and  y-coordinates. We flip the y-coordinates to      *
* avoid an upside down picture.                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  REAL fx,                  /* slope of linear transform for x-values */
       fy;                  /* ditto for y                            */

  for (fx = grafmaxx / (xmax - xmin), fy = grafmaxy / (ymax - ymin);
       n != 0;
       n--, punkt++)
    *polypoints++ = (int)           ((punkt[0][0] - xmin) * fx),
    *polypoints++ = (int)(grafmaxy - (punkt[0][1] - ymin) * fy);
}



#endif                                      /* #ifdef BGIGRAFMOEGLICH */
#ifdef BGIGRAFMOEGLICH
/*--------------------------------------------------------------------*/

static void mein_setcolor
        (
         int  farbe
        )

/***********************************************************************
* set up a new color for plotting, but limit their maximal number.     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* setcolor, getmaxcolor                                                *
***********************************************************************/

{
  int maxfarbe;

  setcolor((farbe > (maxfarbe = getmaxcolor())) ? maxfarbe : farbe);
}



/*--------------------------------------------------------------------*/

static void stuetzpunkte_eintragen
        (
         int  ns,
         int  polypoints[]
        )

/***********************************************************************
* Mark the screen points (polypoints[2*i],polypoints[2*i+1]) from      *
* the [0..2*ns-1] vector polypoints by a cross, and additionally mark  *
* the first point in polypoints by a rectangle.                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* mein_setcolor, LIGHTGREEN, setlinestyle, SOLID_LINE, THICK_WIDTH,    *
* line, LIGHTRED, NORM_WIDTH, rectangle                                *
***********************************************************************/

{
  int i,               /* Loop variable                               */
      x, y;            /* screen coordinates of current point         */


  mein_setcolor(LIGHTGREEN);                        /* crosses thick  */
  setlinestyle(SOLID_LINE, 0, THICK_WIDTH);         /* in light green */


  for (i = 0; i < ns; i++)
  {
    x = *polypoints++;
    y = *polypoints++;
    line(x - 4, y - 4, x + 4, y + 4);
    line(x + 4, y - 4, x - 4, y + 4);

    if (i == 0)                                     /* first point?   */
    {                                               /* mark with a    */
      mein_setcolor(LIGHTRED);                      /* red rectangle  */
      setlinestyle(SOLID_LINE, 0, NORM_WIDTH);
      rectangle(x - 6, y - 6, x + 6, y + 6);
      mein_setcolor(LIGHTGREEN);
      setlinestyle(SOLID_LINE, 0, THICK_WIDTH);
    }
  }
}



#elif defined(MCGRAFMOEGLICH)
/*--------------------------------------------------------------------*/

static void stuetzpunkte_eintragen
        (
         int  ns,
         int  polypoints[]
        )

/***********************************************************************
* Mark the screen points (polypoints[2*i],polypoints[2*i+1]) from      *
* the [0..2*ns-1] vector polypoints by a cross, and additionally mark  *
* the first point in polypoints by a rectangle.                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* _setcolor, LIGHTGREEN, _moveto, _lineto, LIGHTRED, _GBORDER,         *
* _rectangle,                                                          *
***********************************************************************/

{
  int i,               /* Loop variable                               */
      x, y;            /* screen coordinates of current point         */


  _setcolor(LIGHTGREEN);                       /* light green crosses */


  for (i = 0; i < ns; i++)
  {
    x = *polypoints++;
    y = *polypoints++;
    _moveto(x - 4, y - 4);
    _lineto(x + 4, y + 4);
    _moveto(x + 4, y - 4);
    _lineto(x - 4, y + 4);

    if (i == 0)                                     /* first point?   */
    {                                               /* mark by a red  */
      _setcolor(LIGHTRED);                          /* rectangle      */
      _rectangle(_GBORDER, x - 6, y - 6,
                           x + 6, y + 6);
      _setcolor(LIGHTGREEN);
    }
  }
}



/*--------------------------------------------------------------------*/

static void drawpoly
        (
         int  n,
         int  polypoints[]
        )

/***********************************************************************
* Draw a polygon in the current line type and color.                   *
* Here n denotes the number of corners and polypoints is the vector of *
* n pairs of coordinates of those points.                              *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* _moveto, _lineto                                                     *
***********************************************************************/

{
  _moveto(polypoints[0], polypoints[1]);
  for (; n != 0; n--, polypoints += 2)
    _lineto(polypoints[0], polypoints[1]);
}
#endif                               /* #elif defined(MCGRAFMOEGLICH) */



#ifdef BGIGRAFMOEGLICH
/*--------------------------------------------------------------------*/

int zeigkrv2         /* plot a function table in R2 ..................*/
        (
         int  nk,                  /* size of function table .........*/
         REAL *kurve[],            /* function table .................*/
         int  ns,                  /* number of nodes ................*/
         REAL *stuetz[]            /* nodes ..........................*/
        )                          /* error code .....................*/

/***********************************************************************
* The points (kurve[i][0],kurve[i][1]), i=0,...,nk-1, form a table of  *
* values of a planar curve. They define a polygonal line which shall   *
* be plotted so that the screen is maximally filled.                   *
* To achieve this we find the extrema of both sets of coordinates      *
* (xmin, xmax, ymin, ymax) and then transform the rectangle            *
* [xmin,xmax] x [ymin,ymax]  into the plot window. There the points    *
* are represented in integer point scale coordinates for drawing by    *
* drawpoly().                                                          *
*                                                                      *
* The [0..ns-1] vector stuetz contains the nodes used to generate the  *
* table of values for the curve such as the nodes of a cubic spline or *
* the de Boor points for B-splines etc. For color monitors these will  *
* will be plotted by a green cross. Clearly the extrema search above   *
* is also extended over this node data.                                *
*                                                                      *
* This program is designed for the BGI (Borland Graphics Interface).   *
* On computers or C compilers that don't support BGI one should use    *
* analoguous functions for initgraph(), closegraph() etc.              *
*                                                                      *
* The graphics drivers that are loaded via initgraph() are expected to *
* be available in the directory whose name - prior to the start of the *
* program - was put into the environment variable BGI (e. g. by using  *
* the MS-DOS command `set BGI=c:\bc\bgi', if the directory  C:\BC\BGI  *
* contains the needed graphics drivers *.BGI).                         *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* 0: all ok                                                            *
* 1: nk < 1                                                            *
* 2: ns < 1                                                            *
* 3: lack of memory                                                    *
* 4: Graphics error                                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, DETECT, initgraph, getenv, getminmax, max, ONE, TWO,           *
* textheight, getmaxx, getmaxy, itoa, mein_setcolor, YELLOW,           *
* outtextxy, moveto, outtext, setviewport, LIGHTCYAN, lineto,          *
* transformieren, WHITE, drawpoly, stuetzpunkte_eintragen,             *
* graphdefaults, LIGHTMAGENTA, bioskey, closegraph, graphresult, grOK, *
* vminit, vmalloc, VVEKTOR, vmcomplete, vmfree, INIT_FONTS             *
***********************************************************************/

{
  REAL  xmin, xmax,          /* extreme values of x- and y-coordinates*/
        ymin, ymax;
  int   fensterx, fenstery,  /* left upper corner of plot             */
        fensterb, fensterh,  /* width and height of plot window       */
        *polypoints,         /* Vector of screen coordinates          */
        grafiktreiber,       /* Number of graphics driver             */
        grafikmodus,         /* modus set by initgraph() for graphics */
                             /* driver                                */
        grafikfehler,        /* error code from BGI                   */
        texthoehe;           /* height of X in points                 */
  char  ns_string[6];        /* ns as an alpha-numeric string         */
  void  *vmblock;            /* list of dynamical vectors             */


  /* -------------- check input errors ------------------------------ */

  if (nk < 1)
    return 1;

  if (ns < 1)
    return 2;


  /* ---- allocate the vector of screen coordinates dynamically  ---- */
  /* ---- of sufficient size to be able to accomodate all points ---- */
  /* ---- on the curve as well as all nodes after transformation ---- */

  vmblock    = vminit();
  polypoints = (int *)vmalloc(vmblock, VVEKTOR,
                              max(nk, ns), 2 * sizeof(*polypoints));
  if (! vmcomplete(vmblock))                     /* lack of memory?   */
    return 3;                                    /* nothing goes!     */


  /* ------------------- initialize graphics ------------------------ */

  grafiktreiber = DETECT;          /* find the optimal graphics driver*/
                                   /* in initgraph()                  */
  initgraph(&grafiktreiber, &grafikmodus, getenv("BGI"));
  if ((grafikfehler = graphresult()) != grOk)      /* not successful? */
  {
    vmfree(vmblock);
    return 4;
  }

  INIT_FONTS;

  texthoehe = textheight("X");             /* Height of text          */

  fensterx = 8;                            /* define plot window      */
  fensterb = getmaxx() - 8 - fensterx;
  fenstery = 2 * texthoehe + 14;
  fensterh = getmaxy() - 8 - fenstery;


  /* ----------- find minimal rectangle that contains all  ---------- */
  /* ----------- nodes and points on the curve             ---------- */

  getminmax(nk, kurve,                     /* find extreme x- and y-  */
            ns, stuetz,                    /* coordinates of points   */
            &xmin, &xmax,                  /* on curve and nodes      */
            &ymin, &ymax
           );
  if (xmin == xmax)                        /* x-interval zero?        */
    xmin = kurve[0][0] - ONE,              /* set to 2 for transform  */
    xmax = xmin        + TWO;
  if (ymin == ymax)                        /* y-interval zero?        */
    ymin = kurve[0][1] - ONE,              /* set to 2 for transform  */
    ymax = ymin        + TWO;


  /* ------------------- make heading ------------------------------- */

  itoa(ns, ns_string, 10);
  mein_setcolor(YELLOW);
  outtextxy(0, 0, "A planar curve given by");
  moveto(0, texthoehe + 3);
  outtext(ns_string);
  outtext(" nodes marked by crosses");


  /* ------------ set up plot window and draw ----------------------- */

  setviewport(fensterx, fenstery,   /* shift screen origin from (0,0) */
              fensterx + fensterb,  /* to (fensterx,fenstery)         */
              fenstery + fensterh,
              0);

  mein_setcolor(LIGHTCYAN);             /* draw a frame around window */
  moveto(-8,           -8);             /* in light cyan (8 points    */
  lineto(fensterb + 8, -8);             /* from proper window)        */
  lineto(fensterb + 8, fensterh + 8);
  lineto(-8,           fensterh + 8);
  lineto(-8,           -8);


  /* ------------- draw points on curve and nodes ------------------- */

  transformieren(nk, kurve,                 /* transform curve points */
                 xmin, xmax,                /* into the plot window   */
                 ymin, ymax,                /* and draw polygon in    */
                 fensterb, fensterh,        /* white                  */
                 polypoints
                );
  mein_setcolor(WHITE);
  drawpoly(nk, polypoints);

  transformieren(ns, stuetz,                /* transform the nodes    */
                 xmin, xmax,                /* into plot window       */
                 ymin, ymax,                /* and draw               */
                 fensterb, fensterh,
                 polypoints
                );
  stuetzpunkte_eintragen(ns, polypoints);


  /* ------------ finish plot --------------------------------------- */

  graphdefaults();                        /* reverse coordinate shift */

  mein_setcolor(LIGHTMAGENTA);            /* print input demand in    */
  outtextxy(getmaxx() - 200, 0,           /* light magenta            */
            "continue with any key...");
  bioskey(0);                             /* wait for key stroke      */
  closegraph();                           /* finish work with the     */
                                          /* graphics driver          */

  if ((grafikfehler = graphresult()) != grOk)      /* Graphics error? */
  {
    vmfree(vmblock);
    return 4;
  }


  grafikfehler = grafikfehler;            /* calm compiler            */
  vmfree(vmblock);                        /* free dynamic allocations */
  return 0;
}



/* ------------------------------------------------------------------ */

/***********************************************************************
*                    3                                                 *
* To plot points in R  on a two-dimensional surface, one needs a       *
*                      3       2                                       *
* transformation from R  into R , which is chosen here as parallel     *
* projection along the vector p = (px,py,pz) onto the x-z-plane.       *
* The corresponding transformation formulas are:                       *
*                                                                      *
*     x  -->  x - y * (px /py)                                         *
*     y  -->  0                                                        *
*     z  -->  z - y * (pz /py)                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ONE, THREE                                                     *
***********************************************************************/

struct { REAL x, y, z; } p = { ONE, -THREE, ONE };



/* ------------------------------------------------------------------ */

static void projizieren     /* project point matrix onto x-z-plane ...*/
        (
         REAL **c[],               /* point matrix ...................*/
         int  nv,                  /* number of rows .................*/
         int  nw                   /* number of columns ..............*/
        )

/***********************************************************************
* project the points (c[i,j,0],c[i,j,1],c[i,j,2]), i=0(1)nv-1,         *
* j = 0(1)nw-1, along p onto the x-z-plane.                            *
* Doing this the y component of the points remains unchanged, as it is *
* not needed any more.                                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, p                                                              *
***********************************************************************/

{
  REAL faktorx,
       faktorz;
  int  i, j;


  faktorx = p.x / p.y;
  faktorz = p.z / p.y;

  for (i = 0; i < nv; i++)
    for (j = 0; j < nw; j++)
      c[i][j][0] -= c[i][j][1] * faktorx,
      c[i][j][2] -= c[i][j][1] * faktorz;
}



/* ------------------------------------------------------------------ */

static void fminmax   /* find enclosing hexahedron of points in R3 ...*/
        (
         REAL **c[],               /* first matrix of points .........*/
         int  nv,                  /* it is a ........................*/
         int  nw,                  /* [0..nv-1,0..nw-1] matrix. ......*/
         REAL **d[],               /* second matrix of points ........*/
         int  m,                   /* it is a ........................*/
         int  n,                   /* [0..m-1,0..n-1] matrix. ........*/
         REAL *r3xmin,             /* enclosing intervals ............*/
         REAL *r3xmax,             /* [r3xmin,r3xmax] x  .............*/
         REAL *r3ymin,             /* [r3zmin,r3zmax] x  .............*/
         REAL *r3ymax,             /* [r3zmin,r3zmax]    .............*/
         REAL *r3zmin,
         REAL *r3zmax
        )

/***********************************************************************
* find the enclosing hexahedron                                        *
*         [r3xmin,r3xmax] x [r3ymin,r3ymax] x [r3zmin,r3zmax]          *
* that is spanned by the points (c[i,j,0],c[i,j,1],c[i,j,2]),          *
* i=0(1)nv-1, j=0(1)nw-1, and (d[i,j,0],d[i,j,1],d[i,j,2]), i=0(1)m-1, *
* j=0(1)n-1                                                            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  REAL mix, max,             /* aux variables for minima und maxima   */
       miy, may,
       miz, maz;
  int  i, j;                 /* loop variables                        */


  mix = c[0][0][0];
  max = mix;
  miy = c[0][0][1];
  may = miy;
  miz = c[0][0][2];
  maz = miz;
  for (i = 0; i < nv; i++)
    for (j = 0; j < nw; j++)
    {
      if (c[i][j][0] < mix)
        mix = c[i][j][0];
      else if (c[i][j][0] > max)
        max = c[i][j][0];
      if (c[i][j][1] < miy)
        miy = c[i][j][1];
      else if (c[i][j][1] > may)
        may = c[i][j][1];
      if (c[i][j][2] < miz)
        miz = c[i][j][2];
      else if (c[i][j][2] > maz)
        maz = c[i][j][2];
    }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      if (d[i][j][0] < mix)
        mix = d[i][j][0];
      else if (d[i][j][0] > max)
        max = d[i][j][0];
      if (d[i][j][1] < miy)
        miy = d[i][j][1];
      else if (d[i][j][1] > may)
        may = d[i][j][1];
      if (d[i][j][2] < miz)
        miz = d[i][j][2];
      else if (d[i][j][2] > maz)
        maz = d[i][j][2];
    }

  *r3xmin = mix;
  *r3xmax = max;
  *r3ymin = miy;
  *r3ymax = may;
  *r3zmin = miz;
  *r3zmax = maz;
}



/* ------------------------------------------------------------------ */

static void ftrafo  /* transform planar points to screen coordinates  */
        (
         REAL **c[],               /* matrix of points ...............*/
         int  nv,                  /* it is a ........................*/
         int  nw,                  /* [0..nv-1,0..nw-1] matrix. ......*/
         REAL xmin,                /* starting rectangle .............*/
         REAL xmax,                /* [xmin,xmax]x[zmin,zmax] ........*/
         REAL zmin,
         REAL zmax,
         int  grafmaxx,            /* destination rectangle ..........*/
         int  grafmaxy,            /* [0,grafmaxx]x[0,grafmaxy] ......*/
         int  *polypoints[]        /* points for screen ..............*/
        )

/***********************************************************************
* transform the points (c[i,j,0],c[i,j,2]), i=0(1)nv-1, j=0(1)nw-1,    *
* linearly from the rectangle  [xmin,xmax] x [zmin,zmax]  to the plot  *
* window  [0,grafmaxx] x [0,grafmaxy]  so that xmin maps to 0, xmax    *
* to grafmaxx, zmin to grafmaxy and zmax to 0.                         *
* The resulting truncated integer values (in point scale for plotting) *
* are stored in the matrix polypoints which contains an alternating    *
* sequence of x- and y-coordinates. We flip the y-coordinates to avoid *
* an upside down picture.                                              *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  REAL fx,                  /* slope of linear transform for x-values */
       fz;                  /* slope of linear transform for z-values */
  int  i, j;                /* loop variables                         */


  fx = grafmaxx / (xmax - xmin);
  fz = grafmaxy / (zmax - zmin);

  for (i = 0; i < nv; i++)
    for (j = 0; j < nw; j++)
      polypoints[i][2 * j]     =            (c[i][j][0] - xmin) * fx,
      polypoints[i][2 * j + 1] = grafmaxy - (c[i][j][2] - zmin) * fz;
}



/* ------------------------------------------------------------------ */

static void mal_stuetz   /* draw polygons pointwise and as grid ......*/
        (
         int  **polypoints,        /* matrix of polygons .............*/
         int  m,                   /* number of points: ..............*/
         int  n,                   /* m * n  .........................*/
         int  voffen,              /* open v curves? .................*/
         int  woffen,              /* open w curves? .................*/
         int  st_gitter            /* with grid? .....................*/
        )

/***********************************************************************
* mark the screen points (polypoints[i,2*j],polypoints[i,2*j+1]) from  *
* the [0..m-1,0..n-1] matrix polypoints by crosses and draw the points *
* both in v direction (first index of polypoints) and in w direction   *
* (second index of polypoints) as polygons. Additionally the points of *
* the first node vector (those belonging to index i = 0), which go in  *
* w direction, are marked by rectangles.                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* moveto, mein_setcolor, LIGHTGREEN, setlinestyle, SOLID_LINE,         *
* THICK_WIDTH, line, YELLOW, NORM_WIDTH, lineto, LIGHTRED, rectangle   *
***********************************************************************/

{
  int x, y,                    /* screen coordinates of current point */
      i, j;                    /* loop variables                      */


  for (i = 0; i < m; i++)                     /* at first: crosses    */
  {                                           /* and w curves         */
    moveto(polypoints[i][0], polypoints[i][1]);
    for (j = 0; j < n; j++)
    {
      x = polypoints[i][2 * j];
      y = polypoints[i][2 * j + 1];
      mein_setcolor(LIGHTGREEN);
      setlinestyle(SOLID_LINE, 0, THICK_WIDTH);
      line(x - 4, y - 4, x + 4, y + 4);
      line(x + 4, y - 4, x - 4, y + 4);
      setlinestyle(SOLID_LINE, 0, NORM_WIDTH);
      if (st_gitter)
      {
        mein_setcolor(YELLOW);
        lineto(x, y);
      }
      if (i == 0)                              /* first               */
      {                                        /* node vector?        */
        mein_setcolor(LIGHTRED);               /* mark points         */
        rectangle(x - 6, y - 6, x + 6, y + 6); /* additionally by     */
        mein_setcolor(LIGHTGREEN);             /* red rectangles      */
      }
    }
    setlinestyle(SOLID_LINE, 0, NORM_WIDTH);
    if (! woffen)                             /* closed w curves?     */
      lineto(polypoints[i][0],                /* => back to the start */
             polypoints[i][1]);               /*    of the w curve    */
  }

  setlinestyle(SOLID_LINE, 0, NORM_WIDTH);    /* and now the          */
  if (st_gitter)                              /* v curves             */
    for (j = 0; j < n; j++)
    {
      moveto(polypoints[0][2 * j],
             polypoints[0][2 * j + 1]);
      for (i = 0; i < m; i++)
        lineto(polypoints[i][2 * j],
               polypoints[i][2 * j + 1]);
      if (! voffen)                           /* closed v curves?     */
        lineto(polypoints[0][2 * j],          /* => back to the start */
               polypoints[0][2 * j + 1]);     /*    of the v curve    */
    }
}



/* ------------------------------------------------------------------ */

int zeigflaeche      /* plot a function table in R3 ..................*/
        (
         REAL **c[],               /* function table .................*/
         int  nv,                  /* The table is a                  */
         int  nw,                  /* [0..nv-1,0..nw-1] matrix. ......*/
         REAL **d[],               /* nodes ..........................*/
         int  m,                   /* The nodes form a                */
         int  n,                   /* [0..m-1,0..n-1] matrix. ........*/
         int  voffen,              /* open v curves? .................*/
         int  woffen,              /* open w curves? .................*/
         int  st_gitter            /* grid through nodes? ............*/
        )                          /* error code .....................*/

/***********************************************************************
* The points (c[i,j,0],c[i,j,1],c[i,j,2]), i=0(1)nv-1, j=0(1)nw-1,     *
*                                                                  3   *
* form a table of values of a surface in three-dimensional space (R ). *
* They define a grid which shall be plotted so that the screen is      *
* maximally filled.
* To achieve this we first project all points (both surface points and *
* node points) onto a fictitious x-z-drawing-plane. Now we find the    *
* extrema (xmin, xmax, zmin, zmax) of the two remaining coordinates of *
* the projected points and then transform them linearly from the       *
* rectangle  [xmin,xmax] x [zmin,zmax]  into the plot window so that   *
* they can be stored in the integer matrix polypoints. This was        *
* defined in such a way that you just have to pass one of its rows to  *
* the function drawpoly() to get a graphical representation of the     *
* corresponding curve on the surface. Thus we receive two grids on the *
* screen, one from the surface points and the other from the node      *
* points (d[i,j,0],d[i,j,1],d[i,j,2]), i=0(1)m-1, j=0(1)n-1 which were *
* used to construct the surface (for instance like the de Boor points  *
* in B splines). Addtionally the node points are marked by crosses to  *
* emphasize the difference between the two grids.                      *
*                                                                      *
* Caution: After the call of this function c and d no longer contain   *
*          the original points, but their projections onto the         *
*          x-z-plane.                                                  *
*                                                                      *
* The input parameters voffen and woffen are needed to be able to draw *
* the polygons defined by the node points open or closed in            *
* correspondence with the v and w curves defined by the surface        *
* points.                                                              *
*                                                                      *
* This solution is designed for the BGI (Borland Graphics Interface).  *
* On computers or C compilers that don't support BGI one should use    *
* analoguous functions for initgraph(), closegraph() etc.              *
*                                                                      *
* The graphics drivers that are loaded via initgraph() are expected to *
* be available in the directory whose name - prior to the start of the *
* program - was put into the environment variable BGI (e. g. by using  *
* the MS-DOS command `set BGI=c:\bc\bgi', if the directory  C:\BC\BGI  *
* contains the needed graphics drivers *.BGI).                         *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* 0: all ok                                                            *
* 1: nv < 1                                                            *
* 2: nw < 1                                                            *
* 3: m  < 1                                                            *
* 4: n  < 1                                                            *
* 5: nw or n to large for memory allocations                           *
* 6: graphics error                                                    *
* 7: lack of memory                                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, vminit, vmalloc, IMATRIX, vmcomplete, vmfree, max, DETECT,     *
* initgraph, getenv, graphresult, grOk, textheight, getmaxx, getmaxy,  *
* fminmax, projizieren, ONE, TWO, itoa, mein_setcolor, YELLOW,         *
* outtextxy, moveto, outtext, setviewport, LIGHTCYAN, lineto, WHITE,   *
* setlinestyle, DASHED_LINE, NORM_WIDTH, SOLID_LINE, ftrafo, drawpoly, *
* mal_stuetz, graphdefaults, LIGHTMAGENTA, bioskey, closegraph,        *
* INT_MAX                                                              *
***********************************************************************/

{
  REAL r3xmin,             /* coordinates of the minimal hexahedron   */
                           /* enclosing all surface and node points   */
       r3xmax,
       r3ymin,
       r3ymax,
       r3zmin,
       r3zmax,
       xmin, xmax,         /* minimal and maximal values of all       */
       zmin, zmax,         /* projected surface and node points       */
       faktor;             /* needed for computing diffx, diffy       */
  int  **polypoints,       /* matrix of pairs of screen coordinates;  */
                           /* its first index corresponds to the      */
                           /* v direction of the surface, the second  */
                           /* index to the w direction                */
       fensterx,           /* left upper corner, width and height     */
       fenstery,           /* of the screen window where surface and  */
       fensterb,           /* node points are drawn                   */
       fensterh,           /* (diagram window, plot window)           */
       grafiktreiber,      /* number of the graphics driver           */
       grafikmodus,        /* modus set by initgraph() for graphics   */
                           /* driver                                  */
       grafikfehler,       /* error code from BGI                     */
       texthoehe,          /* height of an 'X' in the standard font   */
       diffx,              /* (diffx,diffy) = vector by which the     */
       diffy,              /* rear plane of the hexahedron in the     */
                           /* plot window seems to shifted against    */
                           /* the front plane                         */
       i, j;               /* loop variables                          */
  char mn_string[6];       /* m*n as an decimal alpha-numeric string  */
  void *vmblock;           /* list of dynamical vectors/matrices      */


  /* ------------------ check errors in input data ------------------ */

  if (nv < 1)
    return 1;

  if (nw < 1)
    return 2;

  if (m < 1)
    return 3;

  if (n < 1)
    return 4;

  if (max(nw, n) > INT_MAX / 2)   /* danger of overflow in vmalloc()? */
    return 5;

  /* --- allocate the matrix of screen coordinates dynamically    --- */
  /* --- of sufficient size to be able to accomodate all points   --- */
  /* --- on the surface as well as all nodes after transformation --- */

  vmblock    = vminit();
  polypoints = (int **)vmalloc(vmblock, IMATRIX,
                               max(nv, m), 2 * max(nw, n));
  if (! vmcomplete(vmblock))                     /* lack of memory?   */
    return 7;                                    /* nothing goes!     */

  /* ---------------------- initialize graphics --------------------- */

  grafiktreiber = DETECT;        /* find the optimal graphics driver  */
  initgraph(&grafiktreiber,      /* in initgraph()                    */
            &grafikmodus,
            getenv("BGI"));
  grafikfehler = graphresult();  /* get BGI error status              */
  if (grafikfehler != grOk)      /* graphics initializiation failed?  */
  {
    vmfree(vmblock);             /* free dynamic memory               */
    return 6;                    /* report error                      */
  }

  texthoehe = textheight("X");             /* height of a text line   */

  fensterx = 8;                            /* define plot window      */
  fensterb = getmaxx() - 8 - fensterx;
  fenstery = 2 * texthoehe + 14;
  fensterh = getmaxy() - 8 - fenstery;

  /* ------------ find enclosing hexahedron of all points ----------- */

  fminmax(c, nv, nw, d, m, n,           /* find smallest and largest  */
          &r3xmin, &r3xmax,             /* x-, y- and z-coordinate    */
          &r3ymin, &r3ymax,             /* of surface and node points */
          &r3zmin, &r3zmax
         );

  /* ----------- project all points onto the die x-z-plane ---------- */

  projizieren(c, nv, nw);                           /* surface points */
  projizieren(d, m,  n);                            /* node points    */

  /* ----------------- find minimal rectangle that   ---------------- */
  /* ----------------- contains all projected points ---------------- */

  xmin = r3xmin - r3ymin * (p.x / p.y);  /* find smallest and largest */
  zmin = r3zmin - r3ymin * (p.z / p.y);  /* x- and z-coordinate       */
  xmax = r3xmax - r3ymax * (p.x / p.y);  /* of projected points       */
  zmax = r3zmax - r3ymax * (p.z / p.y);

  if (xmin == xmax)                      /* width of rectangle zero?  */
    xmin = c[0][0][0] - ONE,             /* set to two for            */
    xmax = xmin       + TWO;             /* transformation            */
  if (zmin == zmax)                      /* height of rectangle zero? */
    zmin = c[0][0][2] - ONE,             /* set to two for            */
    zmax = zmin       + TWO;             /* transformation            */

  /* ------------------------- print heading ------------------------ */

  itoa(m * n, mn_string, 10);
  mein_setcolor(YELLOW);
  outtextxy(0, 0, "a surface in R3 given by");
  moveto(0, texthoehe + 3);
  outtext(mn_string);
  outtext(" nodes marked by crosses");

  /* ------------------- set plot window and draw ------------------- */

  setviewport(fensterx, fenstery,  /* shift screen origin from (0,0) */
              fensterx + fensterb, /* to (fensterx,fenstery)         */
              fenstery + fensterh, /* (no clipping)                  */
              0);

  mein_setcolor(LIGHTCYAN);             /* draw a frame around window */
  moveto(-8,           -8);             /* in light cyan (8 points    */
  lineto(fensterb + 8, -8);             /* from proper window)        */
  lineto(fensterb + 8, fensterh + 8);
  lineto(-8,           fensterh + 8);
  lineto(-8,           -8);

  /* ------------------------ draw hexahedron ----------------------- */

  faktor = (r3ymin - r3ymax) / p.y;
  diffx  = (int)((faktor * p.x) / (xmax - xmin) * fensterb);
  diffy  = (int)((faktor * p.z) / (zmax - zmin) * fensterh);

  mein_setcolor(WHITE);
  moveto(0,                fensterh);           /* visable edges      */
  lineto(fensterb - diffx, fensterh);           /* solid              */
  lineto(fensterb,         fensterh - diffy);
  lineto(fensterb,         0);
  lineto(fensterb - diffx, diffy);
  lineto(0,                diffy);
  lineto(0,                fensterh);
  moveto(fensterb - diffx, fensterh);
  lineto(fensterb - diffx, diffy);
  moveto(fensterb,         0);
  lineto(diffx,            0);
  lineto(0,                diffy);

  setlinestyle(DASHED_LINE, 0, NORM_WIDTH);     /* invisible edges    */
  moveto(0,                fensterh);           /* dashed             */
  lineto(diffx,            fensterh - diffy);
  lineto(fensterb,         fensterh - diffy);
  moveto(diffx,            fensterh - diffy);
  lineto(diffx,            0);
  setlinestyle(SOLID_LINE, 0, NORM_WIDTH);      /* remaining lines    */
                                                /* solid              */

  /* ----------------- draw surface and node points ----------------- */

  ftrafo(c, nv, nw,                           /* transform surface    */
         xmin, xmax,                          /* points into plot     */
         zmin, zmax,                          /* window and draw      */
         fensterb, fensterh,                  /* curve polygons in    */
         polypoints                           /* white (both in       */
        );                                    /* w direction and in   */
  mein_setcolor(WHITE);                       /* v direction so that  */
  for (i = 0; i < nv; i++)                    /* they form a grid)    */
    drawpoly(nw, polypoints[i]);
  for (j = 0; j < nw; j++)
  {
    moveto(polypoints[0][2 * j],
           polypoints[0][2 * j + 1]);
    for (i = 1; i < nv; i++)
      lineto(polypoints[i][2 * j],
             polypoints[i][2 * j + 1]);
  }

  ftrafo(d, m, n,                             /* transform node       */
         xmin, xmax,                          /* points into plot     */
         zmin, zmax,                          /* window...            */
         fensterb, fensterh,
         polypoints
        );
  mal_stuetz(polypoints, m, n,                /* ... and draw         */
             voffen, woffen, st_gitter);

  /* -------------------------- finish plot ------------------------- */

  graphdefaults();                        /* reverse coordinate shift */

  mein_setcolor(LIGHTMAGENTA);            /* print input demand in    */
  outtextxy(getmaxx() - 200, 0,           /* light magenta            */
            "continue with any key...");
  bioskey(0);                             /* wait for key stroke      */
  closegraph();                           /* finish work with the     */
                                          /* graphics driver          */

  vmfree(vmblock);                        /* free memory of matrox    */
                                          /* of polygons              */

  grafikfehler = graphresult();
  if (grafikfehler != grOk)               /* graphics error?          */
    return 6;                             /* report it                */


  return 0;                               /* zero: no error           */
}



#elif defined(MCGRAFMOEGLICH)
/*--------------------------------------------------------------------*/

int zeigkrv2         /* plot a function table in R2 ..................*/
        (
         int  nk,                  /* size of function table .........*/
         REAL *kurve[],            /* function table .................*/
         int  ns,                  /* number of nodes ................*/
         REAL *stuetz[]            /* nodes ..........................*/
        )                          /* error code .....................*/

/***********************************************************************
* The points (kurve[i][0],kurve[i][1]), i=0,...,nk-1, form a table of  *
* values of a planar curve. They define a polygonal line which shall   *
* be plotted so that the screen is maximally filled.                   *
* To achieve this we find the extrema of both sets of coordinates      *
* (xmin, xmax, ymin, ymax) and then transform the rectangle            *
* [xmin,xmax] x [ymin,ymax]  into the plot window. There the points    *
* are represented in integer point scale coordinates for drawing by    *
* drawpoly().                                                          *
*                                                                      *
* The [0..ns-1] vector stuetz contains the nodes used to generate the  *
* table of values for the curve such as the nodes of a cubic spline or *
* the de Boor points for B-splines etc. For color monitors these will  *
* will be plotted by a green cross. Clearly the extrema search above   *
* is also extended over this node data.                                *
*                                                                      *
* This program is designed for the graphics library of QuickC 2.0.     *
* For other computers or C compilers one should use analogue functions *
* for _setvideomode(), _getvideoconfig() etc.                          *
*                                                                      *
* For PCs with a Hercules graphics card, this plot will only succeed   *
* if - prior to starting the program - the graphics driver msherc.com  *
* has been installed.                                                  *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* 0: all ok                                                            *
* 1: nk < 1                                                            *
* 2: ns < 1                                                            *
* 3: lack of memory                                                    *
* 4: Graphics error                                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, videoconfig, max, _VRES16COLOR, _setvideomode,                 *
* _TEXTMONO, _getvideoconfig, getminmax, ONE, TWO, itoa,               *
* _setcolor, YELLOW, _outtext, _settextposition, _setvieworg,          *
* LIGHTCYAN, _moveto, _lineto, transformieren, WHITE, drawpoly,        *
* stuetzpunkte_eintragen, LIGHTMAGENTA, _bios_keybrd, _KEYBRD_READ,    *
* _DEFAULTMODE, vminit, vmalloc, VVEKTOR, vmcomplete, vmfree           *
***********************************************************************/

{
  REAL  xmin, xmax,          /* extreme x- and y-values for table of  */
        ymin, ymax;          /* values and nodes                      */
  int   fensterx, fenstery,  /* left upper corner of plot window      */
        fensterb, fensterh,  /* Width and height of plot window       */
        *polypoints,         /* Vector with point scales of screen    */
                             /* points                                */
        texthoehe,           /* Hight of a line of text               */
        grafikmodus;         /* modus set up in  _setvideomode() for  */
                             /* the graphics                          */
  char  ns_string[6];        /* ns as an alpha-numeric string         */
  static                     /* desription of current graphics con-   */
  struct videoconfig graf;   /* figuration                            */
  void  *vmblock;            /* list of dynamical vectors             */


  /* ---------------- check for input errors ------------------------ */

  if (nk < 1)
    return 1;

  if (ns < 1)
    return 2;


  /* ----- allocate the vector of screen coordinates dynamically --- */
  /* ----- of sufficient size to be able to accomodate all points -- */
  /* ----- on the curve as well as all nodes after transformation -- */

  vmblock    = vminit();
  polypoints = (int *)vmalloc(vmblock, VVEKTOR,
                              max(nk, ns), 2 * sizeof(*polypoints));
  if (! vmcomplete(vmblock))                     /* lack of memory?   */
    return 3;                                    /* nothing goes!     */


  /* ------------------- inotialize graphics ------------------------ */

  grafikmodus = _VRES16COLOR;             /* find best graphics modus */
  while (! _setvideomode(grafikmodus))
    grafikmodus--;
  if (grafikmodus == _TEXTMONO)                 /* no success?        */
  {
    vmfree(vmblock);
    return 4;
  }

  _getvideoconfig(&graf);
  texthoehe = graf.numypixels / graf.numtextrows + 1;

  fensterx = 8;                            /* define the plot window  */
  fensterb = graf.numxpixels - 9 - fensterx;
  fenstery = 2 * texthoehe + 14;
  fensterh = graf.numypixels - 9 - fenstery;


  /* ----------- find the minimal rectangle for all points ---------- */
  /* ----------- on curve and for nodes                    ---------- */

  getminmax(nk, kurve,        /* find extrema of x- and y-coordinates */
            ns, stuetz,       /* for nodes and curve points           */
            &xmin, &xmax,
            &ymin, &ymax
           );
  if (xmin == xmax)                        /* x-interval length 0 ?   */
    xmin = kurve[0][0] - ONE,         /* set equal to 2 for transform */
    xmax = xmin        + TWO;
  if (ymin == ymax)                        /* y-interval length 0 ?   */
    ymin = kurve[0][1] - ONE,              /* set to 2 for transform  */
    ymax = ymin        + TWO;


  /* ------------------- put out heading ---------------------------- */

  itoa(ns, ns_string, 10);
  _settextcolor(YELLOW);
  _settextposition(1, 1);
  _outtext("A planar curve marked by crosses at the");
  _settextposition(2, 1);
  _outtext(ns_string);
  _outtext(" nodes");


  /* ------------ set up plot window and draw ----------------------- */

  _setvieworg(fensterx, fenstery); /* shift screen coordinates from   */
                                   /* (0,0) to (fensterx,fenstery)    */

  _setcolor(LIGHTCYAN);             /* draw a frame in light cyan     */
  _moveto(-8,           -8);        /* 8 points outside proper window */
  _lineto(fensterb + 8, -8);
  _lineto(fensterb + 8, fensterh + 8);
  _lineto(-8,           fensterh + 8);
  _lineto(-8,           -8);


  /* ------------- plot points on curve and nodes ------------------- */

  transformieren(nk, kurve,  /* transform points on curve into plot   */
                 xmin, xmax, /* and connect with a white polygon      */
                 ymin, ymax,
                 fensterb, fensterh,
                 polypoints
                );
  _setcolor(WHITE);
  drawpoly(nk, polypoints);

  transformieren(ns, stuetz,  /* transform nodes to window coordinates*/
                 xmin, xmax,                      /* and plot         */
                 ymin, ymax,
                 fensterb, fensterh,
                 polypoints
                );
  stuetzpunkte_eintragen(ns, polypoints);


  /* ------------ finish graphics ----------------------------------- */

  _settextcolor(LIGHTMAGENTA);            /* print input request in   */
  _settextposition(1, 56);                /* light magenta and        */
  _outtext("continue with any key...");
  _bios_keybrd(_KEYBRD_READ);             /* wait for key stroke      */
  _setvideomode(_DEFAULTMODE);            /* return to previous       */
                                          /* operating mode           */


  vmfree(vmblock);                       /* free dynamic allocations  */
  return 0;
}

int zeigflaeche(REAL **c[],  /* QuickC, but no Turbo C?               */
                int  nv,     /* The function is presently implemented */
                int  nw,     /* empty.                                */
                REAL **d[],
                int  m,
                int  n,
                int  voffen,
                int  woffen,
                int  st_gitter
               )
{
  return 0;
}



#endif                               /* #elif defined(MCGRAFMOEGLICH) */
#else                        /* no graphics supported or wanted?      */

#include <basis.h>           /*  for  REAL                            */

int zeigkrv2
        (
         int  nk,
         REAL *kurve[],
         int  ns,
         REAL *stuetz[]
        )                    /* define these functions anyway as empty*/
{                            /* for GNU CC 2.0 which does not allow   */
  return 0;                  /* empty ANSI C source files, what a pain*/
}                            /* If desired and possible, the user can */
                             /* try to include his routine for his    */
int zeigflaeche              /* special compiler here instead.        */
        (
         REAL **c[],
         int  nv,
         int  nw,
         REAL **d[],
         int  m,
         int  n,
         int  voffen,
         int  woffen,
         int  st_gitter
        )
{
  return 0;
}
#endif                                            /* #ifdef MITGRAFIK */

/* ------------------------- END zeigkrv2.c ------------------------- */
