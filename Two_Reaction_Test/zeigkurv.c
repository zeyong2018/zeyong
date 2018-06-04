#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE zeigkurv.c ------------------------ */

/***********************************************************************
*                                                                      *
* Module for plotting a table of values for a planar curve with its    *
* -----------------------------------------------------------------    *
* underlying nodes                                                     *
* -----------------                                                    *
*                                                                      *
* Programming language: Turbo C 2.0                                    *
* Compiler:             Turbo C 2.0, PureC 1.0                         *
* Author:               Juergen Dietel, Computer Center, RWTH Aachen   *
* Date:                 1.22.1993 - 2.26.1997                          *
*                                                                      *
***********************************************************************/

#include <cnumgraf.h>         /*  for  MITGRAFIK,                     */
                              /*       BGIGRAFMOEGLICH, GRAFIK_H,     */
                              /*       BIOS_H, INIT_FONTS             */

#if defined(MITGRAFIK) &&     /* graphics supported or wanted?       */\
    defined(BGIGRAFMOEGLICH)

#include GRAFIK_H             /*  for  initgraph, closegraph, grOk    */
                              /*       DETECT, outtextxy, setcolor,   */
                              /*       graphresult, getmaxcolor,      */
                              /*       SOLID_LINE, THICK_WIDTH,       */
                              /*       YELLOW, WHITE, moveto,         */
                              /*       LIGHTGREEN, LIGHTMAGENTA,      */
                              /*       lineto, setlinestyle, outtext, */
                              /*       LIGHTRED, NORM_WIDTH,          */
                              /*       setviewport, graphdefaults     */
#include <stdlib.h>           /*  for  getenv, itoa                   */
#include BIOS_H               /*  for  bioskey                        */
#include <basis.h>            /*  for  REAL, ONE, TWO, min, max       */
#include <vmblock.h>          /*  for  vminit, vmalloc, IMATRIX,      */
                              /*       vmcomplete, vmfree, VVEKTOR    */
#include <zeigkurv.h>         /*  for  zeigkurv                       */



/*--------------------------------------------------------------------*/

static void getminmax(int  n,
                      REAL x[],
                      REAL *min,
                      REAL *max
                     )

/***********************************************************************
* Find minimal and maximal entry of a  [0..n-1] vector x               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  REAL mi,
       ma;

  for (mi = ma = *x, n--, x++; n != 0; n--, x++)
  {
    if (*x < mi)
      mi = *x;
    if (*x > ma)
      ma = *x;
  }
  *min = mi;
  *max = ma;
}



/*--------------------------------------------------------------------*/

static void transformieren(int  n,
                           REAL x[],
                           int  polypoints[],
                           REAL min,
                           REAL max,
                           int  grafmax
                          )

/***********************************************************************
* Transform the components of the  [0..n-1] vector x linearly from the *
* interval [min, max] into the interval [0, grafmax] so that min maps  *
* to 0 and max maps to grafmax. The result is copied to the vector     *
* polypoints. On first call, all elements of polypoints                *
* with even indices are used; on second call all odd indexed elements  *
* of polypoints are used. Moreover on second call min is mapped to     *
* grafmax and max to 0 for the y-coordinates to avoid an upside down   *
* picture.                                                             *
* A third call achieves the same as the first, a fourth the same as    *
* the second, etc.                                                     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  static int offset = 0;     /* start of indexing of polypoints       */
  int        *hilf;          /* Address of polypoint[offset + 2 * i]  */
  REAL       faktor;         /* slope of linear transformation        */

  hilf = polypoints + offset;
  for (faktor = grafmax / (max - min);
       n != 0;
       n--, hilf += 2, x++)
  {
    *hilf = (*x - min) * faktor;
    if (offset)                        /* working on  y-coordinates ? */
      *hilf = grafmax - *hilf;                /* reflect coordinates  */
  }

  offset = ! offset;           /* toggle between x- and y-coordinates */
}



/*--------------------------------------------------------------------*/

static void mein_setcolor(int farbe)

/***********************************************************************
* Set up a new color, be careful not to excced the maximal number of   *
* colors.                                                              *
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

static void stuetzpunkte_eintragen(int ns,
                                   int polypoints[]
                                  )

/***********************************************************************
* Mark the points (polypoints[2*i],polypoints[2*i+1]) of the           *
* [0..2*ns-1] vector polypoints with a cross, and the first point of   *
* polypoints additionally with a rectangle.                            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* mein_setcolor, LIGHTGREEN, setlinestyle, SOLID_LINE, THICK_WIDTH,    *
* line, LIGHTRED, NORM_WIDTH, rectangle,                               *
***********************************************************************/

{
  int i,              /* Loop variable                                */
      x0, y0;         /* screen coordinates of the current point      */


  mein_setcolor(LIGHTGREEN);                   /* crosses light green */
  setlinestyle(SOLID_LINE, 0, THICK_WIDTH);    /* and wide            */


  for (i = 0; i < ns; i++)
  {
    x0 = *polypoints++;
    y0 = *polypoints++;
    line(x0 - 4, y0 - 4, x0 + 4, y0 + 4);
    line(x0 + 4, y0 - 4, x0 - 4, y0 + 4);

    if (i == 0)                                /* first node ?        */
    {                                          /* add a red rectangle */
      mein_setcolor(LIGHTRED);
      setlinestyle(SOLID_LINE, 0, NORM_WIDTH);
      rectangle(x0 - 6, y0 - 6, x0 + 6, y0 + 6);
      mein_setcolor(LIGHTGREEN);
      setlinestyle(SOLID_LINE, 0, THICK_WIDTH);
    }
  }
}



/*--------------------------------------------------------------------*/

int zeigkurv(int  n,
             int  ns,
             REAL x[],
             REAL y[],
             REAL stuetz_x[],
             REAL stuetz_y[]
            )

/***********************************************************************
* The points (x[i],y[i]), i=0,1,...,n-1 form the table of values.      *
* The polygon defined by these points should be drawn on one plot that *
* optimally covers the range of occuring x- and y-values. For this     *
* purpose we find the extreme x- and y-values of the above table of    *
* values and of the nodes of the curve (xmin, xmax, ymin, ymax). Then  *
* the x[i] and y[i] must be transformed linearly from their ranges     *
* [xmin, xmax] and [ymin, ymax] to the coordinates of the plot window. *
* Finally the transformed plot coordinate indices are stored in the    *
* integer vector int polypoints[] in order to be able to draw their    *
* polygon with drawpoly().                                             *
* The nodes are given by two [0..ns-1] vectors stuetz_x and stuetz_y . *
* On color monitors the nodes are shown as green crosses in addition   *
* to the polygon from the curve points in the table of values.         *
* This program is designed for the BGI of Turbo C 2.0 for IBM PCs.     *
* For other computers or C compilers one should use analogue functions *
* for initgraph(), closegraph() etc.                                   *
* The graphics drivers that are loaded via initgraph() are expected to *
* be available in the directory the name of which - prior to the start *
* of the program - was assigned als value to the environment variable  *
* BGI, e.g. use the MS-DOS command "set BGI=c:\tc", if the directory   *
* C:\TC contains the needed graphics drivers *.BGI .                   *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: all is ok                                                         *
* 1: n < 1                                                             *
* 2: ns < 1                                                            *
* 3: lack of available memory                                          *
* 4: BGI graphics error                                                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, DETECT, initgraph, getenv, getminmax, min, max,                *
* ONE, TWO, textheight, getmaxx, getmaxy, itoa, mein_setcolor, YELLOW, *
* outtextxy, moveto, outtext, setviewport, LIGHTCYAN, lineto,          *
* transformieren, WHITE, drawpoly, stuetzpunkte_eintragen,             *
* graphdefaults, LIGHTMAGENTA, bioskey, closegraph, graphresult, grOK, *
* vminit, vmalloc, VVEKTOR, vmcomplete, vmfree, INIT_FONTS             *
***********************************************************************/

{
  REAL  xmin, xmax,                /* min and max of x- and y-values  */
        ymin, ymax,                /* in table                        */
        smin, smax;                /* Min and max values of nodes     */
  int   fensterx, fenstery,        /* left upper corner of window     */
        fensterb, fensterh,        /* width and height of window      */
                                   /* for the plot                    */
        *polypoints,               /* vector of screen coordinates    */
        grafiktreiber = DETECT,    /* find proper graphics driver     */
                                   /* by calling initgraph()          */
        grafikmodus,               /* modus set by initgraph() used   */
                                   /* by the graphics driver          */
        grafikfehler,              /* error code of the BGI           */
        texthoehe;                 /* font size of X in points        */
  char  ns_string[6];              /* ns as alpha-numeric string      */
  void  *vmblock;                  /* Liste der dynamisch             */
                                   /* vereinbarten Vektoren           */


  /* -------------- check input data -------------------------------- */

  if (n < 1)
    return 1;

  if (ns < 1)
    return 2;


  /* ----- allocate the vector of screen coordinates dynamically --- */
  /* ----- of sufficient size to be able to accomodate all points -- */
  /* ----- on the curve as well as all nodes after transformation -- */

  vmblock    = vminit();
  polypoints = (int *)vmalloc(vmblock, VVEKTOR, max(n, ns),
                              2 * sizeof(*polypoints));
  if (! vmcomplete(vmblock))                     /* lack of memory?   */
    return 3;                                    /* nothing goes!     */


  /* ------------------- initialize graphics  ----------------------- */

  initgraph(&grafiktreiber, &grafikmodus, getenv("BGI"));
  if ((grafikfehler = graphresult()) != grOk)         /* no success ? */
  {
    vmfree(vmblock);
    return 4;
  }

  INIT_FONTS;


  getminmax(n, x, &xmin, &xmax);          /* find maximal and minimal */
  getminmax(ns, stuetz_x, &smin, &smax);  /* x-coordinate of points   */
  xmin = min(xmin, smin);                 /* and nodes on curve       */
  xmax = max(xmax, smax);
  if (xmin == xmax)
    xmin = x[0] - ONE,
    xmax = xmin + TWO;
  getminmax(n, y, &ymin, &ymax);          /* find maximal and minimal */
  getminmax(ns, stuetz_y, &smin, &smax);  /* y-coordinates of points  */
  ymin = min(ymin, smin);                 /* on the curve and nodes   */
  ymax = max(ymax, smax);
  if (ymin == ymax)
    ymin = y[0] - ONE,
    ymax = ymin + TWO;

  texthoehe = textheight("X");

  fensterx = 8;                          /* define window for drawing */
  fensterb = getmaxx() - 8 - fensterx;
  fenstery = 2 * texthoehe + 14;
  fensterh = getmaxy() - 8 - fenstery;


  itoa(ns, ns_string, 10);      /* put out diagram headings in yellow */
  mein_setcolor(YELLOW);
  moveto(0, 0);
  outtext("A planar curve based on the ");
  outtext(ns_string);
  moveto(0, texthoehe + 3);
  outtext("nodes marked by crosses");

  setviewport(fensterx, fenstery,  /* set up display window, i.e.,    */
              fensterx + fensterb, /* find the screen coordinates of  */
              fenstery + fensterh, /* (0,0) and shift them to         */
              0);                  /* (fensterx,fenstery)             */


  mein_setcolor(LIGHTCYAN);              /* draw a frame around the   */
  moveto(-8,           -8);              /* whole window in light cyan*/
  lineto(fensterb + 8, -8);              /* (8 points outside the     */
  lineto(fensterb + 8, fensterh + 8);    /* proper dimensions of the  */
  lineto(-8,           fensterh + 8);    /* window                    */
  lineto(-8,           -8);

  transformieren(n, x, polypoints, xmin, xmax,    /* map curve points */
                 fensterb);                       /* into the drawing */
  transformieren(n, y, polypoints, ymin, ymax,    /* window           */
                 fensterh);

  mein_setcolor(WHITE);     /* draw polygon through the points of the */
  drawpoly(n, polypoints);  /* curve in white                         */


  transformieren(ns, stuetz_x, polypoints, xmin,  /* transform nodes  */
                 xmax, fensterb);                 /* to the window of */
  transformieren(ns, stuetz_y, polypoints, ymin,  /* drawing          */
                 ymax, fensterh);                 /* transform ...    */

  stuetzpunkte_eintragen(ns, polypoints);         /* ... and draw     */

  graphdefaults();              /* reverse shift of coordinate system */


  mein_setcolor(LIGHTMAGENTA);            /* write input demand in    */
  outtextxy(getmaxx() - 200, 0,           /* light magenta            */
            "continue with any key...");
  bioskey(0);                             /* wait for key stroke      */


  closegraph();                    /* stop work with graphics driver  */

  if ((grafikfehler = graphresult()) != grOk)     /* Graphics error ? */
  {
    vmfree(vmblock);
    return 4;
  }


  grafikfehler = grafikfehler;  /* stroke the hand of compiler, haha  */
  vmfree(vmblock);              /* free dynamic allocations           */
  return 0;
}



#else                            /* no graphics supported or wanted?  */

#include <basis.h>               /*  for  REAL                        */
int zeigkurv(int  n,             /* define this function as empty     */
             int  ns,            /* to help GNU CC 2.0 which does not */
             REAL *x,            /* allow empty source files in ANSI C*/
             REAL *y,
             REAL *stuetz_x,     /* If desired and possible the above */
             REAL *stuetz_y      /* function could be rewritten for   */
            )                    /* another C compiler.               */
{
  return 0;
}
#endif                                            /* #ifdef MITGRAFIK */

/* ------------------------- END zeigkurv.c ------------------------- */
