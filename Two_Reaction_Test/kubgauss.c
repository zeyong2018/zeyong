#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE kubgauss.c ----------------------- */

#include <basis.h>         /*  for  REAL, ZERO, TWO, ONE, HALF, THREE */
#include <kubatur.h>       /*  for  Kub4GauE, Kub4GauV, Kub3GauN      */
/*.BA*/



/*.FE{C 16.6}{Gau"s Cubature Formulas for Rectangles}
             {Gau"s Cubature Formulas for Rectangles}*/

/*.BE*/
/***********************************************************************
*       Here the constants for cubature over rectangles using          *
*       Newton-Cotes formulas are provided.                            *
*                                                                      *
*       The structure kub4_ta contains a pair of doubles:              *
*         a node t and the corresponding weight a.                     *
*                                                                      *
*       The arrays K_i (i = 0..7) contain (i+1) such pairs             *
*         according to the order of the cubature formula.              *
*                                                                      *
*       The components of the array KubArr point to the K_i so that    *
*         by `KubArr[n][j]' the jth pair of a cubature formula of      *
*         order n can be addressed.                                    *
*                                                                      *
*       Author:        Uli Eggermann, 02.17.1991                       *
***********************************************************************/

typedef struct
{
  REAL t, a;                                   /* node t und weight a */
} kub4_ta;
/*.IX{kub4\unt ta}*/

static kub4_ta K_0 [] =  {{              ZERO,  TWO              }};
static kub4_ta K_1 [] =  {{ -.577350269189626,  ONE              },
                          {  .577350269189626,  ONE              }};
static kub4_ta K_2 [] =  {{ -.774596669241483, .5555555555555556 },
                          {              ZERO, .8888888888888888 },
                          {  .774596669241483, .5555555555555556 }};
static kub4_ta K_3 [] =  {{ -.861136311594053, .347854845137454  },
                          { -.339981043584856, .652145154862546  },
                          {  .339981043584856, .652145154862546  },
                          {  .861136311594053, .347854845137454  }};
static kub4_ta K_4 [] =  {{ -.906179845938664, .236926885056189  },
                          { -.538469310105683, .478628670499366  },
                          {              ZERO, .5688888888888889 },
                          {  .538469310105683, .478628670499366  },
                          {  .906179845938664, .236926885056189  }};
static kub4_ta K_5 [] =  {{ -.9324695142031521,.17132449237917   },
                          { -.661209386466265, .360761573048139  },
                          { -.238619186083197, .467913934572691  },
                          {  .238619186083197, .467913934572691  },
                          {  .661209386466265, .360761573048139  },
                          {  .9324695142031521,.17132449237917   }};
static kub4_ta K_6 [] =  {{ -.949107912342759, .12948496616887   },
                          { -.741531185599394, .279705391489277  },
                          { -.405845151377397, .381830050505119  },
                          {              ZERO, .417959183673469  },
                          {  .405845151377397, .381830050505119  },
                          {  .741531185599394, .279705391489277  },
                          {  .949107912342759, .12948496616887   }};
static kub4_ta K_7 [] =  {{ -.960289856497536, .101228536290376  },
                          { -.7966664774136269,.222381034453374  },
                          { -.525532409916329, .313706645877887  },
                          { -.18343464249565,  .362683783378362  },
                          {  .18343464249565,  .362683783378362  },
                          {  .525532409916329, .313706645877887  },
                          {  .7966664774136269,.222381034453374  },
                          {  .960289856497536, .101228536290376  }};

static kub4_ta *KubArr [] = { K_0, K_1, K_2, K_3, K_4, K_5, K_6, K_7 };
/*.IX{KubArr}*/

#define kub4_max  ((sizeof (KubArr) / sizeof (KubArr[0])) - 1)

/* ------------------------------------------------------------------ */
/*.BA*/

int Kub4GauE (
/*.IX{Kub4GauE}*/
              REAL  a, REAL b, int Nx,
              REAL  c, REAL d, int Ny,
              int   Verf,
              REAL  f (REAL,REAL),
              REAL* Wert,
              int   Schaetzen,
              REAL* FehlerSch
             )
/***********************************************************************
* Cubature over rectangles using Newton-Cotes formulas.                *
*                                                                      *
* Integrate the function  f(x,y) using the summed cubature formula     *
* of Gauss for the rectangle (a,b; c,d).                               *
* The edges of the sub-rectangles have the lengths:  (b-a) / Nx   and  *
* (d-c) / Ny .                                                         *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*   REAL a, b, Nx    left, right x-end points, number of intervals     *
*   REAL c, d, Ny    left, right y-end points, number of intervals     *
*   int  Verf        order of method  (0 <= Verf <= 7)                 *
*   REAL f ()        function                                          *
*   REAL *Wert       value of integral                                 *
*   int  Schaetzen   if nonzero : find estimate                        *
*   REAL *FehlerSch  error estimate, obtained if desired by a second   *
*                    cubature pass using half of the step size         *
*                                                                      *
* Return value :                                                       *
*   0:               o.k.                                              *
*   1:               Nx improper                                       *
*   2:               Ny improper                                       *
*   3:               order incorrect                                   *
*   4:               Integration interval has length zero              *
*                                                                      *
* Author             Uli Eggermann, 03.31.1996                         *
.BA*)
***********************************************************************/
/*.BE*/
{
  int  i, j, k, u, v, Nab, Ncd;
  REAL Hab, Hcd, Wert1 = ZERO;

  if (Nx < 1)               return (1);
  if (Ny < 1)               return (2);
  if (Verf < 0 || Verf > 7) return (3);
  if (a == b || c == d)     return (4);

  for (k = 1; k <= (Schaetzen ? 2 : 1); k++)
  {
    Nab = k * Nx;                            /* number of intervals   */
    Ncd = k * Ny;
    Hab = HALF * (b - a) / Nab;              /* half of x-step size   */
    Hcd = HALF * (d - c) / Ncd;              /* half of y-step size   */

    #define  kub  KubArr [Verf]

    *Wert = ZERO;                            /* initialize            */

    for (i = 0; i < Nab; i++)
      for (j = 0; j < Ncd; j++)
        for (u = 0; u <= Verf; u++)
          for (v = 0; v <= Verf; v++)
          {
            REAL w = Hab * Hcd * kub[u].a * kub[v].a;
            REAL x = a + Hab * (kub[u].t + 2 * i + 1);
            REAL y = c + Hcd * (kub[v].t + 2 * j + 1);
            REAL z = f (x, y);
            *Wert += w * z;
          }
    if (Schaetzen && k == 1) Wert1 = *Wert;       /* store value      */
  }

  if (Schaetzen) {                                     /* estimate    */
    *FehlerSch = (*Wert - Wert1) / THREE;
  }

  return (0);
  #undef kub
}


/* ------------------------------------------------------------------ */
/*.BA*/

int Kub4GauV (
/*.IX{Kub4GauV}*/
              REAL* x, int Nx,
              REAL* y, int Ny,
              int   Verf,
              REAL  f (REAL,REAL),
              REAL* Wert,
              int   Schaetzen,
              REAL* FehlerSch
             )
/***********************************************************************
* Cubature over rectangles using  Newton-Cotes formulas.               *
*                                                                      *
* Integrate the function f(x,y) using the summed Gaussian cubature     *
* formula on the rectangle (a,b) x (c,d).                              *
*                                                                      *
* The edge lengths for the sub-rectangles are not identical, however,  *
* as in Kub4GauE, but are given by thr two vectors  X []  and  Y [] .  *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   REAL X []        vector of x-interval end points:                  *
*                      a = X[0] < X[1] < .. < X[Nx] = b                *
*   REAL Y []        vector of y-interval end points:                  *
*                      c = Y[0] < Y[1] < .. < Y[Ny] = d                *
*   int  Verf        order of method  (0 <= Verf <= 7)                 *
*   REAL f ()        Function                                          *
*   int  Schaetzen   0: no error estimation                            *
*                    1: error estimation using a second pass of        *
*                       cubature for half the step size                *
*   REAL *Wert       value for integral                                *
*   REAL *FehlerSch  error estimate of  Wert                           *
*                                                                      *
* Return value :                                                       *
* --------------                                                       *
*   0:               o.k.                                              *
*   1:               order number improper                             *
*   2:               x-interval of length zero                         *
*   3:               y-interval of length zero                         *
*                                                                      *
* Author:         Uli Eggermann, 03.31.1991                            *
.BA*)
***********************************************************************/
/*.BE*/
{
  int    k, Si, Sj, i, j, u, v;
  REAL   Hx, Hy, Wert1 = ZERO, Tx, Ty;

  if (Verf < 0 || Verf > 7) return 1;             /* check valid order*/
  for (i = 0; i < Nx; i++)                        /* X-interval test  */
    if (x [i+1] <= x [i]) return 2;
  for (j = 0; j < Ny; j++)                        /* Y-interval test  */
    if (y [j+1] <= y [j]) return 3;

  #define KUB KubArr[Verf]

  for (k = 1; k <= (Schaetzen ? 2 : 1); k++)      /* with estimate ?  */
  {                                               /* initialize       */
    for (*Wert = ZERO, i = 0; i < Nx; i++)        /* X-intervals      */
      for (Hx = (x [i+1] - x [i]) / (2*k),        /* halve X-direction*/
           Si = 1; Si < 2*k; Si += 2)
        for (Tx = x [i] + Si * Hx,                /* X-interval center*/
             j = 0; j < Ny; j++)                  /* Y-intervals      */
          for (Hy = (y [j+1] - y [j]) / (2*k),    /* halve Y-direction*/
               Sj = 1; Sj < 2*k; Sj += 2)
            for (Ty = y [j] + Sj * Hy,            /* Y-interval center*/
                 u = 0; u <= Verf; u++)
              for (v = 0; v <= Verf; v++)

                *Wert +=         Hx * KUB[u].a  *    Hy * KUB[v].a  *
                         f (Tx + Hx * KUB[u].t, Ty + Hy * KUB[v].t);

    if (Schaetzen && k == 1) Wert1 = *Wert;       /* store value      */
  } /* for (k) */
  #undef KUB

  if (Schaetzen) {                                /* estimate         */
    *FehlerSch = (*Wert - Wert1) / THREE;
  }

  return 0;

}
/*.BA*/



/*.FE{C 16.7}{Gau"s Cubature Formulas for Triangles}
             {Gau"s Cubature Formulas for Triangles}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int Kub3GauN (
/*.IX{Kub3GauN}*/
              REAL Px, REAL Py,
              REAL Qx, REAL Qy,
              REAL Rx, REAL Ry,
              int n,   int m,
              REAL f(REAL,REAL),
              REAL* Wert
             );
/***********************************************************************
* Cubature over triangular regions using the n-point Gauss formula     *
*                                                                      *
* Integrate the function f (x,y) over the triangle PQR using the summed*
* n-point Gauss cubature formula on m x m subtriangles.                *
* The subtriangles have edges of length 1/m of the original edge       *
* lengths.                                                             *
.BE*)
*                                                                      *
* Input parameters :                                                   *
*   REAL   Px,Py   coordinates of  P                                   *
*   REAL   Qx,Ry   coordinates of  Q                                   *
*   REAL   Rx,Ry   coordinates of  R                                   *
*   int    n       order of method (= number of points in each sub-    *
*                  triangle)                                           *
*   int    m       number of subtriangles along one edge               *
*   REAL   f ()    function                                            *
*   REAL   *Wert   value of integral                                   *
*                                                                      *
* Return value :                                                       *
*   0:             o.k.                                                *
*   1:             m  improper                                         *
*   2:             the corners P, Q and R are collinear                *
*   3:             nth order method not implemented                    *
*                                                                      *
* Author:          Uli Eggermann, 8.1.1990                             *
.BA*)
***********************************************************************/
/*.BE*/

/*   Epsilon serves as a check of collinearity; if the area of P,Q and*/
/*   R has area less than Epsilon/2, we judge the three points to be  */
/*   collinear.                                                       */

#define GauNEpsilon (REAL)0.0001

/***********************************************************************
* Constants for n-point Gaussian cubature:                             *
***********************************************************************/

typedef struct
{
  REAL w, x, y;                         /* weights, x-, y-coordinates */
} Tripel;
/*.IX{Tripel}*/

Tripel Gau1Konst [] = { { 1./2., 1./3., 1./3. }  };
Tripel Gau2Konst [] = { { 1./4., 1./6., 1./2. },
                        { 1./4., 1./2., 1./6. }  };
Tripel Gau3Konst [] = { { 1./6., 1./6., 1./6. },
                        { 1./6., 2./3., 1./6. },
                        { 1./6., 1./6., 2./3. }  };
Tripel Gau7Konst [] =
{
  { 0.1125,              0.3333333333333333,  0.3333333333333333 },
  { 0.0661970763942531,  0.4701420641051151,  0.4701420641051151 },
  { 0.0661970763942531,  0.05971587178976981, 0.4701420641051151 },
  { 0.0661970763942531,  0.4701420641051151,  0.05971587178976981 },
  { 0.06296959027241357, 0.1012865073234563,  0.1012865073234563 },
  { 0.06296959027241357, 0.7974269853530873,  0.1012865073234563 },
  { 0.06296959027241357, 0.1012865073234563,  0.7974269853530873 }
};

static Tripel *GauN [] =
{
  Gau1Konst, Gau2Konst, Gau3Konst, NULL, NULL, NULL, Gau7Konst
};

#define GauNmax  sizeof (GauN) / sizeof (GauN[0])

/***********************************************************************
* Summed cubature over triangle using  Newton-Cotes                    *
***********************************************************************/

int Kub3GauN (REAL Px, REAL Py, REAL Qx, REAL Qy,
              REAL Rx, REAL Ry, int n, int m,
              REAL f(REAL,REAL),  REAL* Wert)
{
  int    d, i, j, k;
  REAL   Fak, Area, Dx, Dy;
  REAL   hPQx, hPQy, hPRx, hPRy;

  #define G  GauN[n-1][k]                            /* Abbreviations */
  #define X  Dx + Fak * (G.x * hPQx + G.y * hPRx)
  #define Y  Dy + Fak * (G.x * hPQy + G.y * hPRy)

  if (m < 1)                   return 1;      /* m  o.k. ?            */

  Area = Px * (Qy - Ry)                       /* Test collinearity    */
       + Qx * (Ry - Py)
       + Rx * (Py - Qy);
  if (FABS(Area) < GauNEpsilon) return 2;

  if (n < 1 || n > GauNmax ||                 /* desired method       */
      GauN [n-1] == NULL)      return 3;      /*  implemented ?       */

  *Wert = ZERO;                               /* initialize           */
  Area /= (m * m) ;                           /* double triangle area */
  hPQx = (Qx - Px) / m;
  hPRx = (Rx - Px) / m;                       /* edge vectors for the */
  hPQy = (Qy - Py) / m;                       /*   m * m              */
  hPRy = (Ry - Py) / m;                       /*   subtriangles       */

  for (d = 0; d < 2; d++)                     /* types of triangles d */
  {
    Fak = d ? -ONE : ONE;                     /* d = 1: reflected     */
    for (j = d; j < m; j++)                   /* j: along   PR        */
      for (i = d; i < m - j + d; i++)         /* i: along   PQ        */
      {
        Dx = Px + i * hPQx + j * hPRx;        /* (Dx,Dy) ist top      */
        Dy = Py + i * hPQy + j * hPRy;        /* corner of subtriangle*/
        for (k = 0; k < n; k++)               /* Sum of weighted      */
         *Wert += G.w * f (X, Y);             /*   function values    */
      }                                       /* (see above for X, Y) */
  }

  *Wert *= Area;        /* multiply by double the subarea             */
  return 0;

  #undef GauNEpsilon
  #undef GauNmax
  #undef G
  #undef X
  #undef Y              /* cancel abbreviations                       */
}

/* ------------------------- END kubgauss.c ------------------------- */
