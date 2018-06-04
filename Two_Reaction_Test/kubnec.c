#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE kubnec.c ------------------------ */

#include <basis.h>      /*  for  REAL, ZERO, FABS, ONE, TWO, POW      */
#include <vmblock.h>
#include <kubatur.h>    /*  for  Kub4NeCn, Kub3NeC3, Kub4RoRi,        */
                        /*       Kub4BuRi, Kub3RoRi                   */

static REAL K4KnotGew (int, int);
static void Kub4RoST  (REAL a, REAL b, int Nx, REAL c, REAL d,
                       int Ny, int nST, REAL f (REAL,REAL), REAL* W);
static int  RoRiExtr  (REAL* RoFo, int nRoFo, int Ordnung,
                       REAL* Wert, REAL* FehlerSch);
static void Kub4BuST  (REAL a, REAL b, int Nx, REAL c, REAL d,
                       int Ny, int j, REAL f (REAL,REAL), REAL* L);
static int  BuRiExtr  (REAL* BuFo, int nBuFo, int Ordnung,
                       REAL* Wert, REAL* FehlerSch);
static int  BuNenner  (int j);
static  int Kub3Nec3n (REAL Px, REAL Py, REAL Qx, REAL Qy, REAL Rx,
                       REAL Ry, int n, REAL f(REAL,REAL), REAL* W);



/*.BA*/
/*.FE{C 16.3}
     {Newton-Cotes Cubature Formulas for Rectangular Regions}
     {Newton-Cotes Cubature Formulas for Rectangular Regions}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int Kub4NeCn (
/*.IX{Kub4NeCn}*/
              REAL  a, REAL b, int Nx,
              REAL  c, REAL d, int Ny,
              int   Verfahren,
              REAL  f (REAL,REAL),
              REAL* Wert,
              int   Schaetzen,
              REAL* FehlerSch
             );
/***********************************************************************
* Cubature over rectangles using Newton-Cotes formulas.                *
*                                                                      *
* Integrate the function f(x,y) over a rectangle (a,b) x (c,d) using   *
* the summed Newton-Cotes cubature formulas for sub-rectangles.        *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   REAL   a,b        left, right x-end points                         *
*   int    Nx         number of x-intervals                            *
*   REAL   c,d        ditto for y                                      *
*   int    Ny                                                          *
*   int    Verfahren  Number of method :                               *
*                       1: trapezoidal rule                            *
*                       2: Simpson's rule                              *
*                       3: 3/8                                         *
*                       4: 4/90                                        *
*                       5: 5/288                                       *
*                       6: 6/840                                       *
*                       7: 7/17280 formula                             *
*   double  f (x,y)   function                                         *
*   REAL   *Wert      value for integral                               *
*   int    Mit_Sch    0: no estimation                                 *
*                     1: estimate error by repeating cubature for half *
*                        the step size                                 *
*   REAL   *Sch       error estimate for Wert                          *
*                                                                      *
* Return value :                                                       *
*   0:              o.k.                                               *
*   1:              Nx improper                                        *
*   2:              Ny improper                                        *
*   3:              method number incorrect                            *
*   4:              Integration interval of length zero                *
*                                                                      *
* Author:           Uli Eggermann, 3.31.1996                           *
.BA*)
***********************************************************************/
/*.BE*/

/**********************************************************************
* Global  constants and variables                                     *
**********************************************************************/
static REAL const K_1 [] = { 1./2, 1./2 };
static REAL const K_2 [] = { 1./6, 4./6, 1./6 };
static REAL const K_3 [] = { 1./8, 3./8, 3./8, 1./8 };
static REAL const K_4 [] = { 7./90, 32./90, 12./90, 32./90, 7./90 };
static REAL const K_5 [] = { 19./288, 75./288, 50./288,
                             50./288, 75./288, 19./288 };
static REAL const K_6 [] = { 41./840, 216./840, 27./840, 272./840,
                             27./840, 216./840, 41./840 };
static REAL const K_7 [] = { 751./17280, 3577./17280, 1323./17280,
                             2989./17280, 2989./17280, 1323./17280,
                             3577./17280, 751./17280 };
static REAL const *KubArr [] = { K_1, K_2, K_3, K_4, K_5, K_6, K_7 };

int KubVer = -1;            /* global variables:  method number       */
int KubX, KubY;             /* number of x-, y-intervals              */

/**********************************************************************
* Cubature over rectangles via  Newton-Cotes                          *
**********************************************************************/

int Kub4NeCn (REAL a, REAL b, int Nx, REAL c, REAL d, int Ny,
              int  Verfahren, REAL f (REAL,REAL),
              REAL* Wert, int Schaetzen, REAL* FehlerSch)
{
  int  i,   j,   k,  Ordnung;
  REAL Hab, Hcd, Wert1 = ZERO;

  if (Nx < 1)                         return (1);
  if (Ny < 1)                         return (2);
  if (Verfahren < 1 || Verfahren > 7) return (3);
  if (a == b || c == d)               return (4);

                                      /* method:     1,2,3,4,5,6,7    */
  Ordnung = (Verfahren / 2) * 2 + 2;  /* order:      2,4,4,6,6,8,8    */

  KubVer = Verfahren;                 /* copy to global variable      */

  for (k = 1; k <= (Schaetzen ? 2 : 1); k++)
  {
    KubX  = k * Nx * KubVer;
    KubY  = k * Ny * KubVer;
    *Wert = ZERO;                                    /* initialize    */
    Hab   = (b - a) / KubX;                          /* step sizes    */
    Hcd   = (d - c) / KubY;

    for (i = 0; i <= KubX; i++)                      /* Cubature      */
      for (j = 0; j <= KubY; j++) {
        *Wert += K4KnotGew (i,j) * f (a + i * Hab, c + j * Hcd);
      }
    *Wert *= Hab * Hcd * KubVer * KubVer;            /* final multipl.*/
    if (Schaetzen && k == 1) Wert1 = *Wert;          /* store value   */
  }

  if (Schaetzen)                                       /* estimate    */
  {
    *FehlerSch = (Wert1 - *Wert) / ((1 << Ordnung) - 1);
  }

  return (0);
}

static REAL K4KnotGew (int i, int j)
/*.IX{K4KnotGew}*/
/***********************************************************************
* Local function for finding node weights                              *
*                                                                      *
* Weight the functional values at the nodes depending on their location*
* (edge, interior, ...) in the summed cubature formulas.               *
***********************************************************************/
{
  REAL f;
  int  k;

  #define Faktor(a,b) ((k == 0 && a > 0 && a < b) ? TWO : ONE)

  /***********************************************************
  * for x-direction :                                        *
  *   1) node at interval end  (k == 0)                      *
  *   2) node not at left end (a > 0) nor at right end (a<b) *
  * for y-direction :                                        *
  *   1) node at interval end  (k == 0)                      *
  *   2) node not at bottom (a > 0) and not at top (a < b)   *
  ***********************************************************/

  k  = i % KubVer;
  f  = KubArr [KubVer-1] [k] * Faktor (i, KubX);
  k  = j % KubVer;
  f *= KubArr [KubVer-1] [k] * Faktor (j, KubY);

  #undef Faktor

  return (f);
}
/*.BA*/



/*.FE{C 16.4}
     {Newton-Cotes Cubature Formulas for Triangles}
     {Newton-Cotes Cubature Formulas for Triangles}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int Kub3NeC3 (
/*.IX{Kub3NeC3}*/
              REAL Px, REAL Py,
              REAL Qx, REAL Qy,
              REAL Rx, REAL Ry,
              int n,
              REAL f (REAL,REAL),
              REAL* Wert
             );
/***********************************************************************
* Cubature over triangles using the  Newton-Cotes formulas             *
*                                                                      *
* This function integrates f (x,y) over the triangle PQR using the     *
* summed 3 point cubature formulas of Newton-Cotes on sub-triangles.   *
* The sub-triangles are obtained by a regular partition of the edges   *
* of PQR into n equal parts.                                           *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*   REAL   Px,Py     coordinates of P                                  *
*   REAL   Qx,Ry     coordinates of Q                                  *
*   REAL   Rx,Ry     coordinates of R                                  *
*   int    n         partion number along edges                        *
*   REAL   f ()      function                                          *
*   REAL   *Wert     value of integral                                 *
*                                                                      *
* Return value :                                                       *
*   0:               o.k.                                              *
*   1:               n improper                                        *
*   2:               corners  P, Q and R are collinear                 *
*                                                                      *
* Author:            Uli Eggermann, 8.1.1990                           *
.BA*)
***********************************************************************/
/*.BE*/

#define Kub3NeC3Epsilon .000001           /* for collinearity test    */

/***********************************************************************
* Cubature over triangle via summed 3-point Newton-Cotes formula       *
***********************************************************************/

int Kub3NeC3 (REAL Px, REAL Py, REAL Qx, REAL Qy, REAL Rx, REAL Ry,
              int n, REAL f (REAL,REAL), REAL* Wert)
{
   REAL hPQx, hPQy, hPRx, hPRy, Area;
   int  i = 0, j;

   if (n < 1)                          return 1;        /* n ok ?     */
   Area =   Px * (Qy - Ry)
          + Qx * (Ry - Py)
          + Rx * (Py - Qy);                             /* P, Q and R */
   if (FABS (Area) < Kub3NeC3Epsilon ) return 2;        /* collinear? */

   n *= 2;                            /* number of halved edges       */
   hPQx = (Qx - Px) / n;
   hPQy = (Qy - Py) / n;              /* halve the vector PQ          */
   hPRx = (Rx - Px) / n;
   hPRy = (Ry - Py) / n;              /* halve the vector PR          */

   #define  auswerten  (i % 2 || j % 2)            /* i or j even     */
   #define  X          Px + hPQx*i + hPRx*j        /* x-argument      */
   #define  Y          Py + hPQy*i + hPRy*j        /* y-argument      */
   #define  wieoft     ((i==0 || j==0 || i==n-j) ? 1 : 2)
                                /*  do not double at the edges        */

   *Wert = ZERO;                                   /* integral  = 0   */

   for (j = 0; j < n; j++)                        /* j moves along PR */
     for (i = 0; i <= n-j; i++)                   /* i moves along PQ */
       if (auswerten)                              /* sum if needed   */
         *Wert += wieoft * f(X,Y);

   *Wert += wieoft * f(X,Y);                       /*   sum           */
   *Wert *= Area / ((REAL)1.5 * n * n);            /* last factor     */
   return 0;

   #undef Kub3NeC3Epsilon
   #undef auswerten
   #undef X
   #undef Y
   #undef wieoft
}
/*.BA*/



/*.FE{C 16.5}
     {Romberg Cubature for Rectangles (and Triangles)}
     {Romberg Cubature for Rectangular Regions
      (and Triangular Regions)}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int Kub4RoRi (
/*.IX{Kub4RoRi}*/
              REAL  a, REAL b, int Nx,
              REAL  c, REAL d, int Ny,
              int   nST,
              REAL  f (REAL,REAL),
              REAL* Wert,
              REAL* FehlerSch
             )
/***********************************************************************
* Cubature over a rectangle using the summed Romberg-Richardson method *
*                                                                      *
* Intecrate the function f(x,y) over a rectangle using the summed      *
* trapezoidal rule for the Romberg sequence of step sizes. Then we     *
* compute an improved value using Richardson extrapolations.           *
*                                                                      *
* The Richardson sequence is :                                         *
*       { h/2, h/4, h/8, h/16, h/32, h/64, h/128, h/256 .. }           *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   REAL   a,b,Nx   left, right x-end points                           *
*   int    Nx       # of x-intervals                                   *
*   REAL   c,d,Ny   ditto for y                                        *
*   int    Ny       # of y-intervals                                   *
*   int    nST      number of summed trapezoidal cubatures             *
*   REAL   f ()     function                                           *
*   REAL   *Wert    value for integral                                 *
*   REAL   *FehlSch error estimate for  Wert                           *
*                                                                      *
* Return value :                                                       *
*   0:              o.k.                                               *
*   1:              Nx improper                                        *
*   2:              Ny improper                                        *
*   3:              nST improper                                       *
*   4:              rectangle has area zero                            *
*   5:              lack of available memory                           *
*                                                                      *
* subroutine used :                                                    *
*   Kub4RoST (lokal)                                                   *
*                                                                      *
* Author            Uli Eggermann, 8.1.1990                            *
.BA*)
***********************************************************************/
/*.BE*/
{
    int ret;
    REAL   *L;
    void   *vmblock;

    if (Nx  < 1)          return (1);
    if (Ny  < 1)          return (2);
    if (nST < 1)          return (3);
    if (a == b || c == d) return (4);

    vmblock = vminit();
    L = (REAL *)vmalloc(vmblock, VEKTOR, nST, 0);
    if (! vmcomplete(vmblock))
      return 5;

                               /* initialize  Romberg column :        */
                               /* Cubature using trapezoidal rule     */

    Kub4RoST (a, b, Nx, c, d, Ny, nST, f, L);

                                       /* Compute  Richardson columns */

    ret = RoRiExtr (L, nST, 2, Wert, FehlerSch) ? 5 : 0;

    vmfree(vmblock);
    return (ret);
}

static void Kub4RoST (
/*.IX{Kub4RoST}*/
               REAL  a, REAL b, int Nx,
               REAL  c, REAL d, int Ny,
               int   nST,
               REAL  f (REAL,REAL),
               REAL* W
              )
/***********************************************************************
* Initialize  Romberg column: Cubature via trapezoidal rule            *
***********************************************************************/
{
  int    Nab, Ncd, i, j, k;
  REAL   Hab, Hcd, Faktor;

  for (j = 0; j < nST; j++)
  {
    i = 1 << j; Nab = Nx * i;  Hab = (b - a) / Nab;        /* i = 2^j */
                Ncd = Ny * i;  Hcd = (d - c) / Ncd;
    W[j] = ZERO;
    for (i=0; i<=Nab; i++)
      for (k=0; k<=Ncd; k++)
      {
                          Faktor  = ONE;       /* at edge: 1          */
        if (i>0 && i<Nab) Faktor  = TWO;       /* in interior: 2 or 4 */
        if (k>0 && k<Ncd) Faktor *= TWO;

        W[j] += Faktor * f (a + i * Hab,  c + k * Hcd);
      }

    W[j] *= Hab * Hcd * (REAL)0.25;
  }
  return;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int Kub4BuRi (
/*.IX{Kub4BuRi}*/
              REAL a, REAL b, int Nx,
              REAL c, REAL d, int Ny,
              int nST,
              REAL f(REAL,REAL),
              REAL* Wert,
              REAL* FehlerSch
             )
/***********************************************************************
* Cubature over rectangles using the summed Bulirsch-Richardson method *
*                                                                      *
* We compute the integral of f(x,y) using the summed trapezoidal rule  *
* for the Bulirsch sequence of step sizes and use Richardson extra-    *
* polations to find a better approximate value.                        *
*                                                                      *
* The step sizes for the interval length h are :                       *
*         { h/2, h/3, h/4, h/6, h/8, h/12, h/16, h/24, ... }           *
.BE*)
*                                                                      *
* Input parameters :                                                   *
*   REAL   a,b,Nx   left and right end points in x, number of intervals*
*   REAL   c,d,Ny   left and right end points in y, number of intervals*
*   int    nST      number of summed trapezoidal cubatures             *
*   REAL   f (x,y)  Function                                           *
*   REAL   *Wert    value for integral                                 *
*   REAL   *FehlSch error estimate for  Wert                           *
*                                                                      *
* Return value :                                                       *
*   0:              o.k.                                               *
*   1:              lack of memory                                     *
*   2:              Interval number Nx improper                        *
*   3:              Interval number Ny improper                        *
*   4:              nST incorrect                                      *
*   5:              Integration over interval of size zero             *
*                                                                      *
* subroutines used :                                                   *
*   BuNenner        local function for finding denominator of step size*
*   Kub4BuST        Perform cubature                                   *
*   BuRiExtr        Bulirsch-Richardson extrapolation                  *
*                                                                      *
* Author            Uli Eggermann, 8.5.1990                            *
.BA*)
***********************************************************************/
/*.BE*/
/**********************************************************************
       Cubatur over ectangles using the summed
       Bulirsch-Richardson method
**********************************************************************/

{
    REAL   *L;
    int     j;
    void   *vmblock;

    if (Nx  < 1)          return (2);            /* Check input       */
    if (Ny  < 1)          return (3);
    if (nST < 1)          return (4);
    if (a == b || c == d) return (5);

                      /* Initialize Bulirsch column                   */
                      /* Cubature using trapezoidal rule              */

    vmblock = vminit();
    L = (REAL *)vmalloc(vmblock, VEKTOR, nST, 0);
    if (! vmcomplete(vmblock))
      return 1;

    for (j=0; j<nST; j++)
      Kub4BuST (a, b, Nx, c, d, Ny, j, f, L);

                               /*  Bulirsch-Richardson extrapolation */

    j = BuRiExtr (L, nST, 2, Wert, FehlerSch);

    vmfree(vmblock);
    return (j);
}

static void Kub4BuST (
/*.IX{Kub4BuST}*/
               REAL a, REAL b, int Nx,
               REAL c, REAL d, int Ny,
               int j,
               REAL f(REAL,REAL),
               REAL* L
              )
/***********************************************************************
* Initialize  jth  component of the array L with cubature values       *
* from trapezoidal rule                                                *
***********************************************************************/
{
  int nn, i, k, Nab, Ncd;
  REAL   Hab, Hcd, Faktor;

  nn = BuNenner (j);
  Nab = nn * Nx;      Hab = (b - a) / Nab;
  Ncd = nn * Ny;      Hcd = (d - c) / Ncd;
  for (i = 0; i <= Nab; i++)
    for (k = 0; k <= Ncd; k++)
    {
                                 Faktor  = ONE;
      if ( i > 0 && i < Nab )    Faktor  = TWO;
      if ( k > 0 && k < Ncd )    Faktor *= TWO;

      L[j] += Faktor * f (a + i * Hab, c + k * Hcd);
    }
                                          /*                    hx hy */
  L[j] *= Hab * Hcd * (REAL)0.25;         /* Multiply by        ----- */
                                          /*                      4   */
  /* printf ("%"LZP"f \n", L(j)); */
}



/* ------------------------------------------------------------------ */
/*.BA*/

int Kub3RoRi (
/*.IX{Kub3RoRi}*/
              REAL  Px, REAL Py,
              REAL  Qx, REAL Qy,
              REAL  Rx, REAL Ry,
              int   n,
              REAL  f (REAL,REAL),
              REAL* Wert,
              REAL* FehlerSch
             )
/***********************************************************************
* cubature over triangular regions using the summed 3 point formula    *
* and Romberg-Richardson extrapolation:                                *
*                                                                      *
* Function f (x,y) is integrated approximately over the triangle PQR   *
* using the summed 3 point cubature formula on similar subtriangles.   *
* The number of cubatures with halved subtriangle edge length each is  *
* defined by n.                                                        *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*   REAL   Px,Py     coordinaten of P                                  *
*   REAL   Qx,Ry     coordinaten of Q                                  *
*   REAL   Rx,Ry     coordinaten of R                                  *
*   int    n         number of cubatures                               *
*   REAL   f (x,y)   function to be integrated                         *
*                                                                      *
* Output parameters:                                                   *
*   REAL   *Wert     approximation of the double integral              *
*   REAL   *FehlSch  error estimate for `Wert'                         *
*                                                                      *
* Return value:                                                        *
*   0:               o.k.                                              *
*   1:               n improper                                        *
*   2:               corners P, Q and R are collinear                  *
*   3:               no more memory for auxiliary vector               *
*                                                                      *
* Subroutines used:                                                    *
*   RoRiExtr      Romberg-Richardson extrapolation                     *
*   Kub3NeC3      computation of cubature value                        *
*   Kub3Nec3n     computation and storing of the n cubature values     *
*                                                                      *
* author         Uli Eggermann, 07.05.1990                             *
.BA*)
***********************************************************************/
/*.BE*/
{
  int     i;
  REAL   *W;
  void   *vmblock;

  if (n < 1)                                            return (1);
  vmblock = vminit();
  W = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
    return 3;

               /* storing Newton-Cotes values in the first column     */
  i = Kub3Nec3n (Px,Py, Qx,Qy, Rx,Ry, n, f, W);

                  /* computing the resulting Richardson columns,      */
  if (i == 0)     /* but only if the first column is ok.              */
    if (n == 1)
    {
      *Wert = *W;
      *FehlerSch = ZERO;
      i = 0;
    }
    else i = RoRiExtr (W, n, 2, Wert, FehlerSch);

  vmfree(vmblock);                                  /* re-allocation  */
  return i;
}

static int Kub3Nec3n (
/*.IX{Kub3Nec3n}*/
               REAL Px, REAL Py,
               REAL Qx, REAL Qy,
               REAL Rx, REAL Ry,
               int n,
               REAL f(REAL,REAL),
               REAL* W
              )
/***********************************************************************
* computing and storing the n cubature values in given array W[].      *
***********************************************************************/

{
  int i, ret = 0;
  if (n < 1) return (1);

  for (i = 0; i < n; i++)
  {
    if ((ret = Kub3NeC3 (Px,Py,Qx,Qy,Rx,Ry, 1<<i, f, &W[i])) != 0)
       break;
    #ifdef _TEST_
      printf ("%f ", W[i]);                             /* test print */
    #endif
  }
  #ifdef _TEST_
    printf ("\n");                                      /* test print */
  #endif
  return (ret);
}

static int BuRiExtr (REAL* BuFo,
/*.IX{BuRiExtr}*/
              int   nBuFo,
              int   Ordnung,
              REAL* Wert,
              REAL* FehlerSch
             )
/***********************************************************************
* Richardson extrapolation for a given Bulirsch sequence.              *
*                                                                      *
* Parameters:                                                          *
*                                                                      *
*   REAL   BuFo []     Bulirsch sequence                               *
*   int    nBuFo       lenght of BuFo                                  *
*   int    Ordnung     order of method                                 *
*   REAL   *Wert       final value of extrapolation                    *
*   REAL   *FehlerSch  error estimate for  Wert                        *
*                                                                      *
* Return value :                                                       *
*                                                                      *
*   0:                 o.k.                                            *
*   1:                 lack of memory                                  *
*   2:                 nBuFo too small                                 *
*                                                                      *
* subroutines used :                                                   *
*                                                                      *
*   BuNenner ()  denominators of the step sizes, common to the Bulirsch*
*                sequence and Richardson extrapolation.                *
*                                                                      *
* author         Uli Eggermann, 03.31.1996                             *
***********************************************************************/
{
  int  j, k;
  REAL P, *RiEx;
  void *vmblock;

  if (nBuFo < 2) return 2;
  vmblock = vminit();
  RiEx = (REAL *)vmalloc(vmblock, VEKTOR, nBuFo, 0);
  if (! vmcomplete(vmblock))
    return 1;
  for (k=0; k<nBuFo; k++)  RiEx[k] = BuFo[k];

  for (k=1; k<nBuFo; k++)
    for (j=0; j<nBuFo-k; j++)
    {
      P = POW ((REAL) BuNenner (j+k) / (REAL) BuNenner (j),
               (REAL) k * Ordnung);
      RiEx[j] = (P * RiEx [j+1] - RiEx [j]) / (P - ONE);
    }

  *Wert      = RiEx [0];
  *FehlerSch = RiEx [0] - RiEx [1];

  vmfree(vmblock);
  return 0;
}

static int BuNenner (int j)
/*.IX{BuNenner}*/
/***********************************************************************
* Compute denominators of step sizes                                   *
* (depends on parity of j)                                             *
***********************************************************************/
{
  if (j == 0)     return     1;
  if (j % 2 == 0) return 3 * 1 << ((j + 2) / 2);        /* j even     */
                  return     1 << ((j + 1) / 2);        /* j odd      */
}

static int RoRiExtr (
/*.IX{RoRiExtr}*/
              REAL* RoFo,
              int   nRoFo,
              int   Ordnung,
              REAL* Wert,
              REAL* FehlerSch
             )
/***********************************************************************
* Richardson extrapolation for a given Romberg sequence                *
*                                                                      *
* Parameter:                                                           *
*   REAL  RoFo []    given  Romberg sequence                           *
*   int   nRoFo      length of  RoFo                                   *
*   int   Ordnung    Order of method                                   *
*   REAL* Wert       final value of extrapolation                      *
*   REAL* FehlerSch  error estimate of  Wert                           *
*                                                                      *
* Return value :                                                       *
*   0:              o.k.                                               *
*   1:              lack of available memory                           *
*   2:              nRoFo too small                                    *
*                                                                      *
* Author          Uli Eggermann, 3.31.1996                             *
***********************************************************************/
{
  int   j, k, p, s;
  REAL* RiEx;
  void *vmblock;

  if (nRoFo < 2) return (2);
  vmblock = vminit();
  RiEx = (REAL *)vmalloc(vmblock, VEKTOR, nRoFo, 0);
  if (! vmcomplete(vmblock))
    return 1;

  for (k = 0; k < nRoFo; k++)                                 /* copy */
    RiEx[k] = RoFo[k];

  p = s = 1 << Ordnung;                                  /* 2 ^ order */

  for (k = nRoFo - 1;   k > 0;   p *= s, k--)
    for (j = 0; j < k; j++)
      RiEx[j] = (p * RiEx[j+1] - RiEx[j]) / (p - 1);

  *Wert      = RiEx[0];
  *FehlerSch = RiEx[0] - RiEx[1];

  vmfree(vmblock);

  return (0);
}

/* -------------------------- END kubnec.c -------------------------- */
