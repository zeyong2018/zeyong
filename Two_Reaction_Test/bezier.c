#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODUL bezier.c ------------------------- */

#include <basis.h>   /*  for  REAL, ZERO, ONE, THREE, sqr, FABS, TWO, */
                     /*       FOUR, SQRT, MACH_EPS                    */
#include <bikub.h>   /*  for  kubbez, valbez, bezier, rechvp, rechwp  */
/*.BA*/



/*.FE{C 12.3.1}{B\'ezier Spline Curves}{B\'ezier Spline Curves}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int kubbez     /* compute Bezier points of a bezier spline curve .....*/
/*.IX{kubbez}*/
        (
         REAL   *b[],                /* weight points ................*/
         REAL   *d[],                /* Bezier points ................*/
         int    m,                   /* number of curve pieces .......*/
         int    dim                  /* 2,3 for planar, spatial curve */
        )                            /* error code ...................*/

/***********************************************************************
* Computes Bezier points using the cubic Bezier method.                *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*    REAL d[][3]           coordinates of weight points                *
*    int  m                number of curve pieces                      *
*    int  dim              = 2: planar curve                           *
*                          = 3: spatial curve                          *
*                                                                      *
* Output parameter:                                                    *
*                                                                      *
*    REAL b[][3]           coordinates of the Bezier points            *
*                                                                      *
* Return value:                                                        *
*                                                                      *
*   Error code. The following values are possible:                     *
*   = 0: all is ok                                                     *
*   = 1: invalid input parameters:                                     *
*        m < 2  or  dim < 2  or  dim > 3                               *
*                                                                      *
* Global names used:                                                   *
*                                                                      *
*   REAL, TWO, THREE, FOUR, SIX                                        *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/

{
  int i, k;

  if (m < 2 || dim < 2 || dim > 3)
    return 1;

  for (i = 0; i < dim; i++)
  {
    for (k = 1; k < m; k++)
    {
      b[3*k-2][i] = (TWO*d[k-1][i] +  d[k][i]                ) / THREE;
      b[3*k]  [i] = (d[k-1][i] + FOUR*d[k][i] +     d[k+1][i]) / SIX;
      b[3*k+2][i] = (                 d[k][i] + TWO*d[k+1][i]) / THREE;
    }
    b[  2  ][i]   = (                 d[0][i] + TWO*d[ 1 ][i]) / THREE;
    b[3*m-2][i]   = (TWO*d[m-1][i] +  d[m][i]                ) / THREE;

    b[  0  ][i] = d[0][i];                /* set up end points        */
    b[ 3*m ][i] = d[m][i];                /* for natural cubic spline */
  }

  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int valbez     /* evaluation of a Bezier spline curve ................*/
/*.IX{valbez}*/
        (
         REAL   t,                   /* parameter value t from [0,1]  */
         int    m,                   /* number of curve pieces        */
         int    dim,                 /* 2,3 for planar, spatial curve */
         REAL   *b[],                /* Bezier points                 */
         REAL   *x,                  /* coordinates of curve point    */
         REAL   *y,
         REAL   *z
        )                            /* error code                    */

/* ================================================================== */
/*   `valbez' computes the cartesian coordinates (x,y,z) of a point   */
/*   on the Bezier curve parametrized by t from [0,1]. t=0 renders    */
/*   the starting point and t=1 the end point of the curve.           */
/*  ================================================================  */
/*.BE*/
/*   Input parameters:                                                */
/*                                                                    */
/*    Name    Type               Meaning                              */
/*   ---------------------------------------------------------------  */
/*    t       REAL               curve parameter from [0,1]           */
/*    m       int                number of curve segments             */
/*    dim     int                = 2: planar curve                    */
/*                               = 3: spatial curve                   */
/*    b       REAL **            Bezier points (output from kubbez()) */
/*                                                                    */
/*   Output parameters:                                               */
/*                                                                    */
/*    Name          Type         Meaning                              */
/*   ---------------------------------------------------------------  */
/*    x, y, z       REAL         cartesian coordinates of the point   */
/*                                                                    */
/*   Return value:                                                    */
/*     = 0 : all ok, coordinates computed                             */
/*     = 1 : m < 2                                                    */
/*     = 2 : t not from the interval [0, 1]                           */
/*                                                                    */
/*    Global names used:                                              */
/*      REAL, ZERO, ONE, sqr, THREE                                   */
/*                                                                    */
/* ================================================================== */

{
  int  k3;
  REAL tt, v;

  if (m < 2 || dim < 2 || dim > 3)
    return 1;

  if (t < ZERO || t > ONE)
    return 2;

  tt = THREE * t * (REAL)m;
  k3 = (int)(tt / THREE) * 3;
  if (k3 == 3 * m)               /* curve piece m???                  */
    k3 -= 3;                     /* doesn't exist, take previous one! */
  t = (tt - k3) / THREE;
  v = ONE - t;

  *x   = sqr(v) * (v * b[k3][0]  +  THREE * t * b[k3+1][0]) +
         sqr(t) * (THREE * v * b[k3+2][0] + t * b[k3+3][0]);
  *y   = sqr(v) * (v * b[k3][1]  +  THREE * t * b[k3+1][1]) +
         sqr(t) * (THREE * v * b[k3+2][1] + t * b[k3+3][1]);
  if (dim == 3)
    *z = sqr(v) * (v * b[k3][2]  +  THREE * t * b[k3+1][2]) +
         sqr(t) * (THREE * v * b[k3+2][2] + t * b[k3+3][2]);

  return 0;
}
/*.BA*/



/*.FE{C 12.3.2}{B\'ezier Spline Surfaces}
               {B\'ezier Spline Surfaces}*/

/*.BE*/
/* ------------------------------------------------------------------ */

static void b_point (REAL*** b,
/*.IX{b\unt point}*/
                     REAL*** d,
                     int     m,
                     int     n)
/***********************************************************************
* Computes the unknown Bezier points for a bicubic Bezier spline.      *
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*   REAL  b [3*m+1][3*n+1][3]  coordinates of the Bezier points:       *
*                              (see Figure 12.4)                       *
*   REAL  d [m+1][n+1][3]      coordinates of the weights              *
*   int   m                    number of patches in one direction      *
*   int   n                    ditto in other direction                *
*                                                                      *
* Output parameter:                                                    *
*                                                                      *
*   REAL  b [3*m+1][3*n+1][3]  coordinates of all Bezier points        *
***********************************************************************/
{
  int i, j, k;

  for (k=0; k<3; k++)
  {
    for (i=1; i<=m; i++)
      for (j=1; j<=n; j++)
        b [3*i-2][3*j-2][k] = (4.*d[i-1][j-1][k] + 2.*d[i-1][j][k] +
                               2.*d[ i ][j-1][k] +    d[ i ][j][k])/9.;
    for (i=0; i<=m-1; i++)
      for (j=1; j<=n; j++)
        b [3*i+2][3*j-2][k] = (4.*d[i+1][j-1][k] + 2.*d[i][j-1][k] +
                               2.*d[i+1][ j ][k] +    d[i][ j ][k])/9.;
    for (i=1; i<=m; i++)
      for (j=0; j<=n-1; j++)
        b [3*i-2][3*j+2][k] = (4.*d[i-1][j+1][k] + 2.*d[i-1][j][k] +
                               2.*d[ i ][j+1][k] +    d[ i ][j][k])/9.;
    for (i=0; i<=m-1; i++)
      for (j=0; j<=n-1; j++)
        b [3*i+2][3*j+2][k] = (4.*d[i+1][j+1][k] + 2.*d[i][j+1][k] +
                               2.*d[i+1][ j ][k] +    d[i][ j ][k])/9.;
    for (i=1; i<=m; i++)
      for (j=1; j<=n-1; j++)
        b [3*i-2][3*j][k] = (2.*d[i-1][j-1][k] + 8.*d[i-1][ j ][k] +
                                d[ i ][j-1][k] + 2.*d[i-1][j+1][k] +
                             4.*d[ i ][ j ][k] +    d[ i ][j+1][k])/18.;
    for (i=1; i<=m-1; i++)
      for (j=1; j<=n; j++)
        b [3*i][3*j-2][k] = (2.*d[i-1][j-1][k] + 8.*d[ i ][j-1][k] +
                                d[i-1][ j ][k] + 2.*d[i+1][j-1][k] +
                             4.*d[ i ][ j ][k] +    d[i+1][ j ][k])/18.;
    for (i=1; i<=m-1; i++)
      for (j=0; j<=n-1; j++)
        b [3*i][3*j+2][k] = (2.*d[i-1][j+1][k] + 8.*d[ i ][j+1][k] +
                                d[i-1][ j ][k] + 2.*d[i+1][j+1][k] +
                             4.*d[ i ][ j ][k] +    d[i+1][ j ][k])/18.;
    for (i=0; i<=m-1; i++)
      for (j=1; j<=n-1; j++)
        b [3*i+2][3*j][k] = (2.*d[i+1][j-1][k] + 8.*d[i+1][ j ][k] +
                                d[ i ][j-1][k] + 2.*d[i+1][j+1][k] +
                             4.*d[ i ][ j ][k] +    d[ i ][j+1][k])/18.;
    for (i=1; i<=m-1; i++)
      for (j=1; j<=n-1; j++)
        b [3*i][3*j][k] = (     d[i-1][j-1][k] + 4.*d[ i ][j-1][k] +
                                d[i+1][j-1][k] + 4.*d[i-1][ j ][k] +
                            16.*d[ i ][ j ][k] + 4.*d[i+1][ j ][k] +
                                d[i-1][j+1][k] + 4.*d[ i ][j+1][k] +
                                d[i+1][j+1][k]                    )/36.;
  }
}

static void intpol (REAL*   diff,
/*.IX{intpol}*/
                    int     i,
                    int     j,
                    REAL*** b,
                    int     m,
                    int     n)
/***********************************************************************
* Changes the interpolation points after the Bezier method has computed*
* the bicubic spline surface.                                          *
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*    REAL  diff [3]             coordinates of the difference vector,  *
*                               used to shift the Bezier surface       *
*    int   i, j                 Indices for the patch around which the *
*                               Bezier surface is altered              *
*    int   m                    number of patches in first direction   *
*    int   n                    diito for other direction              *
*    REAL  b [3*m+1][3*n+1][3]  coordinated of the Bezier points       *
*                                                                      *
* Output parameters:                                                   *
*                                                                      *
*    REAL  b [3*m+1][3*n+1][3]  Coordinates of the Bezier points       *
***********************************************************************/
{
  static REAL gewicht[7][7] =
                {
                 { 0.0625, 0.125, 0.25, 0.25, 0.25, 0.125, 0.0625 },
                 {  0.125,  0.25,  0.5,  0.5,  0.5,  0.25,  0.125 },
                 {   0.25,   0.5,  1.0,  1.0,  1.0,   0.5,   0.25 },
                 {   0.25,   0.5,  1.0,  1.0,  1.0,   0.5,   0.25 },
                 {   0.25,   0.5,  1.0,  1.0,  1.0,   0.5,   0.25 },
                 {  0.125,  0.25,  0.5,  0.5,  0.5,  0.25,  0.125 },
                 { 0.0625, 0.125, 0.25, 0.25, 0.25, 0.125, 0.0625 }
                };
  REAL        linker_rand [3][7],            /* saved boundary points */
              rechter_rand[3][7],            /* of the bezier matrix  */
              unterer_rand[3][7],
              oberer_rand [3][7],
              tmp;
  int         k1, k2, l,                     /* loop variables        */
              i3, j3;                        /* 3*i and 3*j           */

  i3 = 3 * i;
  j3 = 3 * j;

  if (i == 1 || i == m-1 ||       /* boundary points to be destroyed? */
      j == 1 || j == n-1)
    for (l = 0; l < 3; l++)                        /* save            */
      for (k1 = -3; k1 <= 3; k1++)                 /* all 28 possibly */
        unterer_rand[l][3+k1] = b[i3+k1][0][l],    /* affected        */
        oberer_rand [l][3+k1] = b[i3+k1][3*n][l],  /* boundary points */
        linker_rand [l][3+k1] = b[0]   [j3+k1][l],
        rechter_rand[l][3+k1] = b[3*m] [j3+k1][l];

  for (l = 0; l < 3; l++)                                /* move      */
    for (tmp = diff[l], k1 = -3; k1 <= 3; k1++)          /* 45 points */
      for (k2 = -3; k2 <= 3; k2++)
        b[i3+k1][j3+k2][l] += tmp * gewicht[3+k1][3+k2];

  if (i == 1 || i == m-1 ||             /* boundary points destroyed? */
      j == 1 || j == n-1)
    for (l = 0; l < 3; l++)                        /* restore         */
      for (k1 = -3; k1 <= 3; k1++)                 /* all 28 possibly */
        b[i3+k1][0][l]    = unterer_rand[l][3+k1], /* destroyed       */
        b[i3+k1][3*n][l]  = oberer_rand [l][3+k1], /* boundary points */
        b[0]   [j3+k1][l] = linker_rand [l][3+k1],
        b[3*m] [j3+k1][l] = rechter_rand[l][3+k1];
}
/*.BA*/

int bezier (REAL*** b,
/*.IX{bezier}*/
            REAL*** d,
            int     modified,
            int     m,
            int     n,
            REAL    eps)

/***********************************************************************
* The bicubic and modified bicubic Bezier spline surfaces.             *
* This algorithm computes interpolation points from the input data for *
* a spline surface according to the bicubic Bezier spline method.      *
* When using the modified method, the given interpoation points are    *
* used as weight points first, for which the pseudo-interpolation      *
* points are computed. These are then altered until they agree with    *
* true interpolation points to within eps.                             *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*    REAL  b [3*m+1][3*n+1][3]  Coordinates of the Bezier points:      *
*                               The following entries must be used:    *
*                                 b [i][j][k] for k=0,1,2 and          *
*                                  i=0,...,3*m   and j=0;              *
*                                  i=0           and j=0 (1) 3*n;      *
*                                  i=0,...,3*m   and j=      3*n;      *
*                                  i=      3*m   and j=0 (1) 3*n.      *
*                               In addition we need the entries        *
*                               i=3,6,...,3*m-6,3*m-3 and j=3,6,...,   *
*                                                      .. ,3*n-6,3*n-3 *
*                               for the modified method.               *
*    REAL  d [m+1][n+1][3]      modified  = 0: coordinates of weight   *
*                                              points                  *
*                               modified != 0: not used                *
*    int   modified             modified  = 0: Bezier method           *
*                               modified != 0: modified Bezier method  *
*    int   m                    number of patches in first direction   *
*    int   n                    ditto in second direction              *
*    REAL  eps                  modified = 1 : accuracy bound for      *
*                                              interpolation           *
*                                                                      *
* Output parameters:                                                   *
*                                                                      *
*  REAL    b [3*m+1][3*n+1][3]  Coordinates of the Bezier points       *
*                               b [i][j][k]                            *
*                               for i=0,...,3*m, j=0,...,3*n, k=0,1,2  *
*                                                                      *
*                                                                      *
* Return value:                                                        *
*                                                                      *
*  = 0: no error                                                       *
*  = 1: m < 2  or  n < 2                                               *
*  = 2: eps too small (only modified method)                           *
*                                                                      *
* Subroutines used:         intpol, b_point                            *
*                                                                      *
* Macros used:              FABS                                       *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i, j, l;
  REAL diff [3];

  if (m < 2 || n < 2)
    return 1;

  if (modified)
  {

    int okay = 0;

    if (eps < (REAL)128.0 * MACH_EPS)
      return 2;

    for (l=0; l<3; l++)
      for (i=0; i<=m; i++)
        for (j=0; j<=n; j++)
          d [i][j][l] = b [3*i][3*j][l];
    b_point (b, d, m, n);

    while (!okay)
    {
      for (i=1; i<=m-1; i++)
        for (j=1; j<=n-1; j++)
        {
          for (l=0; l<3; l++)
            diff [l] = d [i][j][l] - b [3*i][3*j][l];
          intpol (diff, i, j, b, m, n);
        }
      for (okay=1,i=1; i<=m-1; i++)
        for (j=1; j<=n-1; j++)
          for (l=0; l<3; l++)
            if (FABS (d [i][j][l] - b [3*i][3*j][l]) > eps)
              okay = 0;
    } /* while (!okay) */

  }    /* (modified) */

  else /* (not modified) */
    b_point (b, d, m, n);

  return 0;
}

static void rechp (REAL*** b, int m, int n,
/*.IX{rechp}*/
            REAL vp, REAL wp, REAL* point)
/***********************************************************************
* Compute the spacial coordinates of the point of intersection of two  *
* parameter lines of a Bezier spline.                                  *
* ==================================================================== *
*                                                                      *
*   Input parameters:                                                  *
*   -----------------                                                  *
*                                                                      *
*    Name    Type/size                 Meaning                         *
*   ------------------------------------------------------------------ *
*    b       REAL  /[3*m+1][3*n+1][3]  coordinates of Bezier points    *
*                                      b is a pointer array            *
*    m       int/---                   number of patches in first      *
*                                      direction                       *
*    n       int/---                   ditto for second direction      *
*                                      Richtung                        *
*    vp,wp   REAL  /---                parameter lines of the problem  *
*                                                                      *
*                                                                      *
*   Output parameters:                                                 *
*   ------------------                                                 *
*                                                                      *
*    Name    Type/size                 Meaning                         *
*   ------------------------------------------------------------------ *
*    point   REAL  /[3]                coordinates of intersection     *
*                                      point on  Bezier surface        *
*                                                                      *
***********************************************************************/
{
  int  i, j, k;
  REAL h, h1, h2, h3, h4, h5, h6, h7, h8, v, w, vv, ww;

  vv = vp * (3 * n);                ww = wp * (3 * m);
  i  = (int) (vv / 3.) * 3;         j  = (int) (ww / 3.) * 3;
  if (i >= 3*n) i = 3 * (n-1);      if (j >= 3*m) j = 3 * (m-1);
  v  = (vv - i) / 3.;               w  = (ww - j) / 3.;

  h  = 1 - v;                       h1 =      h * h * h;
                                    h2 = 3. * h * h * v;
                                    h3 = 3. * h * v * v;
                                    h4 =      v * v * v;

  h  = 1 - w;                       h5 =      h * h * h;
                                    h6 = 3. * h * h * w;
                                    h7 = 3. * h * w * w;
                                    h8 =      w * w * w;

  for (k=0; k<=2; k++)
  {
    point [k]  = (b [ j ][ i ][k] * h1 + b [ j ][i+1][k] * h2 +
                  b [ j ][i+2][k] * h3 + b [ j ][i+3][k] * h4  ) * h5;
    point [k] += (b [j+1][ i ][k] * h1 + b [j+1][i+1][k] * h2 +
                  b [j+1][i+2][k] * h3 + b [j+1][i+3][k] * h4  ) * h6;
  }
  for (k=0; k<=2; k++)
  {
    point [k] += (b [j+2][ i ][k] * h1 + b [j+2][i+1][k] * h2 +
                  b [j+2][i+2][k] * h3 + b [j+2][i+3][k] * h4  ) * h7;
    point [k] += (b [j+3][ i ][k] * h1 + b [j+3][i+1][k] * h2 +
                  b [j+3][i+2][k] * h3 + b [j+3][i+3][k] * h4  ) * h8;
  }
  return;
}
/*.BA*/

int rechvp (REAL*** b, int m, int n, REAL vp,
/*.IX{rechvp}*/
            int num, REAL *points[])

/***********************************************************************
* Compute the coordinates of num points of a Bezier surface that lie on*
* the parameter line vp.                                               *
* (vp=0  if i=0; vp=1 if i=3*n; i.e. vp acts as a scale of the (m x n) *
* patches in the second direction  n).                                 *
.BE*)
* ==================================================================== *
*                                                                      *
*   Input parameters:                                                  *
*   -----------------                                                  *
*                                                                      *
*    Name    Type/size                 Meaning                         *
*   ------------------------------------------------------------------ *
*    b       REAL  /[3*m+1][3*n+1][3]  coordinates of Bezier points, an*
*                                      array of pointers               *
*    m       int/---                   number of patches in first      *
*                                      direction                       *
*    n       int/---                   ditto in second direction       *
*    vp      REAL  /---                defines the parameter line      *
*    num     int/---                   number of points desired        *
*                                                                      *
*                                                                      *
*   Output parameters:                                                 *
*   ------------------                                                 *
*                                                                      *
*    Name    Type/size                 Meaning                         *
*   ------------------------------------------------------------------ *
*    points  REAL  /[num][3]           coordinates of the computed     *
*                                      points                          *
*                                                                      *
*                                                                      *
*   Return value:                                                      *
*   -------------                                                      *
*                                                                      *
*   = 0: all is ok                                                     *
*   = 1: m < 2  or  n < 2                                              *
*   = 2: vp not from [0,1]                                             *
*   = 3: num < 2                                                       *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   subroutine used :                 rechp                            *
*   -----------------                                                  *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int    i, k;
  REAL   step, h, point [3];

  if (m < 2 || n < 2)
    return 1;

  if (vp < ZERO || vp > ONE)
    return 2;

  if (num < 2)
    return 3;

  h = (REAL)(num - 1);
  for (i = 0; i <= num-1; i++)
  {
    step = i / h;
    rechp (b, m, n, vp, step, point);
    for (k = 0; k <= 2; k++)
      points [i][k] = point [k];
  }
  return 0;
}
/*.BA*/

int rechwp (REAL*** b, int m, int n, REAL wp,
/*.IX{rechwp}*/
            int num, REAL *points[])

/***********************************************************************
* Compute the coordinates of num points of a Bezier surface that lie on*
* the parameter line wp.                                               *
* (wp=0  if j=0; wp=1 if j=3*m; i.e. wp acts as a scale of the (m x n) *
* patches in the first direction  m).                                 *
.BE*)
* ==================================================================== *
*                                                                      *
*   Input parameters:                                                  *
*   -----------------                                                  *
*                                                                      *
*    Name    Type/size                 Meaning                         *
*   ------------------------------------------------------------------ *
*    b       REAL  /[3*m+1][3*n+1][3]  coordinates of Bezier points, an*
*                                      array of pointers               *
*    m       int/---                   number of patches in first      *
*                                      direction                       *
*    n       int/---                   ditto in second direction       *
*    wp      REAL  /---                defines the parameter line      *
*    num     int/---                   number of points desired        *
*                                                                      *
*                                                                      *
*   Output parameters:                                                 *
*   ------------------                                                 *
*                                                                      *
*    Name    Type/size                 Meaning                         *
*   ------------------------------------------------------------------ *
*    points  REAL  /[num][3]           coordinates of the computed     *
*                                      points                          *
*                                                                      *
*                                                                      *
*   Return value:                                                      *
*   -------------                                                      *
*                                                                      *
*   = 0: all is ok                                                     *
*   = 1: m < 2  or  n < 2                                              *
*   = 2: wp not from [0,1]                                             *
*   = 3: num < 2                                                       *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   subroutine used :                 rechp                            *
*   -----------------                                                  *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int i, k;
  REAL   step, h, point [3];

  if (m < 2 || n < 2)
    return 1;

  if (wp < ZERO || wp > ONE)
    return 2;

  if (num < 2)
    return 3;

  h = (REAL  ) (num - 1);
  for (i=0; i<=num-1; i++)
  {
    step = (REAL  ) (i) / h;
    rechp (b, m, n, step, wp, point);
    for (k=0; k<=2; k++)
      points [i][k] = point [k];
  }
  return 0;
}
/*.BA*/



/*.FE{C 12.3.3}
     {Modified Interpolating Cubic B\'ezier Splines}
     {Modified Interpolating Cubic B\'ezier Splines}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int mokube     /* compute Bezier points on a interpolating curve .....*/
/*.IX{mokube}*/
        (
         REAL   *b[],                /* weight points ................*/
         REAL   *d[],                /* Bezier points ................*/
         int    m,                   /* number of spline segments ....*/
         int    dim,                 /* 2,3 for planar, spatial curve */
         REAL   eps                  /* accuracy of interpolation ....*/
        )                            /* error code ...................*/

/***********************************************************************
* The modified cubic  Bezier method :                                  *
* The given interpolation points are used as weight points for which   *
* pseudo-interpolation points are computed.                            *
* These are altered until they agree with the interpolation points to  *
* within the desired accuracy.                                         *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*    REAL d[m+1][3]     coordinates of the Bezier points               *
*                         (used as weights)                            *
*    int  m             number of spline segments                      *
*    int  dim           = 2: planar curve                              *
*                       = 3: spatial curve                             *
*    REAL eps           desired accuracy                               *
*                                                                      *
* Output parameter:                                                    *
*                                                                      *
*    REAL b[3*m+1][3]   coordinates of the Bezier points               *
*                                                                      *
* Return value:                                                        *
*                                                                      *
*   Error code. The following values are possible:                     *
*   = 0: all is ok                                                     *
*   = 1: invalid input parameters:                                     *
*        m < 2  or  dim < 2  or  dim > 3                               *
*   = 2: eps too small                                                 *
*                                                                      *
* Global names used:                                                   *
*                                                                      *
*   REAL, MACH_EPS, kubbez, TWO, FOUR, SQRT                            *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i, k, okay = 0, fehler;

  if (eps < (REAL)128.0 * MACH_EPS)
    return 2;

  fehler = kubbez(b, d, m, dim);     /* compute missing Bezier points */

  if (fehler)
    return fehler;

  while (! okay)      /* the computed Bezier points are adjusted      */
                      /* until they match the given Bezier points to  */
                      /* within eps                                   */
  {
     for (i = 0; i < dim; i++)
     {
       for (k = 3; k <= 3 * m - 3; k += 3)
       {
         REAL diff  = d[k/3][i] - b[k][i],
              diff2 = diff / TWO,
              diff4 = diff / FOUR;
         if (k != 3)     b[k-3][i] += diff4;
                         b[k-2][i] += diff2;
                         b[k-1][i] += diff;
                         b[ k ][i] += diff;
                         b[k+1][i] += diff;
                         b[k+2][i] += diff2;
         if (k != 3*m-3) b[k+3][i] += diff4;
       }
     }
    for (k = 1; k < m; k++)
    {
      REAL diff = ZERO;
      for (okay = 1, i = 0; i < dim; i++)
        diff +=   (d[k][i] - b[3*k][i])
                * (d[k][i] - b[3*k][i]);
      if (SQRT(diff) > eps)
      {
        okay = 0;
        break;
      }
    }
  }

  return 0;
}

/* --------------------------- END bezier.c ------------------------- */
