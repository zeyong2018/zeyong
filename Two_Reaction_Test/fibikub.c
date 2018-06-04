#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 16.8}
     {Riemann Double Integrals using Bicubic Splines}
     {Riemann Double Integrals using Bicubic Splines}*/

/*.BE*/
/* ------------------------ MODULE fibikub.c ------------------------ */

#include <basis.h>                   /*  for  REAL, ZERO, ONE, mat4x4 */
#include <bikub.h>                   /*  for  xyintv                  */
#include <kubatur.h>                 /*  for  fibiku, fibik2          */

static REAL fibik1 (mat4x4** a, int i, int j, REAL xi, REAL yj);
/*.BA*/

REAL fibiku (int      n,
/*.IX{fibiku}*/
             int      m,
             mat4x4** a,
             REAL*    x,
             REAL*    y
            )
/***********************************************************************
* Computes the Riemann surface integral over hte domain of a bicubic   *
* spline  (x [0] to x [n], y [0] to y [n]).                            *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*   int  n                   number of x-intervals                     *
*   int  m                   number of y-intervals                     *
*   REAL a [n+1][m+1][4][4]  aray of pointers:                         *
*                            the coefficients of the bicubic spline    *
*   REAL x [n+1]             end points of  x-intervals                *
*   REAL y [m+1]             end points of  y-intervals                *
*                                                                      *
* Return value :             value of Riemann integral                 *
*                                                                      *
* subroutines used :         fibik1                                    *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int    i, j;
  REAL   value;

  for (value=ZERO, i=0; i<=n-1; i++)
    for (j=0; j<=m-1; j++)
      value += fibik1 (a, i, j, x[i+1] - x[i], y[j+1] - y[j]);

  return (value);
}
/*.BA*/

int fibik2 (int      n,
/*.IX{fibik2}*/
            int      m,
            mat4x4** a,
            REAL*    x,
            REAL*    y,
            REAL     xlow,
            REAL     ylow,
            REAL     xup,
            REAL     yup,
            REAL*    value
           )
/***********************************************************************
* Computes the double integral for a spline over the rectangle         *
* xlow to xup, ylow to yup.                                            *
* The individual partitions neen not coincide with the nodes for the   *
* bicubic splin.                                                       *
.BE*)
*                                                                      *
* Input parameters:    see  fibiku                                     *
*                                                                      *
* Output parameter:                                                    *
*                                                                      *
*   REAL value         integral value                                  *
*                                                                      *
* Return value :                                                       *
*                                                                      *
*     = 0 : no error                                                   *
*    != 0 : error in xyintv                                            *
*                                                                      *
* subroutines used :   xyintv, fibik1                                  *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int   ilow,  iup,  jlow,  jup, error, i, j;
  REAL xilow, xiup, yjlow, yjup, factor, xi;

  error = xyintv (n, m, x, y, &ilow, &jlow, &xilow, &yjlow, xlow, ylow);
  if (error) return error;

  error = xyintv (n, m, x, y, &iup, &jup, &xiup, &yjup, xup, yup);
  if (error) return error+1;

  factor = ONE;

  if (ilow > iup)
  {
     i     =  iup;  iup =  ilow;   ilow =  i;
    xi     = xiup; xiup = xilow;  xilow = xi;
    factor = -factor;
  }
  if (jlow > jup)
  {
     i     =  jup;  jup =  jlow;   jlow =  i;
    xi     = yjup; yjup = yjlow;  yjlow = xi;
    factor = -factor;
  }
  *value  =   fibik1 (a, ilow, jlow, xilow, yjlow)
            - fibik1 (a, ilow, jup,  xilow, yjup)
            - fibik1 (a, iup,  jlow, xiup,  yjlow)
            + fibik1 (a, iup,  jup,  xiup,  yjup);

  for (i=ilow; i<=iup-1; i++)
    *value +=   fibik1 (a, i, jup,  x[i+1]-x[i], yjup)
              - fibik1 (a, i, jlow, x[i+1]-x[i], yjlow);

  for (j=jlow; j<=jup-1; j++)
  {
    *value +=   fibik1 (a, iup,  j, xiup,  y[j+1]-y[j])
              - fibik1 (a, ilow, j, xilow, y[j+1]-y[j]);
    for (i=ilow; i<=iup-1; i++)
      *value += fibik1 (a, i, j, x[i+1]-x[i], y[j+1]-y[j]);
  }

  *value *= factor;

  return 0;

} /* fibik2 */

static REAL fibik1 (mat4x4** a,
/*.IX{fibik1}*/
             int      i,
             int      j,
             REAL     xi,
             REAL     yj
            )
/***********************************************************************
* Computes the double integral for a spline function on a rectangle    *
* 0 to xi and 0 to yj. (xi and yj are relative coordinates)            *
*                                                                      *
* Input parameters:           see  fibiku                              *
*                                                                      *
* Return value :              integral value                           *
*                                                                      *
***********************************************************************/
{
  int  k, l;
  REAL xip [5], yjp [5], value = ZERO;

  for (xip [0] = yjp [0] = ONE, k=1; k<=4; k++)
  {
    xip [k] = xip [k-1] * xi;
    yjp [k] = yjp [k-1] * yj;
  }
  for (k=0; k<=3; k++)
    for (l=0; l<=3; l++)
      value += a [i][j][k][l] * xip [k+1] * yjp [l+1] / ((k+1)*(l+1));
  return value;
}

/* -------------------------- END fibikub.c ------------------------- */
