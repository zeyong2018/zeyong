#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Several test functions for cubature                                  *
* Author:         Uli Eggermann, 6.24.1990                             *
***********************************************************************/
#ifndef F_BEISP_C
#define F_BEISP_C

#include <basis.h>
#include <math.h>

int function_calls = 0;

/***********************************************************************
*                           f(x,y) = 1                                 *
***********************************************************************/
REAL one (REAL x, REAL y)
{
  function_calls++;
  x = x; y = y;
  return ONE;
}

/***********************************************************************
*                                      - (x^2 + y^2)                   *
*                           f(x,y) = e                                 *
***********************************************************************/
REAL exp2 (REAL x, REAL y)
{
  function_calls++;
  return EXP (-(x*x + y*y));
}

/***********************************************************************
*                          upper half ball                             *
***********************************************************************/
REAL hkug (REAL x, REAL y)
{
  REAL f2;
  const REAL r2 = TEN * TEN;                     /* radius squared    */
  function_calls++;
  f2 = r2 - x*x - y*y;                     /* function value squared  */
  return (f2 < ZERO ? ZERO : SQRT (f2));          /* zero or sq. root */
}

/***********************************************************************
*                           f(x,y) = xy                                *
***********************************************************************/
REAL xmaly (REAL x, REAL y)
{
  function_calls++;
  return x * y;
}

/***********************************************************************
*                         f(x,y) = cos x                               *
***********************************************************************/
REAL x_welle (REAL x, REAL y) {
  y = y;                                         /* independent of  y */
  function_calls++;
  return COS (x);
}

/***********************************************************************
*                         f(x,y) = cos y                               *
***********************************************************************/
REAL y_welle (REAL x, REAL y) {
  x = x;                                         /* independent of  x */
  function_calls++;
  return COS (y);
}

#endif
