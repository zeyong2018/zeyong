#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*****************************************************************/
/* Test program for cubature over rectangles using  Gauss        */
/*                                                               */
/* We test the function Kub4GauE                                 */
/*                                                               */
/* Author          Uli Eggermann                                 */
/* Date            8.30.1991                                     */
/*****************************************************************/

#include <basis.h>
#include <kubatur.h>

int main (void)
{
  REAL   a, b, c, d;
  REAL   W, F;
  int    Nx = 4, Ny = 4;
  int    Verfahren;
  int    mSCH;                                  /* with eror estimate */

  REAL   exp2 (REAL,REAL);        /* Function declared in F_Beisp.C   */
  extern int function_calls;      /* Variable declared in F_Beisp.C   */

  a = c = ZERO;
  b = d = TEN;

  for (mSCH = 0; mSCH <= 1; mSCH++)
    for (Verfahren = 0; Verfahren <= 7; Verfahren++)
    {
      int i;

      function_calls = 0;

      i = Kub4GauE (a, b, Nx, c, d, Ny, Verfahren,
                    exp2,                           /* from F_Beisp.C */
                    &W, mSCH, &F);

      printf ("%d. (%2d) value=%12.8"LZP"f",   Verfahren, i, W);

      if (mSCH) printf (" error estimate =%10.3"LZP"e", F);
      else      printf ("%27c",' ');

      printf (" %5d function evaluations\n", function_calls);
    }

  return 0;
}
