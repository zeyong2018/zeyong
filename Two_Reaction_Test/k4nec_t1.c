#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*****************************************************************/
/* Test program for cubature over rectangles using Newton-Cotes  */
/*                                                               */
/* We test the function Kub4NeCn                                 */
/*                                                               */
/*****************************************************************/
/*****************************************************************/
/* Author         Uli Eggermann                                  */
/*****************************************************************/

#include <basis.h>
#include <kubatur.h>

int main (void)
{
  REAL   a, b, c, d;
  REAL   W, F;
  int    Nx = 10, Ny = 10;
  int    Verfahren;
  int    mSCH;                                    /* with estimation  */
  REAL   exp2 (REAL,REAL);        /* Function declared   in F_Beisp.C */

  a = c = ZERO;
  b = d = TEN;

  for (mSCH = 0; mSCH <= 1; mSCH++)
    for (Verfahren = 1; Verfahren <= 7; Verfahren++)
    {
      int i, nF;

      i = Kub4NeCn (a, b, Nx, c, d, Ny, Verfahren,
                    exp2,
                    &W, mSCH, &F);

      printf ("%d. (%2d) value=%13.10"LZP"f",   Verfahren, i, W);

      if (mSCH) printf (" error estimate =%10.3"LZP"e", F);
      else      printf ("%27c",' ');

      nF = (Nx * Verfahren + 1) * (Ny * Verfahren + 1);

      if (mSCH) nF += (2 * Nx * Verfahren + 1)*(2 * Ny * Verfahren + 1);

      printf (" %5d function evaluations .\n", nF);
    }

  return 0;
}
