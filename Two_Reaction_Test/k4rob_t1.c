#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/**********************************************************************
* Test program for cubature over rectangles using Richardson          *
* extrapolation and                                                   *
*           a) the Romberg sequence  1, 2, 4, 8, 16, 32, 64 ..., and  *
*           a) the Bulirsch sequence 1, 2, 3, 4, 6, 8, 12, 16 ...   . *
*                                                                     *
* We test the functions Kub4BuRi and Kub4RoRi                         *
*                                                                     *
* Author         Uli Eggermann                                        *
* Date           8.23.1991                                            *
**********************************************************************/

#include <basis.h>
#include <kubatur.h>

REAL   exp2 (REAL,REAL);         /* Function declared   in F_Beisp.C */

int main (void)
{
  REAL   a = ZERO, b = TEN, c = ZERO, d = TEN, W, F;
  int    i, Nx = 4, Ny = 3, NN;

  for (NN = 3; NN < 6; NN ++)
  {
    printf ("\n");

    W = F = ZERO;
    i = Kub4BuRi (a, b, Nx, c, d, Ny, NN,
                  exp2,
                  &W, &F);
    printf ("Bulirsch: %2d: (%d) W = %.10"LZP"e F = %.3"LZP"e\n",
            NN, i, W, F);

    W = F = 0.0;
    i = Kub4RoRi (a, b, Nx, c, d, Ny, NN,
                  exp2,
                  &W, &F);
    printf ("Romberg : %2d: (%d) W = %.10"LZP"e F = %.3"LZP"e\n",
            NN, i, W, F);
  }

  return 0;
}
