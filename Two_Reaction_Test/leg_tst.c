#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test program for the  Legendre functions of module STGEW.C:          *
*                                                                      *
*    LegendreCoeff  to fill in coefficient vector                      *
*    LegPolWert     to evaluate function values                        *
*                                                                      *
*                                                                      *
* We compute the coefficients of  Legendre polynomials of degrees 2 to *
* 10 and plot function values in the interval [0,1].                   *
*                                                                      *
* Uli Eggermann, 9.28.1991                                             *
***********************************************************************/

#include <basis.h>
#include <gax.h>

int main (void)
{
  #define MaxGrad 10

  int Grad;
  REAL x, C [MaxGrad+1][MaxGrad+1];

  for (Grad = 1; Grad <= MaxGrad; Grad++)
    LegendreCoeff (C [Grad], Grad, NULL);

  printf("Functional values of the Legendre polynomial of degree n:"
         " Pn (x)\n");

  printf ("  x  |");
  for (Grad = 2; Grad <= MaxGrad; Grad++)
    printf (" P%d(x)%s|", Grad, Grad < 10 ? " " : "");
  printf ("\n");
  printf ("-----+");
  for (Grad = 2; Grad <= MaxGrad; Grad++)
    printf ("-------%c", Grad < MaxGrad ? '+' : '|');
  printf ("\n");

  for (x = 0.0; x < 1.01; x+=0.05) {
    printf ("%.2"LZP"f |", x);
    for (Grad = 2; Grad <= MaxGrad; Grad++) {
      printf ("%+7.4"LZP"f|", LegPolWert (x, C[Grad], Grad));
    }
    printf ("\n");
  }

  return 0;
}
