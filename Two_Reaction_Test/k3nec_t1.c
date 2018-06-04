#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test program for cubature over triangles using 3-point Newton-Cotes  *
* and Romberg-Richardson extrapolation                                 *
*                                                                      *
* We test the functions Kub3RoRi and Kub3NeC3                          *
*                                                                      *
* Author          Uli Eggermann                                        *
* Date            8.23.1991                                            *
***********************************************************************/

#include <basis.h>
#include <kubatur.h>

int main (void)
{
  REAL   W, F;
  int    i, j;
  REAL   xmaly (REAL, REAL);                        /* from F_Beisp.C */

  for (j = 1; j <= 7; j++) {

    /*            Px,Py, Qx,Qy, Rx,Ry */
    i = Kub3RoRi ( ZERO, ZERO, TEN, ZERO, ZERO, TEN, j,
                   xmaly,                           /* from F_Beisp.C */
                   &W, &F);

    printf ("%2d Newton-Cotes: Error=%2d, ", j, i);
    printf (" value = %12.8"LZP"f, error = %.5"LZP"e\n", W, F);
  }

  return 0;
}
