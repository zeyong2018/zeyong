#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test programm for cubature over triangles via summed Gaussian n-point*
* formula                                                              *
*                                                                      *
* We test the funktion Kub3GauN                                        *
*                                                                      *
* Author          Uli Eggermann                                        *
* Date            8.23.1991                                            *
***********************************************************************/

#include <basis.h>
#include <kubatur.h>

int main (void) {
  int    i, Verfahren, Kantenteilung;
  REAL   WertDerKubatur;
  REAL   Px, Py, Qx, Qy, Rx, Ry;
  REAL   exp2 (REAL, REAL);                /* Function from F_Beisp.C */

  Px = ZERO;  Py = ZERO;
  Qx =  TEN;  Qy = ZERO;
  Rx = ZERO;  Ry =  TEN;


  for (Verfahren = 1; Verfahren <= 7; Verfahren++)
    if (Verfahren < 4 || Verfahren > 6)
    {
      printf (" Method %d:\n", Verfahren);
      for (Kantenteilung = 1; Kantenteilung <= 5; Kantenteilung++)
      {
        printf ("%2d Subtriangle%s", Kantenteilung * Kantenteilung,
                Kantenteilung * Kantenteilung > 1 ? "s:" : ": ");

        i = Kub3GauN (Px, Py, Qx, Qy, Rx, Ry, Verfahren,
                      Kantenteilung,
                      exp2,                /* Function from F_Beisp.C */
                      &WertDerKubatur);

        printf (" value =%18.10"LZP"f (Error=%d)\n", WertDerKubatur, i);

      } /* for Kantenteilung */
    } /* if Verfahren */

  return 0;
} /* main */
