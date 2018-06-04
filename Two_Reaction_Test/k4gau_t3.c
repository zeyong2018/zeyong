#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test program for cubature over rectangles using  Gauss               *
*                                                                      *
* We test the function Kub4GauV                                        *
*                                                                      *
* Author         Uli Eggermann                                         *
* Date           9.2.1991                                              *
***********************************************************************/

#include <basis.h>
#include <vmblock.h>
#include <kubatur.h>

typedef enum { LINEAR, HALBQUADRATISCH, QUADRATISCH } temptyp;
char *LHQ = "LHQ";

REAL   y_welle (REAL,REAL);       /* Function declared   in F_Beisp.C */
extern int function_calls;        /* Variable deklared   in F_Beisp.C */

void Zerlegung (int  Zerlegungsart,
                REAL* x, int Nx, REAL a, REAL b,
                REAL* y, int Ny, REAL c, REAL d);

int main (void)
{
  REAL   a, b, c, d;
  REAL   W, F;
  REAL   *x, *y;
  int    Zerl;
  int    Nx = 4, Ny = 4;
  int    Verfahren;
  int    mSCH;                                    /* with estimation */
  const REAL   pi_halbe = (REAL)1.570796327;
  void   *vmblock;


  vmblock = vminit();
  x = (REAL *)vmalloc(vmblock, VEKTOR, Nx + 1, 0);
  y = (REAL *)vmalloc(vmblock, VEKTOR, Ny + 1, 0);
  if (! vmcomplete(vmblock))
    return 1;

  a = -pi_halbe; b = -a;
  c = a;         d = -c;

  for (mSCH = 0; mSCH <= 1; mSCH++)
   for (Verfahren = 0; Verfahren <= 7; Verfahren++)
   {
    for (Zerl = LINEAR; Zerl <= QUADRATISCH; Zerl++)
    {
      int i;
      function_calls = 0;         /* Variable declared   in F_Beisp.C */
      Zerlegung (Zerl, x,Nx,a,b, y,Ny,c,d);
      i = Kub4GauV (x, Nx, y, Ny, Verfahren,
                    y_welle,
                    &W, mSCH, &F);
      if (i == 0) {
        printf ("%d. %c (%2d) value=%12.8"LZP"f",
                Verfahren, LHQ [Zerl], i, W);
        if (mSCH) printf (" error estimate =%10.3"LZP"e", F);
        else      printf ("%27c",' ');
        printf ("%5d Function evaluations\n", function_calls);
      }
      else printf ("method  %d: error %d\n", Verfahren, i);
    }
    printf ("\n");
   }

  return 0;
}

void Zerlegung (int Zerlegungsart,
                REAL* Ax, int Nx, REAL a, REAL b,
                REAL* Ay, int Ny, REAL c, REAL d)
{
  int i, j;
  #define x ((REAL  ) i / (REAL  ) Nx)
  #define y ((REAL  ) j / (REAL  ) Ny)
  switch (Zerlegungsart) {
  case LINEAR:
    for (i = 0; i <= Nx; i++)  Ax[i] = a + (b - a) * x;
    for (j = 0; j <= Ny; j++)  Ay[j] = c + (d - c) * y;
    break;
  case HALBQUADRATISCH:
    for (i = 0; i <= Nx; i++)  Ax[i] = a + (b - a) * (x * x);
    for (j = 0; j <= Ny; j++)  Ay[j] = c + (d - c) * (y * y);
    break;
  case QUADRATISCH:
    #define h      HALF
    #define g(k,N) (k <= (N/2) ? -FOUR : FOUR)
    for (i = 0; i <= Nx; i++)
      Ax[i] = h * (a+b) + h * (b-a) * (x-h)*(x-h) * g(i,Nx);
    for (j = 0; j <= Ny; j++)
      Ay[j] = h * (c+d) + h * (d-c) * (y-h)*(y-h) * g(j,Ny);
    break;
  }
}
