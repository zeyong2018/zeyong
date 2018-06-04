#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test-Programm for GAX                                 Egg, 9.29.1991 *
***********************************************************************/

#include <basis.h>
#include <vmblock.h>
#include   <gax.h>

REAL quadratic_hyperbola (REAL x)
{
  return 1.0 / (x* x);
}

REAL bell_curve (REAL x)
{
  return EXP (- x* x);
}

REAL test_function (REAL x)
{
  REAL w = x - 0.5;
  return -100 * w * EXP (-50 * w * w) - 2 * EXP (-2 * x);
}

void prinfos (
              int       Methode,
              IntInfo** Knoten,
              int       nMax,
              REAL      Qwert,
              REAL      AbsFehl,
              int       GXFehl
             )
{
  int i;
  printf ("\n n Error est. Quadrature est.       left  right\n");

  for (i = 0; i < nMax; i++)
  {
      printf ("%2d: %10.3"LZP"e %+-20.7"LZP"g %-6.4"LZP"g "
              "%-6.4"LZP"g\n", i,
              Knoten[i]->FSch,  Knoten[i]->QRes,
              Knoten[i]->Links, Knoten[i]->Rechts  );
  }
  printf("Method = %d: Integral = %-10.5"LZP"g, Error :A=%.3"LZP"e "
         "R=%.3"LZP"e (%d)\n",
          Methode, Qwert, AbsFehl,
          FABS (Qwert) < MACH_EPS ? ZERO : AbsFehl / Qwert, GXFehl);
}

int main (void)
{
  #define Max  40
  IntInfo **Knoten;                   /* Information for subintervals */
  REAL    Qwert;                              /* approximate integral */
  REAL    AbsFehl;                         /* Absolute error of Qwert */
  int     GXFehl;                     /* error return from  GAX: 0=ok */

  /*                                    Prepare intervals from A to B */

  #define  percent       * 0.01
  #define  TestFunction  4

  #if TestFunction == 0
    REAL Zerlegung[] = {1,2,4,6,8,10};         /* improper partition  */

  #elif TestFunction == 1
    REAL Zerlegung[] = { 1, 1.5, 2, 3, 5, 10 };
    REAL (*f) (REAL) = quadratic_hyperbola;
    REAL RelErr = 1 percent;

  #elif TestFunction == 2
    REAL Zerlegung[] = { -1, -0.3, -0.1, 0.1, 0.3, 1 };
    REAL (*f) (REAL) = bell_curve;
    REAL RelErr = 1 percent;

  #elif TestFunction == 3
    REAL Zerlegung[] = { 0, 0.1, 0.2, 0.5 };
    REAL (*f) (REAL) = bell_curve;
    REAL RelErr = 1 percent;

  #elif TestFunction == 4
    REAL Zerlegung[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
    REAL (*f) (REAL) = test_function;
    REAL RelErr = 0.01 percent;

  #elif TestFunction == 5
    REAL Zerlegung[] = { 0.0, 1.0 };
    REAL (*f) (REAL) = test_function;
    REAL RelErr = 1E-4;
  #endif
  void    *vmblock = NULL;

  int Nanf = (int)(sizeof(Zerlegung)/
                   sizeof(Zerlegung[0]))-1; /* number of subintervals */

  int nMax;                /* maximal allowed number of intervals     */
  int Methode;             /* Quadrature method                       */
  int MethPar[5];          /* Secondary parameter for method          */

  char *MethStr [] = { "TRAPEZ", "GAUSSJORDAN", "CLENSHAWCURTIS",
                       "ROMBERG", "NEWTONCOTES" };

  MethPar [GxTRAPEZ        ] =  0;
  MethPar [GxGAUSSJORDAN   ] =  7;
  MethPar [GxCLENSHAWCURTIS] =  2;
  MethPar [GxROMBERG       ] =  0;
  MethPar [GxNEWTONCOTES   ] =  7;

  printf ("Test program for Gax: Integral depending on method\n"
          "with error estimate controlled interval halving\n\n"
          "Test function is "
          "f(x) = -100 * w * exp (-50 * w * w) - 2 * exp (-2 * x)"
          "\ninitial partition = "
          "{ 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 }"
          "\n"
          );

  for (Methode = GxTRAPEZ; Methode <= GxNEWTONCOTES; Methode++)
  {
    Qwert  = AbsFehl = ZERO;
    nMax   = Max;

    vmblock = vminit();

    printf ("\nMethod %s:\n", MethStr [Methode]);

    GXFehl = GaxT (Zerlegung, Nanf,     RelErr,  f,
                   &nMax,     Methode,  MethPar [Methode],
                   &Qwert,    &AbsFehl, &Knoten, vmblock);

    prinfos (Methode, Knoten, nMax, Qwert, AbsFehl, GXFehl);

    GaxF (vmblock);                             /* allocated in  GaxT */
    vmblock = NULL;

  }

  return 0;
}
