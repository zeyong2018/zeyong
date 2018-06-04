#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>     /*  for  REAL, freopen, stdout, NULL, fprintf, */
                       /*       stderr, fehler_melden, PI, SIN, ZERO, */
                       /*       COS, ONE, TWO, max, FABS, printf, LZP */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vminit, VVEKTOR  */
#include <fft.h>       /*  for  complex, fourn                        */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )


/***********************************************************************
*         This program tests the function fourn() with the values of   *
*         the function  F(X)  defined as                               *
*               0.5          for  X = 0                                *
*              (PI-X)/PI     for  X  in (0, PI]                        *
*               0            otherwise,                                *
*         whose (nonperiodic) Fourier transform  F^(V)  is given by    *
*                                                                      *
*         (1-COS(2*PI*PI*V)+I*(SIN(2*PI*PI*V)-2*PI*PI*V)/(4*PI**3*V*V).*
*                                                                      *
*         Here I denotes the imaginary unit (I^2 = -1).                *
*         A second test function is  G(X)  defined as                  *
*              SIN(X)        for  X  in [0, PI]                        *
*               0            otherwise,                                *
*         whose (nonperiodic) Fourier transform G^(V)  is              *
*                                                                      *
*            (COS(2*PI*PI*V)+1-I*SIN(2*PI*PI*V))/(1-4*PI*PI*V*V).      *
*                                                                      *
*         In both cases we use equidistant nodes. For increasing       *
*         numbers of nodes, the accuracy must increase; specifically   *
*         it must do so more rapidly for G(X) than for F(X).           *
***********************************************************************/

{
#define ANZAHLBEISPIELE     7
#define maxM             1366

#define fre(v)  (ONE-COS(TWO*pi2*(v)))/(FOUR*pi2*PI*(v)*(v))
#define fim(v)  (SIN(TWO*pi2*(v))-TWO*pi2*(v))/(FOUR*pi2*PI*(v)*(v))
#define gre(v)  (COS(TWO*pi2*(v))+ONE)/(ONE-FOUR*pi2*(v)*(v))
#define gim(v)  -SIN(TWO*pi2*(v))/(ONE-FOUR*pi2*(v)*(v))

  static int M[ANZAHLBEISPIELE] = { 16, 64, 288, 512, 798, 1024, maxM };

  complex    *F;        /* [0..maxM-1] vector of f function values    */
  complex    *G;        /* [0..maxM-1] vector with g function values  */
  complex    fdiff;     /* distance betwen the exact abd the computed */
                        /* functional value                           */
  complex    gdiff;     /* ditto for g                                */
  REAL       Pi;        /* 3.14...                                    */
  REAL       pi2;       /* `Pi*Pi'                                    */
  REAL       fakt;      /* step size                                  */
  REAL       xj;        /* current node for f and g                   */
  REAL       tk;        /* node of Fourier transform                  */
  int        fehler;    /* error code from fourn()                    */
  int        i;         /* inter variable                             */
  int        j;         /* ditto                                      */
  int        k;         /* ditto                                      */
  void       *vmblock;  /* List of dynamically allocated vectors      */


  /* --- assign output file to standard output                     -- */

  if (argc >= 2)                           /* at least one argument?  */
    if (freopen(argv[1], "w", stdout) == NULL)        /* open output  */
    {                                                 /* file      */
      fprintf(stderr, "Error in opening %s!\n", argv[1]);
      return 1;
    }


  vmblock = vminit();                 /* initialize storage           */
  F = (complex *)vmalloc(vmblock, VVEKTOR, maxM, sizeof(*F));
  G = (complex *)vmalloc(vmblock, VVEKTOR, maxM, sizeof(*G));
  if (! vmcomplete(vmblock))   /* Storage demand not                  */
  {                            /*          satisfied?                 */
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }


  Pi  = PI;
  pi2 = sqr(Pi);

    /* find real and imaginary parts of functional values */

  for (i = 0; i < ANZAHLBEISPIELE; i++)
  {
    fakt = PI / (REAL)M[i];
    for (j = 0; j < M[i]; j++)
    {
      xj     = (REAL)j * fakt;
      F[j].x = (PI -  xj) / PI;
      F[j].y = ZERO;
      G[j].x = SIN(xj);
      G[j].y = ZERO;
    }
    F[0].x = HALF;

    /* Find the Fourier transform               */

    fehler = fourn(M[i], F, ZERO, fakt)  |
             fourn(M[i], G, ZERO, fakt);

    switch (fehler)
    {
      case 0:
        break;
      case 1:
        fehler_melden("fourn(): M < 1", 10 + fehler,
                      __FILE__, __LINE__);
        break;
      case 2:
        fehler_melden("fourn(): M too large",
                      10 + fehler, __FILE__, __LINE__);
        break;
      case 3:
        fehler_melden("fourn(): lack of memory",
                      10 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("fourn(): unknown error",
                      10 + fehler, __FILE__, __LINE__);
    }

    /* Compare with exact Fourier transformn               */

    fdiff.x = fdiff.y = gdiff.x = gdiff.y = ZERO;
    for (k = -M[i] / 2; k < 0; k++)
    {
      tk      = (REAL)k / Pi;
      fdiff.x = max(fdiff.x, FABS(fre(tk) - F[k + M[i]].x));
      fdiff.y = max(fdiff.y, FABS(fim(tk) - F[k + M[i]].y));
      gdiff.x = max(gdiff.x, FABS(gre(tk) - G[k + M[i]].x));
      gdiff.y = max(gdiff.y, FABS(gim(tk) - G[k + M[i]].y));
    }
    tk      = ZERO;
    fdiff.x = max(fdiff.x, FABS(HALF * PI - F[0].x));
    fdiff.y = max(fdiff.y, FABS(F[0].y));
    gdiff.x = max(gdiff.x, FABS(gre(tk) - G[0].x));
    gdiff.y = max(gdiff.y, FABS(gim(tk) - G[0].y));
    for (k = 1; k < M[i] / 2; k++)
    {
      tk      = (REAL)k / Pi;
      fdiff.x = max(fdiff.x, FABS(fre(tk) - F[k].x));
      fdiff.y = max(fdiff.y, FABS(fim(tk) - F[k].y));
      gdiff.x = max(gdiff.x, FABS(gre(tk) - G[k].x));
      gdiff.y = max(gdiff.y, FABS(gim(tk) - G[k].y));
    }
    printf("M = %4d  error with regard to f(x):  %10.3"LZP"e"
           ",%10.3"LZP"e\n", M[i], fdiff.x, fdiff.y);
    printf("M = %4d  error with regard to g(x):  %10.3"LZP"e"
           ",%10.3"LZP"e\n", M[i], gdiff.x, gdiff.y);
    printf("\n");
  }


  return 0;
}
