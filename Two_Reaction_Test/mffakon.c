#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>     /*  for  REAL, freopen, stdout, NULL, fprintf, */
                       /*       stderr, fehler_melden, PI, SIN, ZERO, */
                       /*       COS, ONE, TWO, max, FABS, printf, LZP */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vminit, VVEKTOR  */
#include <fft.h>       /*  for  complex, ffakon                       */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/*****************************************************************
*  This program tests the function ffakon()   with values of the *
*  functions                                                     *
*       SIN(X)       for  X in [0, PI]            (1)            *
*        0           otherwise                                   *
*  and                                                           *
*        0.5         for  X = 0                   (2)            *
*       (PI-X)/PI    for  X in (0, PI]                           *
*        0           otherwise,                                  *
*  whose (nonperiodic) convolution is given by                   *
*       1 - X/PI - COS(X) + SIN(X)/PI   for  X  in [0, PI]       *
*       2 - X/PI + SIN(X)/PI            for  X  in (PI, 2*PI]    *
*       0                               otherwise,               *
*  while their (nonperiodic) correlation is                      *
*       1 + X/PI + SIN(X)/PI            for  X  in [-PI, 0]      *
*       X/PI + COS(X) + SIN(X)/PI       for  X  in (0, PI]       *
*       0                               otherwise.               *
*  We use equidistant nodes. An increase in the number of nodes  *
*  must increase the accuracy of the approximation.              *
*****************************************************************/

{
#define ANZAHLBEISPIELE     8
#define maxM             1024

  static int M[ANZAHLBEISPIELE]   ={10, 10, 41, 41, 128, 128, 500, 500};
  static int tau[ANZAHLBEISPIELE] ={ 5,  5,  7,  7,   9,   9,  10,  10};

  int        faltung;   /* Flag for convolution or correlation        */
  int        fehler;    /* error code of  ffakon()                    */
  int        i;         /* variable                                   */
  int        j;         /* variable                                   */
  complex    *F;        /* [0..maxM-1] vector with values for test    */
                        /* function (1)                               */
  complex    *H;        /* [0..maxM-1] vector with values for test    */
                        /* function (2)                               */
  complex    ffalt;
  complex    fkorr;
  REAL       xj;
  REAL       dx;
  REAL       fufal;
  REAL       fukor;
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
  H = (complex *)vmalloc(vmblock, VVEKTOR, maxM, sizeof(*H));
  if (! vmcomplete(vmblock))   /* Storage demand not                  */
  {                            /*          satisfied?                 */
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

  faltung = 0;

    /* find real and imaginary parts of functional values */

  for (i = 0; i < ANZAHLBEISPIELE; i++)
  {
    dx = PI / (REAL)M[i];
    for (j = 0; j <= M[i]; j++)
    {
      xj     = (REAL)j * dx;
      F[j].x = SIN(xj);
      F[j].y = ZERO;
      H[j].x = (PI - xj) / PI;
      H[j].y = ZERO;
    }
    H[0].x = HALF;

    /* ompute convolution (faltung=1) or correlation (faltung=0)     */

    faltung = ! faltung;

    fehler = ffakon(M[i], F, M[i], H, tau[i], dx, faltung);

    switch (fehler)
    {
      case 0:
        break;
      case 1:
        fehler_melden("ffakon(): M < 1  or N < 1", 10 + fehler,
                      __FILE__, __LINE__);
        break;
      case 2:
        fehler_melden("ffakon(): tau < 1  or  tau too large",
                      10 + fehler, __FILE__, __LINE__);
      default:
        fehler_melden("ffakon(): unknown error",
                      10 + fehler, __FILE__, __LINE__);
    }

    /* Compare to exact values                        */

    ffalt.x = ffalt.y = fkorr.x = fkorr.y = ZERO;
    if (faltung)
    {
      for (j = 0; j <= 2 * M[i]; j++)
      {
        xj = (REAL)j * dx;
        if (j <= M[i])
          fufal = ONE - xj / PI - COS(xj) + SIN(xj) / PI;
        else
          fufal = TWO - xj / PI + SIN(xj) / PI;
        ffalt.x = max(ffalt.x, FABS(F[j].x - fufal));
        ffalt.y = max(ffalt.y, FABS(F[j].y));
      }
      printf("M = N = %3d  Error of convolution:    "
             "  %10.3"LZP"e,%10.3"LZP"e\n",
             M[i], ffalt.x, ffalt.y);
    }
    else
    {
      for (j = 0; j <= M[i]; j++)
      {
        xj = - PI + (REAL)j * dx;
        if (j <= M[i])
          fukor = ONE + xj / PI + SIN(xj) / PI;
        else
          fukor = xj / PI + COS(xj) + SIN(xj) / PI;
        fkorr.x = max(fkorr.x, FABS(F[j].x - fukor));
        fkorr.y = max(fkorr.y, FABS(F[j].y));
      }
      printf("M = N = %3d  Error for correlation:"
             "  %10.3"LZP"e,%10.3"LZP"e\n",
             M[i], fkorr.x, fkorr.y);
    }
  }


  return 0;
}
