#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>     /*  for  REAL, freopen, stdout, NULL, fprintf, */
                       /*       stderr, fehler_melden, PI, FABS, SIN, */
                       /*       ZERO, HALF, COS, printf, LZP          */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vminit, VVEKTOR  */
#include <fft.h>       /*  for  complex, ffako                        */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
*         This program performs a test of the function ffako()  using  *
*         values for the PI-periodic functions                         *
*              | SIN X |                                               *
*         and                                                          *
*               0.5         for  X = 0                                 *
*               (PI-X)/PI   for  0 < X < PI  ,                         *
*         whose cyclic convolution is                                  *
*              2*(PI-X)/PI**2 - COS(X)/PI ,                            *
*         and whose cyclic correlation is                              *
*              2*X/PI**2 + COS(X)/PI .                                 *
*         We use equidistant nodes in the period interval [0, PI).     *
*         Note that the matching accuracy must improve for increasing  *
*         numbers of nodes.                                            *
***********************************************************************/

{
#define ANZAHLBEISPIELE    5
#define maxM             512

  static int M[ANZAHLBEISPIELE] = { 16, 41, 64, 188, maxM };

  complex    *F;        /* [0..maxM-1] vector with values for test    */
                        /* function (1)                               */
  complex    *H;        /* [0..maxM-1] vector with values for test    */
                        /* function (2)                               */
  complex    ffalt;
  complex    fkorr;
  REAL       fufal;
  REAL       fukor;
  REAL       fakt;
  REAL       xj;
  int        fehler;    /* error code from  fdicht()                  */
  int        i;         /* variable                                   */
  int        j;         /* variable                                   */
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

    /* find real and imaginary parts of functional values */

  for (i = 0; i < ANZAHLBEISPIELE; i++)
  {
    fakt = PI / (REAL)M[i];
    for (j = 0; j < M[i]; j++)
    {
      xj     = (REAL)j * fakt;
      F[j].x = FABS(SIN(xj));
      F[j].y = ZERO;
      H[j].x = (PI - xj) / PI;
      H[j].y = ZERO;
    }
    H[0].x = HALF;

    /* ompute discrete cyclic convolution and correlation      */

    fehler = ffako(M[i], F, H);

    switch (fehler)
    {
      case 0:
        break;
      case 1:
        fehler_melden("ffako(): M < 2", 10 + fehler,
                      __FILE__, __LINE__);
        break;
      case 2:
        fehler_melden("ffako(): M too large",
                      10 + fehler, __FILE__, __LINE__);
      case 3:
        fehler_melden("ffako(): memory exhausted",
                      10 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("ffako(): unknown error",
                      10 + fehler, __FILE__, __LINE__);
    }

    /* Compare with exact values                                   */

    ffalt.x = ffalt.y = fkorr.x = fkorr.y = ZERO;
    for (j = 0; j < M[i]; j++)
    {
      xj      = (REAL)j * fakt;
      fufal   = TWO * (PI - xj) / (PI * PI) - COS(xj) / PI;
      fukor   = TWO * xj / (PI * PI) + COS(xj) / PI;
      ffalt.x = max(ffalt.x, FABS(F[j].x - fufal));
      ffalt.y = max(ffalt.y, FABS(F[j].y));
      fkorr.x = max(fkorr.x, FABS(H[j].x - fukor));
      fkorr.y = max(fkorr.y, FABS(H[j].y));
    }
    printf("\n"
           "M = %3d  Error in convolution:   %10.3"LZP"e,%10.3"LZP"e\n"
           "M = %3d  Error in correlation:   %10.3"LZP"e,%10.3"LZP"e\n",
           M[i], ffalt.x, ffalt.y, M[i], fkorr.x, fkorr.y);
  }


  return 0;
}
