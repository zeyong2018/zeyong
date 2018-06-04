#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, NULL, LZS, */
                         /*       fprintf, stderr, REAL, ZERO, LZP,   */
                         /*       fehler_melden                       */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <kubsplin.h>    /*  for  spltrans                            */
#include <spliwert.h>    /*  for  strwert                             */
#include <splintab.h>    /*  for  strtab                              */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define TABELLENLAENGE  100



/***********************************************************************
*                                                                      *
* Test program: interpolating transformed parametric polynomial        *
*               splines of degree 3                                    *
*               Compute coefficients, tabulate and plot for Turbo-C    *
*               compilers                                              *
*                                                                      *
* Input is read from the first file in the command line.               *
* If none is given input comes from the keyboard.                      *
* Output other than an error code appears in the second file of the    *
* command line.                                                        *
* For calls with less than 2 files in the command line, all output     *
* appears on the screen.                                               *
* For Turbo C compiler the spline and its nodes are plotted.           *
* This can be voided by adding a third entry of "n" on the command line*
*                                                                      *
* Construction of input files:                                         *
* n                 (number of nodes)                                  *
* x1  y1            (1st node: x-value, y-value)                       *
* x2  y2            (2nd node: ditto)                                  *
* ...               ...                                                *
* xn  yn            (nth node: ditto)                                  *
* mv                > 0: coordinates of the transformation vector given*
*                   = 0: no tranformation (px = py = 0)                *
*                   < 0: the transformation vector shall be computed   *
*                        inside  spltrans()                            *
* px  py            (transformation coordinates) (mv > 0)              *
*                                                                      *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL *x, *y,
       *phin,
       *a, *b, *c, *d,
       *xtab, *ytab,
       px,
       py,
       phid,
       panf,
       pend,
       laenge,           /* Length of interpolation interval          */
       phi0,
       ablei[4],
       c1,
       ckr,
       xk,
       yk;
  int  n,
       i,
       fehler,
       mv,
       lentab;
  void *vmblock;         /* List of dynamically allocated vectors     */


  if ((fehler = umleiten(argc, argv))   /* assign an input/output file*/
      != 0)                             /* to the standard file       */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "number of nodes:                 ");
#endif
  scanf("%d", &n);

  vmblock = vminit();                          /* initialize storage  */
  x    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  phin = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  a    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  b    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  c    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  d    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  xtab = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  ytab = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "node (x[%2d],y[%2d]):            ", i, i);
#endif
    scanf("%"LZS"f%"LZS"f", x + i, y + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "type of transformation (-1...1): ");
#endif
  scanf("%d", &mv);

  if (mv > 0)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "transformation vector (px,py):   ");
#endif
    scanf("%"LZS"f%"LZS"f", &px, &py);
  }
  else
    px = py = ZERO;


  /* ------------ print input data --------------------------------- */

  printf("\n"
         "Transformed parametric interpolating spline\n"
         "===========================================\n"
        );

  printf("\n\number of nodes: %d (= n+1)\n", n);

  if (mv < 0)
    printf("transformation vector to be computed in spltrans().\n");
  else if (mv == 0)
    printf("no transformation.\n");
  else
    printf("transformation vector (px, py) = (%"LZP"f, %"LZP"f))\n",
           px, py);

  printf("\n\n"
         "nodes:\n\n"
         " i           x[i]                  y[i]\n"
         "------------------------------------------------\n"
        );
  for (i = 0; i < n; i++)
    printf("%2d%23.13"LZP"g%23.13"LZP"g\n", i, x[i], y[i]);


  /* ----------------- compute spline ------------------------------- */

  /* -------- call spltrans() to compute spline coefficients  ------- */

  fehler = spltrans(n, x, y, mv, &px, &py, a, b, c, d, phin, &phid);

  switch (fehler)
  {
    case  0:
      break;
    case -1:
      fehler_melden("spltrans(): n < 3",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case -3:
      fehler_melden("spltrans(): angles phin[i] are not monotonically"
                    "increasing", 10 + fehler, __FILE__, __LINE__);
      if (mv < 0)
        fprintf(stderr, "the transformation vector from spltrans()\n"
                        "    (px, py) = (%"LZP"f, %"LZP"f)\n"
                        "cannot be used. Please choose different "
                        "point.\n", px, py);
      break;
    case -4:
      fehler_melden("spltrans(): x[n-1] != x[0] or y[n-1] != y[0]",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("spltrans()", 10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* --------- print output ----------------------------------------- */

  printf("\n\n"
         "Spline coefficients:\n\n"
          " i           phin[i]                     a[i]               "
          "        b[i]                       c[i]                     "
          "  d[i]\n"
          "-------------------------------------------------------"
          "-------------------------------------------------------"
          "----------------------------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %27.18"LZP"g%27.18"LZP"g%27.18"LZP"g%27.18"LZP"g"
           "%27.18"LZP"g\n", i, phin[i], a[i], b[i], c[i], d[i]);
  printf("%2d %27.18"LZP"g%27.18"LZP"g                           "
         "%27.18"LZP"g\n", i, phin[n - 1], a[n - 1], c[n - 1]);
  printf("\ncoordinate shift to\n"
         "         (px, py) = (%"LZP"f, %"LZP"f)\n"
         "coordinate system rotated by %"LZP"f degrees\n", px, py,
         phid * (REAL)180.0 / PI);


  /* ------ compute several points and slopes ---------------------- */

  printf("\n\nSome function and slope values for the spline:\n\n"
         "  phi         S(phi)         S'(phi)        S''(phi)       "
         "S'''(phi)      SX(phi)        SY(phi)        derivative     "
         "curvature\n"
         "-----------------------------------------------------------"
         "------------------------------------------------------------"
         "------------\n");

  for (phi0 = phin[0]; phi0 <= phin[n - 1];
       phi0 += (phin[n - 1] - phin[0]) / (REAL)15.0)
  {
    fehler = strwert(phi0, n, phin, a, b, c, d, phid, px, py,
                     ablei, &xk, &yk, &c1, &ckr);
    switch (fehler)
    {
      case 0: printf("%10.2"LZP"e %15.5"LZP"e%15.5"LZP"e"
                     "%15.5"LZP"e%15.5"LZP"e%15.5"LZP"e%15.5"LZP"e"
                     "%15.5"LZP"e%15.5"LZP"e\n",
                     phi0, ablei[0], ablei[1], ablei[2], ablei[3],
                     xk, yk, c1, ckr);
              break;
      case 1:
      case 2: printf("%8.2"LZP"e %13.5"LZP"e %13.5"LZP"e"
                     " %13.5"LZP"e %13.5"LZP"e %13.5"LZP"e "
                     "%13.5"LZP"e   ?          %13.5"LZP"e\n",
                     phi0, ablei[0], ablei[1], ablei[2], ablei[3],
                     xk, yk, ckr);
              break;
      case 3:
      case 6: printf("%8.2"LZP"e %13.5"LZP"e %13.5"LZP"e "
                     "%13.5"LZP"e %13.5"LZP"e %13.5"LZP"e "
                     "%13.5"LZP"e %12.4"LZP"e   ?\n",
                     phi0, ablei[0], ablei[1], ablei[2], ablei[3],
                     xk, yk, c1);
              break;
      case 4:
      case 5:
      case 7:
      case 8: printf("%8.2"LZP"e %13.5"LZP"e %13.5"LZP"e "
                     "%13.5"LZP"e %13.5"LZP"e %13.5"LZP"e "
                     "%13.5"LZP"e   ?            ?\n",
                     phi0, ablei[0], ablei[1], ablei[2], ablei[3],
                     xk, yk);
    }
  }


  /* tabulate spline from phin[0] to phin[n-1] with minimal step size */

  laenge = phin[n - 1] - phin[0];
  panf   = phin[0]     - (REAL)0.015 * laenge;
  pend   = phin[n - 1] + (REAL)0.015 * laenge;

  fehler = strtab(n - 1, panf, pend,
                  phin, a, b, c, d,
                  phid, px, py, x, y,
                  TABELLENLAENGE - 1, &lentab, xtab, ytab);

  if (fehler != 0)
  {
    fehler_melden("sptab()", 20 + fehler, __FILE__, __LINE__);
    return 20 + fehler;
  }

  printf("\n\n"
        "table of values for spline from about phin[0] to phin[n] with "
        "minimal step size:\n\n"
         " i       xtab[i]              ytab[i]\n"
         "--------------------------------------------\n"
        );
  for (i = 0; i <= lentab; i++)
    printf("%2d %20.14"LZP"f %20.14"LZP"f\n", i, xtab[i], ytab[i]);


  /* ------------ plot spline if desired ---------------------------- */

  if (argc <= 3 || *argv[3] != 'n')
  {
    fehler = zeigkurv(lentab + 1, n, xtab, ytab, x, y);
    switch (fehler)
    {
      case 0:
        break;
      case 3:
        fehler_melden("zeigkurv(): lack of memory",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 4:
        fehler_melden("zeigkurv(): BGI graphics error",
                      30 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("zeigkurv(): other error",
                      30 + fehler, __FILE__, __LINE__);
    }
    if (fehler != 0)
      return 30 + fehler;
  }


  return 0;
}
