#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, fprintf,   */
                         /*       stderr, NULL, REAL, fehler_melden,  */
                         /*       LZS, LZP, ZERO                      */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <hrmsplin.h>    /*  for  hermit                              */
#include <spliwert.h>    /*  for  hmtwert                             */
#include <splintab.h>    /*  for  hmtab                               */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define TABELLENLAENGE  100



/***********************************************************************
*                                                                      *
* Testprogramm: non periodic interpolating Hermite splines of degree 5:*
*               Find coefficients and table of values,                 *
*               plot if Turbo C compiler                               *
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
* x1  y1 y1'        (1st node: x-value, y-value, value of derivative)  *
* x2  y2 y2'        (2nd node: ditto)                                  *
* ...               ...                                                *
* xn  yn yn'        (nth node: ditto)                                  *
* mass              (type of end point data)                           *
* alfa beta         (left end data, right end data)                    *
*                                                                      *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL *x, *y,
       *yi,
       *c, *d, *e, *f,
       *xtab, *ytab,
       alfa, beta,
       xmit,
       hilf,
       ausg[5],
       laenge;           /* Length of interpolation interval          */
  int  n,
       i,
       fehler,
       mass,
       lentab;
  void *vmblock;         /* List of dynamically allocated vectors     */


  if ((fehler = umleiten(argc, argv))   /* assign input/output file   */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nnumber of nodes:                     ");
#endif
  scanf("%d", &n);

  vmblock = vminit();                 /* initialize storage           */
  x    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  y    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  yi   = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  c    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  d    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  e    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  f    = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
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
    fprintf(stderr, "x-, y-coordinates, 1st derivative: ");
#endif
    scanf("%"LZS"f%"LZS"f%"LZS"f", x + i, y + i, yi + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "type of end point condition (1...5): ");
#endif
  scanf("%d", &mass);

  if (mass < 6 && mass > 2)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "boundary conditions:               ");
#endif
    scanf("%"LZS"f%"LZS"f", &alfa, &beta);
  }
  else
    alfa = beta = ZERO;                      /* for   IBM C Set/2 1.0 */


  /* ------------ print out input data ------------------------------ */

  printf("\n"
         "Compute and tabulate an interpolating "
         "Hermite spline\n"
         "======================================"
         "==============\n"
        );

  printf("\n\n"
         "number of nodes:              %d (= n+1)\n"
         "kind of end point conditions: %d\n",
         n, mass
        );

  if (mass < 6 && mass > 2)
  {
    printf("end point conditions:         alfa = %20.13"LZP"f\n"
           "                              beta = %20.13"LZP"f\n",
           alfa, beta
          );
  }

  printf("\n\n"
         "nodes:\n\n"
         " i          x[i]                 y[i]                 y'[i]\n"
         "-----------------------------------------------------------"
         "------\n"
        );
  for (i = 0; i < n; i++)
    printf("%2d %20.13"LZP"f %20.13"LZP"f %20.13"LZP"f\n",
           i, x[i], y[i], yi[i]);


                           /* call  hermit() to compute coefficients  */

  fehler = hermit(n, x, y, yi, mass, alfa, beta, 0, c, d, e, f);

  if (fehler != 0)
  {
    fehler_melden("hermit()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }


  /* ---------- print output ---------------------------------------- */

  printf("\n\n"
         "Spline coefficients:\n\n"
          " i        a[i]               b[i]               "
          "c[i]\n"
          "------------------------------------------------"
          "----------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f\n",
           i, y[i], yi[i], c[i]);

  printf("\n"
          " i        d[i]               e[i]               "
          "f[i]\n"
          "------------------------------------------------"
          "----------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f\n",
           i, d[i], e[i], f[i]);


  /* --- put out function and derivative values for 5 points ------- */

  printf("\n\n"
         "a sample of function and derivative values of spline:\n\n"
         "     x0         y          "
         "y'         y''        y'''       y''''     y'''''\n"
         "--------------------------------------"
         "--------------------------------------\n"
        );
  xmit = (x[n - 1] - x[0]) / (REAL)7.0;
  for (i = 1; i <= 5; i++)
  {
    hilf = hmtwert(n - 1, x[0] + i * xmit, y, yi, c, d, e, f, x,
                   ausg);
    printf("%10.5"LZP"f %10.5"LZP"f %10.5"LZP"f"
           " %10.5"LZP"f %10.5"LZP"f %10.5"LZP"f"
           " %10.5"LZP"f\n", x[0] + i * xmit,
           hilf, ausg[0], ausg[1], ausg[2], ausg[3], ausg[4]
          );
  }


  /* ----- as a control compute values at the nodes ----------------- */

  printf("\n\n"
         "Function and derivative values for the spline "
         "at the nodes:\n\n"
         "     x[i]       y          "
         "y'         y''        y'''       y''''     y'''''\n"
         "--------------------------------------"
         "--------------------------------------\n"
        );
  for (i = 0; i < n; i++)
  {
    hilf = hmtwert(n - 1, x[i], y, yi, c, d, e, f, x, ausg);
    printf("%10.5"LZP"f %10.5"LZP"f %10.5"LZP"f"
           " %10.5"LZP"f %10.5"LZP"f %10.5"LZP"f"
           " %10.5"LZP"f\n",
           x[i], hilf, ausg[0], ausg[1], ausg[2], ausg[3], ausg[4]
          );
  }


  /* tabulate the spline from  x[0] to  x[n-1] with minimal step size,*/
  /* i.e., so that the table is optimally filled -------------------- */

  laenge = x[n - 1] - x[0];
  alfa = x[0] -     (REAL)0.015 * laenge;
  beta = x[n - 1] + (REAL)0.015 * laenge;
  laenge = beta - alfa;
  fehler = hmtab(n - 1, alfa, beta, laenge / (TABELLENLAENGE - 1 - n),
                 TABELLENLAENGE - 1, x, y, yi, c, d, e, f, xtab, ytab,
                 &lentab);
  if (fehler != 0)
  {
    fehler_melden("hmtab()", 20 + fehler, __FILE__, __LINE__);
    return 20 + fehler;
  }
  printf("\n\n"
         "tabulate spline from about x[0] to x[n] with "
         "minimal step size:\n\n"
         " i       xtab[i]              ytab[i]\n"
         "--------------------------------------------\n"
        );
  for (i = 0; i <= lentab; i++)
    printf("%2d %20.14"LZP"f %20.14"LZP"f\n",
           i, xtab[i], ytab[i]);


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
