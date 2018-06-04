#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, abl_mat1, fprintf,        */
                         /*       scanf, printf, stderr, NULL, REAL,  */
                         /*       fehler_melden, LZS, LZP             */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <kubsplin.h>    /*  for  parspl                              */
#include <spliwert.h>    /*  for  pspwert                             */
#include <splintab.h>    /*  for  partab                              */
#include <zeigkurv.h>    /*  for zeigkurv                             */



#define TABELLENLAENGE  91



/***********************************************************************
*                                                                      *
* Test program: parametric interpolating polynomial splines of degree 3*
*               Compute coefficients and a table of values;            *
*               plot if  Turbo C compiler                              *
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
* jt                (1: put in parameter nodes, 0: compute them)       *
* t1                (1st parameter node) (for jt = 1 only)             *
* t2                (2nd parameter node) (ditto)                       *
* ...                                                                  *
* tn                (ditto )                                           *
* mass              (type of end point condition)                      *
* alfax alfay       (left end  point condition)  (mass != 4)           *
* betax betay       (rirht hand condition) (mass != 4)                 *
* alfa  beta        (left and right end point conditions) (mass = 4)   *
*                                                                      *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL     *x, *y,
           *b, *c, *d,
           *by, *cy, *dy,
           *t,
           *xtab, *ytab,
           alfa[2],
           beta[2],
           xmit,
           hilf1,
           hilf2,
           laenge;       /* Length of interpolation interval          */
  int      n,
           i,
           fehler,
           mass,
           jt,
           lentab;
  abl_mat1 ausp;
  void     *vmblock;     /* List dof dynamically allocated vectors    */


  if ((fehler = umleiten(argc, argv))   /* assign in/output to        */
      != 0)                             /* standard files             */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input file --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "number of nodes:                         ");
#endif
  scanf("%d", &n);

  vmblock = vminit();                            /* initialize memory */
  x      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  y      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  b      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  c      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  d      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  by     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  cy     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  dy     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  t      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  xtab   = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  ytab   = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "node (x[%2d],y[%2d]):                    ", i, i);
#endif
    scanf("%"LZS"f%"LZS"f", x + i, y + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr,
          "If parameter values are not put in, set a zero,\n"
          "otherwise take any nonzero integer: "
         );
#endif
  scanf("%d", &jt);
  if (jt != 0)
    for (i = 0; i < n; i++)
    {
#ifdef INTERAKTIV
      fprintf(stderr, "Parameter nodes t[%2d]:               ", i);
#endif
      scanf("%"LZS"f", t + i);
    }

#ifdef INTERAKTIV
  fprintf(stderr, "type of boundary condition (0..4):       ");
#endif
  scanf("%d", &mass);

  if (mass == 1 || mass == 2)   /* prescribed 1st or 2nd derivative   */
  {                             /* wrt. t ?                           */
#ifdef INTERAKTIV
    fprintf(stderr, "end point conditions\n"
                    "(alfa[0], alfa[1], beta[0], beta[1]):  ");
#endif
    scanf("%"LZS"f%"LZS"f%"LZS"f%"LZS"f",
          alfa, alfa + 1, beta, beta + 1);
  }
  else if (mass == 4)           /* prescribed end derivatives dy/dx ? */
  {
#ifdef INTERAKTIV
    fprintf(stderr, "end derivatives (alfa[0], beta[0]):   ");
#endif
    scanf("%"LZS"f%"LZS"f", alfa, beta);
  }


  /* ------------- print input -------------------------------------- */

  printf("\n"
         "Compute and tabulate an interpolating "
         "parametric cubic spline\n"
         "======================================"
         "=======================\n"
        );

  printf("\n\n"
         "number of nodes:        %d (= n+1)\n"
         "type of end conditions: %d\n",
         n, mass
        );
  if (mass != 0 &&                       /* not  not-a-node spline?   */
      mass != 3)                         /* not periodic ?            */
    if (mass != 4)                       /* no derivivatives wrt. x ? */
      printf("end conditions:          alfax = %20.13"LZP"f\n"
             "                         alfay = %20.13"LZP"f\n"
             "                         betax = %20.13"LZP"f\n"
             "                         betay = %20.13"LZP"f\n",
             alfa[0], alfa[1], beta[0], beta[1]
            );
    else
      printf("end conditions:          alfax = %20.13"LZP"f\n"
             "                         betax = %20.13"LZP"f\n",
             alfa[0], beta[0]
            );

  printf("\n\n"
         "nodes:\n\n"
         " i          x[i]                y[i]\n"
         "--------------------------------------------\n"
        );
  for (i = 0; i < n; i++)
    printf("%2d %20.13"LZP"f %20.13"LZP"f\n", i, x[i], y[i]);

  if (jt != 0)
  {
    printf("\n\n"
           "Parameter nodes:\n\n"
           " i          t[i]       \n"
           "-----------------------\n"
          );
    for (i = 0; i < n; i++)
      printf("%2d %20.13"LZP"f\n", i, t[i]);
  }


  /* call parspl() to compute the coefficients and possibly           */
  /* ( if jt = 0) to find the parameter nodes                         */

  fehler = parspl(n, x, y, mass, alfa, beta, jt, t, b, c, d,
                  by, cy, dy);

  if (fehler != 0)
  {
    fehler_melden("parspl()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }


  /* --------- print output ---------------------------------------- */

  printf("\n\n"
         "Spline coefficients for SX:\n\n"
          " i        ax[i]              bx[i]              "
          "cx[i]              dx[i]\n"
          "------------------------------------------------"
          "-----------------------------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f"
           "  %17.13"LZP"f\n", i, x[i], b[i], c[i], d[i]);

  printf("\n\n"
         "Spline coefficients for SY:\n\n"
          " i        ay[i]              by[i]              "
          "cy[i]              dy[i]\n"
          "------------------------------------------------"
          "-----------------------------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f"
           "  %17.13"LZP"f\n", i, y[i], by[i], cy[i], dy[i]);


  /* --- print function and derivative values for 5 points -------- */

  printf("\n\n"
         "some function and derivative values:\n");
  xmit = (t[n - 1] - t[0]) / (REAL)7.0;
  for (i = 1; i <= 5; i++)
  {
    pspwert(n - 1, t[0] + i * xmit, t, x, b, c, d, y, by, cy, dy,
            &hilf1, &hilf2, ausp);
    printf("\nvalue of spline at t = %20.13"LZP"f:\n  "
           "(%20.13"LZP"f,  %20.13"LZP"f)\n",
           t[0] + i * xmit, hilf1, hilf2);
    printf("derivatives for the spline at t = \n");
    printf("%20.13"LZP"f %20.13"LZP"f %20.13"LZP"f\n",
           ausp[1][0], ausp[2][0], ausp[3][0]);
    printf("%20.13"LZP"f %20.13"LZP"f %20.13"LZP"f\n",
           ausp[1][1], ausp[2][1], ausp[3][1]);
  }


  /* tabulate the spline from t[0] to t[n-1] with minimal step size   */

  laenge = t[n - 1] - t[0];
  hilf1 = t[0]     - (REAL)0.015 * laenge;
  hilf2 = t[n - 1] + (REAL)0.015 * laenge;

  fehler = partab(n - 1, hilf1, hilf2,
                  (hilf2 - hilf1) / (TABELLENLAENGE - 1 - n),
                  TABELLENLAENGE - 1,
                  t, x, b, c, d, y, by, cy, dy, xtab, ytab, &lentab);

  if (fehler != 0)
  {
    fehler_melden("partab()", 20 + fehler, __FILE__, __LINE__);
    return 20 + fehler;
  }

  printf("\n\n"
         "table of values for spline from about t[0] to t[n] "
         "with minimal step size:\n\n"
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
        fehler_melden("zeigkurv(): lack of storage space",
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
