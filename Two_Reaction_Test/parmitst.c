#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, fprintf,   */
                         /*       stderr, NULL, REAL, fehler_melden,  */
                         /*       LZS, LZP                            */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, VEKTOR */
#include <hrmsplin.h>    /*  for  parmit                              */
#include <spliwert.h>    /*  for  pmtwert                             */
#include <splintab.h>    /*  for  pmtab                               */
#include <zeigkurv.h>    /*  for  zeigkurv                            */



#define TABELLENLAENGE  140



/***********************************************************************
*                                                                      *
* Test program: parametric interpolating polynomial spline of degree 5:*
*               compute coefficients and tabulate function values;     *
*               plot for Turbo C compiler                              *
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
* richtung          (control for the kind of prescribed derivatives:   *
*                    1: Tangent vectors                                *
*                    2: Normal vectors                                 *
*                    3: derivative dy/dx)                              *
*                                                                      *
* x1 y1  xt1 yt1    (nodes and tangent vectors                         *
* x2 y2  xt2 yt2     ...                      )                        *
* ...    ...        ...                                                *
* xn yn  xtn ytn    (ditto )                                           *
*                or                                                    *
* x1 y1  xn1 yn1    (nodes and normal vectors                          *
* x2 y2  xn2 yn2     ...                     )                         *
* ...    ...        ...                                                *
* xn yn  xnn ynn    (ditto )                                           *
*                or                                                    *
* x1 y1  y1'        (nodes and first derivative wrt. x                 *
* x2 y2  y2'         ...                              )                *
* ...    ...        ...                                                *
* xn yn  yn'        (ditto )                                           *
*                                                                      *
* mass              (type of end point condition :                     *
*                    1: periodic                                       *
*                    2: natural                                        *
*                    3: 2nd end point derivatives of y wrt x           *
*                    4: 2nd end point derivatives of x and y wrt. t    *
*                    5: curvature radii at end points                  *
*                    6: 3rd end point derivative of x and y wrt t)     *
* alfax alfay       (left end condition)  (for mass = 4,6 only)        *
* betax betay       (right end condition) (for mass = 4,6 only)        *
* alfa  beta        (left and right conditions) (for mass = 3,5 only)  *
*                                                                      *
***********************************************************************/

int main(int argc, char *argv[])

{
  REAL     *x, *y,
           *xricht, *yricht,
           *xt, *cx, *dx, *ex, *fx,
           *yt, *cy, *dy, *ey, *fy,
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
           richt,
           lentab;
  abl_mat2 ausp;
  void     *vmblock;     /* List of dynamically allocated vectors    */


  if ((fehler = umleiten(argc, argv))   /* assign  input/output files */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "number of nodes:                               ");
#endif
  scanf("%d", &n);

  vmblock = vminit();                          /* initialize storage  */
  x      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  y      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  xricht = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  yricht = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  xt     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  cx     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  dx     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  ex     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  fx     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  yt     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  cy     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  dy     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  ey     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  fy     = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  t      = (REAL *)vmalloc(vmblock, VEKTOR, n,              0);
  xtab   = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  ytab   = (REAL *)vmalloc(vmblock, VEKTOR, TABELLENLAENGE, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "Control number for type of derivative\n"
                  "(1 = Tangent vector, 2 = Normal vector,\n"
                  "3 = derivative y'(x)):                         ");
#endif
  scanf("%d", &richt);
  if (richt == 1 || richt == 2)
    for (i = 0; i < n; i++)
    {
#ifdef INTERAKTIV
      fprintf(stderr, "x-, y-coordinates, x-, y-components of the\n"
                      "tangent or normal vector:                  ");
#endif
      scanf("%"LZS"f%"LZS"f%"LZS"f%"LZS"f",
            x + i, y + i, xricht + i, yricht + i);
    }
  else if (richt == 3)
    for (i = 0; i < n; i++)
    {
#ifdef INTERAKTIV
      fprintf(stderr, "x-, y-coordinates, derivative:             ");
#endif
      scanf("%"LZS"f%"LZS"f%"LZS"f", x + i, y + i, yricht + i);
    }
  else
  {
    fehler_melden("wrong value for control number",
                  0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "type of end point condition (1...6):           ");
#endif
  scanf("%d", &mass);

  if (mass >= 4 && mass <= 6)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "end point conditions (alfa[0],alfa[1],\n"
                    "                      beta[0],beta[1]):      ");
#endif
    scanf("%"LZS"f%"LZS"f%"LZS"f%"LZS"f",
          alfa, alfa + 1, beta, beta + 1);
  }
  else if (mass == 3)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "boundary conditions (alfa[0],beta[0]):       ");
#endif
    scanf("%"LZS"f%"LZS"f", alfa, beta);
  }


  /* ------------ print out input ---------------------------------- */

  printf("\n"
         "Compute and tabulate an interpolating "
         "parametric Hermite spline\n"
         "======================================"
         "=========================\n"
        );

  printf("\n\n"
         "number of nodes:              %d (= n+1)\n"
         "type of end point conditions: %d\n",
         n, mass
        );
  if (mass == 4 ||        /* prescribed second or third end point     */
      mass == 6)          /* derivative wrt. t ?                      */
    printf("end point conditions:    alfax = %20.13"LZP"f\n"
           "                         alfay = %20.13"LZP"f\n"
           "                         betax = %20.13"LZP"f\n"
           "                         betay = %20.13"LZP"f\n",
           alfa[0], alfa[1], beta[0], beta[1]
          );
  else if (mass > 2)            /* not a periodic or natural spline ? */
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

  printf("\n\nControl number for type of prescribed derivative: %d\n",
         richt);
  if (richt == 3)
  {
    printf("\n"
           "derivative of y wrt. x at nodes:\n\n"
           " i         y'[i]       \n"
           "-----------------------\n"
          );
    for (i = 0; i < n; i++)
      printf("%2d %20.13"LZP"f\n", i, yricht[i]);
  }
  else if (richt == 1)
  {
    printf("\n"
           "Tangent vectors at the nodes:\n\n"
           " i        xtang[i]             ytang[i]\n"
           "--------------------------------------------\n"
          );
    for (i = 0; i < n; i++)
      printf("%2d %20.13"LZP"f %20.13"LZP"f\n",
             i, xricht[i], yricht[i]);
  }
  else if (richt == 2)
  {
    printf("\n"
           "Normal vectors at the nodes:\n\n"
           " i        xnorm[i]             ynorm[i]\n"
           "--------------------------------------------\n"
          );
    for (i = 0; i < n; i++)
      printf("%2d %20.13"LZP"f %20.13"LZP"f\n",
             i, xricht[i], yricht[i]);
  }


  /* ------------- call parmit() to compute coefficients ------------ */

  fehler = parmit(n, x, y, richt, xricht, yricht, mass, alfa, beta,
                  cx, dx, ex, fx, cy, dy, ey, fy, t, xt, yt);

  if (fehler != 0)
  {
    fehler_melden("parmit()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }


  /* --------- print output --------------------------------------- */

  printf("\n\n"
         "Spline coefficients for SX:\n\n"
          " i        ax[i]              bx[i]              "
          "cx[i]\n"
          "------------------------------------------------"
          "----------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f\n",
           i, x[i], xt[i], cx[i]);
  printf("\n"
          " i        dx[i]              ex[i]              "
          "fx[i]\n"
          "------------------------------------------------"
          "----------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f\n",
           i, dx[i], ex[i], fx[i]);

  printf("\n\n"
         "Spline coefficients for SY:\n\n"
          " i        ay[i]              by[i]              "
          "cy[i]\n"
          "------------------------------------------------"
          "----------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f\n",
           i, y[i], yt[i], cy[i]);
  printf("\n"
          " i        dy[i]              ey[i]              "
          "fy[i]\n"
          "------------------------------------------------"
          "----------\n"
        );
  for (i = 0; i < n - 1; i++)
    printf("%2d %17.13"LZP"f  %17.13"LZP"f  %17.13"LZP"f\n",
           i, dy[i], ey[i], fy[i]);


  /* -------------- print parameter nodes --------------------------- */

  printf("\n\n"
         "Parameter nodes:\n\n"
         " i          t[i]       \n"
         "-----------------------\n"
        );
  for (i = 0; i < n; i++)
    printf("%2d %20.13"LZP"f\n", i, t[i]);


  /* --- print functional and derivative values at 5 places --------- */

  printf("\n\n"
         "some function and derivative values for the spline:\n");
  xmit = (t[n - 1] - t[0]) / (REAL)7.0;
  for (i = 1; i <= 5; i++)
  {
    pmtwert(n - 1, t[0] + i * xmit, t, x, xt, cx, dx, ex, fx,
            y, yt, cy, dy, ey, fy, &hilf1, &hilf2, ausp);
    printf("\nvalue of spline at t0 = %10.5"LZP"f is\n"
           " %20.13"LZP"f  %20.13"LZP"f\n",
           t[0] + i * xmit, hilf1, hilf2);
    printf("derivatives of spline = \n");
    printf("%20.13"LZP"f %20.13"LZP"f %20.13"LZP"f\n",
           ausp[1][0], ausp[2][0], ausp[3][0]);
    printf("%20.13"LZP"f %20.13"LZP"f\n", ausp[4][0], ausp[5][0]);
    printf("\n");
    printf("%20.13"LZP"f %20.13"LZP"f %20.13"LZP"f\n",
           ausp[1][1], ausp[2][1], ausp[3][1]);
    printf("%20.13"LZP"f %20.13"LZP"f \n", ausp[4][1], ausp[5][1]);
  }


  /* tabulate spline from t[0] to t[n-1] optimally ------------------ */

  laenge = t[n - 1] - t[0];
  hilf1  = t[0]     - (REAL)0.015 * laenge;
  hilf2  = t[n - 1] + (REAL)0.015 * laenge;

  fehler = pmtab(n - 1, hilf1, hilf2,
                 (hilf2 - hilf1) / (TABELLENLAENGE - 1 - n),
                 TABELLENLAENGE - 1,
                 t, x, xt, cx, dx, ex, fx, y, yt, cy, dy, ey, fy,
                 xtab, ytab, &lentab);

  if (fehler != 0)
  {
    fehler_melden("pmtab()", 20 + fehler, __FILE__, __LINE__);
    return 20 + fehler;
  }

  printf("\n\n"
         "table for spline from about t[0] to t[n] with "
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
