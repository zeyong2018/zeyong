#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>      /*  for  REAL, umleiten, fprintf, stderr,     */
                        /*       scanf, LZS, readln, fehler_melden,   */
                        /*       printf, LZP                          */
#include <vmblock.h>    /*  for  vminit, vmalloc, VEKTOR, vmcomplete, */
                        /*       PMATRIX, vmfree                      */
#include <shepard.h>    /*  for  shepard                              */
#include <zeigkrv2.h>   /*  for  zeigflaeche                          */



/*--------------------------------------------------------------------*/

int main
        (
         int  argc,
         char *argv[]
        )

/***********************************************************************
* Test program for the function shepard() for computing one point on   *
* Shepard interpolation surface.                                       *
*                                                                      *
* Operation of the program:                                            *
* =========================                                            *
* Input from file stdin, output to file stdout. If the first command   *
* line parameter exists, it is assigned to stdin. Analogously for the  *
* second command line parameter and stdout.                            *
* Calls for input and error messages are routed to stderr.             *
* If the macro INTERAKTIV is not defined, we issue no calls for input. *
*                                                                      *
* The input data is put out, then shepard() computes interpolation     *
* points for a rectangular grid. For a Turbo C compiler, we produce    *
* graphics output of the surface, unless the third command line        *
* parameter has been given the value "n".                              *
*                                                                      *
* Some test data sets are listed in mshepard.ei*. This data should     *
* produce the same output as in mshepard.au* .                         *
*                                                                      *
* Input file:                                                          *
* ===========                                                          *
* n         number of nodes - 1                                        *
* mue       Shepard parameter                                          *
* fx[0]     x-coordinate of  0th node                                  *
* ...       ...                                                        *
* fx[n]     x-coordinate of  nth node                                  *
* fy[0]     y-coordinate of  0th node                                  *
* ...       ...                                                        *
* fy[n]     y-coordinate of  nth node                                  *
* fz[0]     0th function value                                         *
* ...       ...                                                        *
* fz[n]     nth function value                                         *
* xa  \     [xa,xb] interpolation interval in x direction.             *
* xb  /                                                                *
* xs        Step size for [xa,xb]                                      *
* ya  \     [ya,yb] interpolation interval in y direction.             *
* yb  /                                                                *
* ys        Step size for [ya,yb]                                      *
* methode   Number of desired variante of Shepard interpolation        *
* rr        Shepard radius for the local ethod                         *
***********************************************************************/

{
  int  n;
  REAL mue;
  REAL *fx;
  REAL *fy;
  REAL *fz;
  REAL xa;
  REAL xb;
  REAL xs;
  REAL ya;
  REAL yb;
  REAL ys;
  int  methode;
  REAL rr;
  REAL x;
  REAL y;
  REAL z;
  int  fehler;      /* error code from umleiten(), shepard()          */
  void *vmblock;    /* List of dynamically allocated Vectors          */
  int  i;
  REAL ***d;
  REAL ***flaeche;
  REAL ymin;
  REAL ymax;
  REAL xmin;
  REAL xmax;
  int  xtablen;
  int  ytablen;
  int  j;


  if ((fehler = umleiten(argc, argv))/* assign input/output files     */
      != 0)                          /* to standard in/output files   */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input data -------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
          "number of nodes - 1:                          ");
#endif
  scanf("%d", &n);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr,
          "Shepard parameter mue:                                 ");
#endif
  scanf("%"LZS"f", &mue);
  readln();


  /* ------------------ allocate dynamic vectors ----------------- */

  vmblock = vminit();                 /* initialize  storage         */
  fx = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  fy = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  fz = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  if (! vmcomplete(vmblock))                    /* lack of memory ?  */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i <= n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "x coordinate of %2d. node:"
                    "                     ", i);
#endif
    scanf("%"LZS"f", fx + i);
  }

  for (i = 0; i <= n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "y coordinate of %2d. node:"
                    "                     ", i);
#endif
    scanf("%"LZS"f", fy + i);
  }

  for (i = 0; i <= n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "%2d. function value:                   "
                    "                     ", i);
#endif
    scanf("%"LZS"f", fz + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr,
          "left endpoint of x-interpolation interval:            ");
#endif
  scanf("%"LZS"f", &xa);
  readln();

#ifdef INTERAKTIV
  fprintf(stderr,
          "right endpoint of x-interpolation interval:           ");
#endif
  scanf("%"LZS"f", &xb);
  readln();

#ifdef INTERAKTIV
  fprintf(stderr,
          "Step size in x direction:                            ");
#endif
  scanf("%"LZS"f", &xs);
  readln();

#ifdef INTERAKTIV
  fprintf(stderr,
          "left endpoint of y-interpolation interval:            ");
#endif
  scanf("%"LZS"f", &ya);
  readln();

#ifdef INTERAKTIV
  fprintf(stderr,
          "right endpoint of y-interpolation interval:           ");
#endif
  scanf("%"LZS"f", &yb);
  readln();

#ifdef INTERAKTIV
  fprintf(stderr,
          "Step size in y direction:                            ");
#endif
  scanf("%"LZS"f", &ys);
  readln();

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
          "global (0), local (1), Franke-Little weights (2):     ");
#endif
  scanf("%d", &methode);
  readln();

  if (methode == 1 || methode == 2)
  {
#ifdef INTERAKTIV
    fprintf(stderr,
            "Shepard radius:                                        ");
#endif
    scanf("%"LZS"f", &rr);
    readln();
  }
  else
    rr = ZERO;


  /* ------------ put out input for checking ----------- */

  printf("\n\n"
         "Surface representation using Shepard interpolation\n"
         "==================================================\n\n"
         "Input data:\n"
         "===========\n\n"
         "number of nodes - 1:             %15d\n"
         "Shepard parameter mue:           %15.2"LZP"g\n\n"
         "nodes and function values:\n"
         "  i    fx[i]       fy[i]       fz[i]\n"
         "---------------------------------------\n",
         n, mue
        );

  for (i = 0; i <= n; i++)
    printf("%3d%12.7"LZP"g%12.7"LZP"g%12.7"LZP"g\n",
           i, fx[i], fy[i], fz[i]);

  printf("\n"
         "x-interval [xa,xb]:              [%6.2"LZP"g,%6.2"LZP"g]\n"
         "x step size:                     %15.2"LZP"g\n"
         "y-interval [xa,xb]:              [%6.2"LZP"g,%6.2"LZP"g]\n"
         "y step size:                     %15.2"LZP"g\n"
         "Method:                          %15d\n",
         xa, xb, xs, ya, yb, ys, methode
        );


  /* --------------- compute interpolating surface -------------- */

  printf("\n"
         "Output data:\n"
         "===========\n\n"
         "Points on the interpolating surface (z = PHI(x,y)):\n"
         "       x           y           z\n"
         "------------------------------------\n");
  for (x = xa; x <= xb; x += xs)
    for (y = ya; y <= yb; y += ys)
    {
      fehler = shepard(x, y, fx, fy, fz, n, mue, methode, rr, &z);
      if (fehler)
        fehler_melden("shepard()", 10 + fehler, __FILE__, __LINE__);
      else
        printf("%12.7"LZP"g%12.7"LZP"g%12.7"LZP"g\n", x, y, z);
    }


  /* ----------- graphic display of surface ---------- */

  if (argc <= 3 || *argv[3] != 'n')     /* Graphic output desired? */
  {

    xtablen = 40;
    ytablen = 20;

    d       = (REAL ***)vmalloc(vmblock, PMATRIX, 1,       n + 1);
    flaeche = (REAL ***)vmalloc(vmblock, PMATRIX, xtablen, ytablen);
    if (! vmcomplete(vmblock))
    {
      fehler_melden("lack of memory", 0, __FILE__, __LINE__);
      return 3;
    }

    xmin = fx[0];
    xmax = fx[0];
    ymin = fy[0];
    ymax = fy[0];
    for (i = 1; i <= n; i++)
      xmin = (fx[i] < xmin) ? fx[i] : xmin,
      xmax = (fx[i] > xmax) ? fx[i] : xmax,
      ymin = (fy[i] < ymin) ? fy[i] : ymin,
      ymax = (fy[i] > ymax) ? fy[i] : ymax;

    xmin -= (xmax - xmin) / TEN;
    xmax += (xmax - xmin) / TEN;
    ymin -= (ymax - ymin) / TEN;
    ymax += (ymax - ymin) / TEN;

    for (i = 0; i <= n; i++)
      d[0][i][0] = fx[i],
      d[0][i][1] = fy[i],
      d[0][i][2] = fz[i];

    for (x = xmin, i = 0; i < xtablen;
         i++, x += (xmax - xmin) / (REAL)xtablen)
      for (y = ymin, j = 0; j < ytablen;
           j++, y += (ymax - ymin) / (REAL)ytablen)
      {
        flaeche[i][j][0] = x;
        flaeche[i][j][1] = y;
        fehler = shepard(x, y, fx, fy, fz, n, mue, methode, rr,
                         &flaeche[i][j][2]);
        if (fehler)
          fehler_melden("shepard()", 20 + fehler, __FILE__, __LINE__);
      }

    fehler = zeigflaeche(flaeche, xtablen, ytablen, d, 1, n + 1,
                         1, 1, 0);

    switch (fehler)
    {
      case 0:
        break;
      case 1:
        fehler_melden("zeigflaeche(): nv too small",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 2:
        fehler_melden("zeigflaeche(): nw too small",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 3:
        fehler_melden("zeigflaeche(): m too small",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 4:
        fehler_melden("zeigflaeche(): n too small",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 5:
        fehler_melden("zeigflaeche(): nw or n too large",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 6:
        fehler_melden("zeigflaeche(): Graphics error",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 7:
        fehler_melden("zeigflaeche(): lack of memory",
                      30 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("zeigflaeche(): other error",
                      30 + fehler, __FILE__, __LINE__);
    }
  }


  vmfree(vmblock);

  return 0;
}
