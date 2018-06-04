#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vminit, PMATRIX */
#include <bikub.h>
#include <zeigkrv2.h>   /*  for  zeigflaeche                          */

/*
   Test program: two-dimensional  Bezier splines:
                 standard and modified method
*/

#define get_vector(x) scanf ("%"LZS"f %"LZS"f %"LZS"f", \
                             &(x)[i][j][0],&(x)[i][j][1],&(x)[i][j][2])

int main (int argc, char *argv[])
{
  int  n, m, i, j, typ;
  REAL ***b, ***d, eps,
       ***c,      /* [0..nv-1,0..nw-1,0..2] matrix with points        */
                  /* on the spline surface                            */
       v;
  int  nv,        /* number of v curve points to be computed          */
       nw;        /* number of w curve points to be computed          */
  int  fehler;
  void *vmblock;  /* list of dynamically allocated vectors and        */
                  /* matrices                                         */

  if ((fehler = umleiten(argc, argv))    /* assign input/output files */
      != 0)                              /* to standard ones          */
    return fehler;  /* 1 or 2 */

  printf ("Bezier Splines (bezier, intpol, bezpkt)\n");

  nv = 20;
  nw = 50;

  printf ("modified or standard method?\n");
  printf ("  1      or    0      ?:\n");
  scanf ("%d", &typ);

  printf ("number of patches in x- and y-direction?\n");
  printf ("Put in m and n:\n");
  scanf ("%d %d", &m,&n);

  printf ("m = %d, n = %d\n", m, n);

  /* ---------- allocate memory for the matrices of points ---------- */

  vmblock = vminit();                      /* initialize memory block */
  d = (REAL ***)vmalloc(vmblock, PMATRIX, m + 1,     n + 1);
  b = (REAL ***)vmalloc(vmblock, PMATRIX, 3 * m + 1, 3 * n + 1);
  c = (REAL ***)vmalloc(vmblock, PMATRIX, nv,        nw);
  if (! vmcomplete(vmblock))   /* allocations partially unsuccessful? */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  printf ("coordinates in R^3 of the Bezier points"
          " b[k][i][j] (k=0,1,2):\n");

  printf ("for i=0,...,%d and j=0:\n", 3*m);
  for (i=0,  j=0; i <= 3*m; i++) get_vector (b);

  printf ("for i=0,...,%d and j=%d:\n", 3*m, 3*n);
  for (i=0, j=3*n; i<=3*m; i++) get_vector (b);

  printf ("for i=0 and j=0,...,%d\n", 3*n);
  for (i=0, j=0; j<=3*n; j++) get_vector (b);

  printf ("for i=%d and j=0,...,%d\n", 3*m, 3*n);
  for (i=3*m, j=0; j<=3*n; j++) get_vector (b);

  if (typ == 1)
  {
    for (i=3; i<=3*m-3; i+=3)
    {
      printf ("for i=%d and j=3,6,...,%d,%d\n", i, 3*n-6, 3*n-3);
      for (j=3; j<=3*n-3; j+=3) get_vector (b);
    }
    printf ("accuracy bound for modification = ");
    scanf ("%"LZS"f",&eps);
  }
  else
  {
    printf ("coordinates of the weight points:\n");
    for (i=0; i<=m; i++)
      for (j=0; j<=n; j++)
        get_vector (d);
  }

  i = bezier (b, d, typ, m, n, eps);
  printf ("\nReturn value of \"bezier\": %d", i);

  printf ("\nOutput:");
  for (i=0; i<=3*m; i++)
  {
    printf ("\n\n");
    for (j=0; j<=3*n; j++)
#if 0
      printf ("%10.7"LZP"f  %10.7"LZP"f  %10.7"LZP"f\n",
               b[i][j][0], b[i][j][1], b[i][j][2]);
#else
      printf ("%10.7"LZP"f  %10.7"LZP"f  %10.7"LZP"f    "
              "L = %-10.4"LZP"g\n",
               b[i][j][0], b[i][j][1], b[i][j][2],
         SQRT (b[i][j][0] * b[i][j][0] +
               b[i][j][1] * b[i][j][1] +
               b[i][j][2] * b[i][j][2]));
#endif
  }

  /* ------------- compute points on the spline surface ------------- */

  for (i = 0; i < nv; i++)
  {
    v = (REAL)i / (REAL)(nv - 1);
    rechvp(b, m, n, v, nw, c[i]);
  }

  /* ------------ plot spline surface if desired          ----------- */

  fehler = 0;
  if (argc <= 3 || *argv[3] != 'n')     /* plot not suppressed?       */
    fehler = zeigflaeche(c, nv, nw, d, m + 1, n + 1, 1, 1, 1);

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
      fehler_melden("zeigflaeche(): graphics error",
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

  return 0;
}
