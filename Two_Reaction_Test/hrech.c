#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vminit, PMATRIX */
#include <bikub.h>

/*
   Test program: Compute arbitrary points on a Bezier surface:
                 standard and modified methods
*/

#define get_vector(x) scanf ("%"LZS"f %"LZS"f %"LZS"f", \
                             &(x)[i][j][0],&(x)[i][j][1],&(x)[i][j][2])

int main (int argc, char *argv[])
{
  #define pMax 30
  int  n, m, i, j, typ, num;
  REAL ***b, ***d, **points, wp, vp, eps;
  int  fehler;
  void *vmblock;  /* liste of dynamically allocated vectors and       */
                  /* matrices                                         */

  if ((fehler = umleiten(argc, argv))    /* assign input/output files */
      != 0)                              /* to standard ones          */
    return fehler;  /* 1 or 2 */

  printf ("Bezier Splines (Points)\n");

  printf ("modified or standard method?\n");
  printf ("   1     or    0     ?:\n");
  scanf ("%d", &typ);

  printf ("number of patches in x- and y-direction?\n");
  printf ("put in m and n :\n");
  scanf ("%d %d", &m,&n);

  printf ("m = %d, n = %d\n", m, n);

  /* ---------- allocate memory for the matrices of points ---------- */

  vmblock = vminit();                      /* initialize memory block */
  d      = (REAL ***)vmalloc(vmblock, PMATRIX, m + 1,     n + 1);
  b      = (REAL ***)vmalloc(vmblock, PMATRIX, 3 * m + 1, 3 * n + 1);
  points = (REAL **) vmalloc(vmblock, MATRIX,  pMax,      3);
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
    printf ("Read in error bound for desired accuracy\n");
    scanf ("%"LZS"f",&eps);
  }
  else
  {
    printf ("Coordinates of weight pointsd:\n");
    for (i=0; i<=m; i++)
      for (j=0; j<=n; j++)
        get_vector (d);
  }


  bezier (b, d, typ, m, n, eps);

  do {
    printf ("Number of points and number of"
            " parameterline (first direction): \n");
    scanf  ("%d %"LZS"f", &num, &wp);
    if (num < 0 || num > pMax)
      printf ("number does not exceed %d!\n", pMax);
  } while (num < 0 || num > pMax);


  rechwp (b,m,n,wp,num,points);

  printf ("for num = %d  and  wp = %5.2"LZP"f:\n",num,wp);
  for (i=0; i<=num-1; i++)
  printf ("p[%d] = %10.5"LZP"f  %10.5"LZP"f  %10.5"LZP"f\n",
          i, points [i][0], points [i][1], points [i][2]);

  do {
    printf ("number of points and number of "
            "parameter line (second direction):\n");
    scanf ("%d %"LZS"f",&num,&vp);
    if (num < 0 || num > pMax)
      printf ("number cannot exceed %d!\n", pMax);
  } while (num < 0 || num > pMax);

  rechvp (b,m,n,vp,num,points);

  printf ("for num = %d  and  vp = %5.2"LZP"f:\n",num,vp);
  for (i=0; i<=num-1; i++)
    printf ("p[%d] = %10.5"LZP"f  %10.5"LZP"f  %10.5"LZP"f\n",
             i, points [i][0], points [i][1], points [i][2]);

  return 0;
}
