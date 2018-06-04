#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test program for polynomial interpolation            Egg, 11.17.1991 *
***********************************************************************/

#include <basis.h>
#include <stdio.h>
#include <math.h>
#include <newtip.h>

int main (int argc, char *argv[]) {

  #define ARR_MAX  10
  int    n, i, error;
  REAL   x [ARR_MAX+1], y [ARR_MAX+1], x0, b [ARR_MAX+1];
  char   ch[4];

  /* --- assign the input file to the standard input ---------------- */
  if (argc >= 2)                           /* at least one entry ?    */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error opening  %s!\n", argv[1]);
      return 2;
    }

  printf ("degree of the algebraic interpolation polynomial (max %d)",
           ARR_MAX);
  scanf  ("%d", &n);

  printf ("%d put in nodes in the form:  x y [return]\n", n+1);

  for (i=0; i<=n; i++)
  {
    printf ("%2d. nodes: ", i);
    scanf  ("%"LZS"f %"LZS"f", &x[i], &y[i]);
  }

  do
  {
    printf ("where do you want to evaluate the interpolation ? ");
    scanf ("%"LZS"f", &x0);

    printf ("n  = %d\n",  n);
    printf ("x0 = %"LZP"f\n", x0);
    for (i=0; i<=n; i++)
      printf ("%"LZP"f    %"LZP"f \n", x[i], y[i]);

    error = newtip (n, x, y, b);

    if (error == 0)
    {
      for (i=0; i<=n; i++)  printf ("b[%d] = %"LZP"f \n", i, b[i]);
      printf ("Value for polynomial interpolant at %"LZP"f is "
              "%"LZP"f\n", x0, valnip (x0, x, b, n));
    }
    else printf ("error code = %d\n", error);
  }
  while (printf ("another evaluation wanted ? "),
         scanf ("%s", ch),
         *ch == 'j' || *ch == 'J');

  return 0;
}


/*                 T E S T  E X A M P L E
Input:
3
0. 2.
1. 5.
2. 4.
3. 8.
2.3
j
2.4
j
2.5
n

Output:
n  = 3
x0 = 2.300000
0.000000    2.000000
1.000000    5.000000
2.000000    4.000000
3.000000    8.000000
b[0] = 2.000000
b[1] = 3.000000
b[2] = -2.000000
b[3] = 1.500000
Value for interpolant at x0 = 4.265500
*/
