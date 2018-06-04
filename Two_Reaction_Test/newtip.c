#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE newtip.c ------------------------ */

#include <basis.h>
#include <newtip.h>
/*.BA*/

int newtip (
/*.IX{newtip}*/
            int   n,
            REAL* x,
            REAL* y,
            REAL* b
           )
/***********************************************************************
* Computes the coefficients of the interpolating polynomial in Newton's*
* notation.                                                            *
* Subsequently the function  valnip can be used to evaluate the        *
* polynomial at specified values.                                      *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   int    n        degree of the Newton interpolating polynomial      *
*                   (= number of nodes - 1)                            *
*   REAL   x[],y[]  the n+1 nodes given by their in x-values, y-values *
*   REAL   b[]      the interpolating polynomial coefficients :        *
*                     p[t] = b[0] +                                    *
*                            + b[1] * (t-x[0]) +                       *
*                            + b[2] * (t-x[0]) * (t-x[1]) +            *
*                            + ...                                     *
*                            + b[n] * (t-x[0]) *...* (t-x[n-1]) .      *
*                                                                      *
* Return value :                                                       *
*   0:              no error                                           *
*   1:              n is negative                                      *
*   2:              two x-values coincide                              *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i, k;
  REAL h;

  if (n < 0) return (1);                  /* number of nodes negative */

  for (i = 0; i <= n; i++)
    b [i] = y [i];

  for (i = 1; i <= n; i++)
    for (k = n; k >= i; k--)
    {
      h = x [k] - x [k-i];
      if (h == ZERO)                           /* x-values distinct ? */
        return (2);
      b [k] = (b [k] - b [k-1]) / h;
    }
  return (0);
}
/*.BA*/

REAL valnip (
/*.IX{valnip}*/
             REAL  t,
             REAL* x,
             REAL* b,
             int   n
            )
/***********************************************************************
* REAL valnip (t, x, b, n) computes the function value of an interpol- *
* ating polynomial with the coefficients b[i]  for the nodes x[i],     *
* i=0,...,n, at t via the generalized Horner scheme.                   *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   REAL t       place where polynomial shall be evaluated             *
*   REAL x[]     x-values of nodes x[i], i=0,...,n                     *
*   REAL b[]     coefficients of polynomial in Newton's notation:      *
*                   p[t] = b[0] +                                      *
*                          + b[1] * (t-x[0]) +                         *
*                          + b[2] * (t-x[0]) * (t-x[1]) +              *
*                          + ...  +                                    *
*                          + b[n] * (t-x[0]) *...* (t-x[n-1]) .        *
*   int n        degree of ther polynomial                             *
*                                                                      *
* Reurn value : value of the polynomial at t                           *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i;
  REAL v = b [n];
  for (i = n - 1; i >= 0; i--)
    v = v * (t - x [i]) + b [i];
  return (v);
}

/* -------------------------- END newtip.c -------------------------- */
