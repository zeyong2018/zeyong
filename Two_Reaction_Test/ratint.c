#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE ratint.c ------------------------- */

#include <basis.h>
#include <vmblock.h>
#include <ratint.h>
/*.BA*/

int ratint (
/*.IX{ratint}*/
            int n,
            int num,
            REAL* x,
            REAL* y,
            int* md,
            REAL eps
           )
/***********************************************************************
* This program attempts to find a rational interpolation for a given   *
* set of nodes (x[i], y[i]), i=0,...,n, and the prescribed numerator   *
* degree  num.                                                         *
* If the given data allows a rational interpolant, then it is uniquely *
* determined for the numerator degree.                                 *
* Here, as in all interpolation problems, the x[i] must be distinct.   *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   int    n       n+1 is the number of nodes (n > 1)                  *
*   int    num     degree of numerator polynomial (num <= n)           *
*   REAL   x[]     nodes : x-values                                    *
*   REAL   y[]             y-values                                    *
*   int    md[]    Multiplication or division when evaluating via      *
*                  Horner scheme                                       *
*   REAL   eps     accuracy for interpolation                          *
*                                                                      *
* Return value :                                                       *
*   0: no error                                                        *
*   1: num > n  or  n < 1                                              *
*   2: two x-values coincide                                           *
*   3: lack of memory for aux arrays                                   *
*   4: Interpolating function does not exist :                         *
*        try with a numerator degree  num > n / 2                      *
*   5: number of nodes and degree of denominator polynomial has become *
*        negative: increase numerator degree                           *
*   6: degree of denominator polynomial negative: change numerator deg.*
*   7: not all nodes were used for the interpolation, change numerator *
*        degree                                                        *
*                                                                      *
* Locally used subroutines:                                            *
*   ratval, sel_ymin                                                   *
*                                                                      *
* Constants used:                                                      *
*   NULL, MACH_EPS                                                     *
*                                                                      *
* Using the notation                                                   *
*              u[i] := (t - x[i])                                      *
* the rational function has the 'continued fraction' representation :  *
*                                                                      *
*                               u[n]*u[n-1]                            *
*   f(t) = y[n] + u[n]*y[n-1] + -------------------------------------- *
*                                                        u[n-2]*u[n-3] *
*                               y[n-2] + u[n-2]*y[n-3] + ------------- *
*                                                        y[n-4] + ...  *
.BA*)
***********************************************************************/
/*.BE*/

#define MULTIPLY   0
#define DIVIDE     1

{
  int    i, j, j1, denom, nend, ret;
  REAL xj, yj, y2, x2;
  REAL *x1, *y1, *z;
  void   sel_ymin (int n, REAL* xj, REAL* yj, REAL* x, REAL* y);
  REAL   ratval   (int n, REAL  x0, REAL* x,  REAL* y, int*  md);
  void   *vmblock;

  if (num > n || n < 1) return (1);

  for (i = 0; i < n; i++)
    for (j = i + 1; j <= n; j++)
      if (x[i] == x[j]) return (2);

  vmblock = vminit();
  x1 = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  y1 = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  z  = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 3;
  }

  eps = max (eps, (REAL)128.0 * MACH_EPS);

  for (i = 0; i <= n; i++)          /* copy nodes                     */
  {
    x1 [i] = x [i];
    y1 [i] = y [i];
    md [i] = MULTIPLY;              /* initialize mult/divide info    */
  }
  nend  = n;
  denom = n - num;        /* denom = degree of denominator polynomial */
  if (num < denom)
  {
    for (i = 0; i <= n; i++)
      if (y [i] != ZERO)   y [i] = ONE / y [i];
      else                 { ret = 4; goto FreeAll; }
    md [n] = DIVIDE;
    j = num; num = denom; denom = j;
  }

  while (nend > 0)
  {
    for (i = 1; i <= num-denom; i++)      /* numerator degree exceeds */
    {                                     /* denominator degree :     */
      xj = x [nend];                      /*  use divided differences */
      yj = y [nend];
      for (j = 0; j < nend; j++)
        y [j] = (y [j] - yj) / (x [j] - xj);
      nend--;
    }
    if (nend < 0 && denom < 0) { ret = 5; goto FreeAll; }
    if (nend > 0)
    {
      sel_ymin (nend, &xj, &yj, x, y);
      for (j1 = j = 0; j < nend; j++)
      {
        y2 = y [j] - yj;
        x2 = x [j] - xj;
        if (FABS (y2) <= FABS(x2) * eps)    /* automatically          */
          z [j1++] = x [j];                 /* interpolated nodes     */
        else
        {
          y [j-j1] = x2 / y2;               /* Interpolate by forming */
          x [j-j1] = x [j];                 /* inverse divided        */
        }                                   /*          differences   */
      }
      for (j = 0; j < j1; j++)
      {
         x [nend-1] = z [j];                /* join automatically     */
         y [nend-1] = ZERO;                 /* interpolated nodes with*/
         for (i = 0; i < nend; i++)         /* those from divided     */
           y[i] *= x[i] - x[nend];          /* differences            */
         nend--;
      }
      if (nend > 0)
      {
        md [--nend] = DIVIDE;            /* Compute new final index,  */
        num   = denom;                   /* numerator and denominator */
        denom = nend - num;              /* degrees                   */
      }
      if (denom < 0 && nend < 0) { ret = 6; goto FreeAll; }
    }
  }
  y2 = FABS (y[n]);                /* At the end check whether all    */
  for (i = 0; i < n; i++)          /* nodes were used                 */
    y2 += FABS (y[i]);
  for (i = 0; i <= n; i++)
  {
    x2 = ratval (n, x1 [i], x, y, md);
    if (FABS (x2 - y1[i]) > n * eps * y2) { ret=7; goto FreeAll; }
  }
  ret = 0;
FreeAll:
  vmfree(vmblock);
  return (ret);
}

void sel_ymin (
/*.IX{sel\unt ymin}*/
                int   nend,
                REAL* xj,
                REAL* yj,
                REAL* x,
                REAL* y
              )
/***********************************************************************
* We find the node with the smallest magnitude y-value and store it in *
* (xj, yj). This node, indexed by nend is the next node to be inter-   *
* polated in retint.                                                   *
***********************************************************************/
{
  int j, k;
  *yj = y [nend];                    /* search for minimal modulus yj */
  j = nend;                          /* and swap with  y[nend]        */
  for (k = 0; k < nend; k++)
    if (FABS (y [k]) < FABS (*yj)) *yj = y [j = k];
    *xj = x [j]; x [j] = x [nend]; x [nend] = *xj;
                 y [j] = y [nend]; y [nend] = *yj;
}
/*.BA*/

REAL ratval (
/*.IX{ratval}*/
             int   n,
             REAL  x0,
             REAL* x,
             REAL* y,
             int*  md
            )
/***********************************************************************
* Evaluates a rational interpolating function from its coefficients    *
* as determined in  ratint.                                            *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   int    n        n+1 is the number of nodes (n > 1)                 *
*   REAL   t        place for the evaluation                           *
*   REAL   x[]      x-values of nodes                                  *
*   REAL   y[]      coefficients of the interpolating rat. function    *
*   int    md[]     information on mult/divide in Horner scheme        *
*                                                                      *
* Return value :                                                       *
*   functional value at t,                                             *
*   this is set equal to  MAXROOT, if division by zero should occur    *
*                                                                      *
* constants used :                                                     *
*   EPSQUAD, MAXROOT                                                   *
.BA*)
***********************************************************************/
/*.BE*/
{
  int i;
  REAL res = y [0];
  for (i = 1; i <= n; i++)
  {
    if (md [i-1] == MULTIPLY)
           res = y [i] + (x0 - x [i]) * res;
    else if (FABS (res) > EPSQUAD)
           res = y [i] + (x0 - x [i]) / res;
         else
           return (MAXROOT);
  }
  if (md [n] == MULTIPLY) return (res);
  return ((FABS (res) > EPSQUAD) ? (ONE / res) : MAXROOT);
}

#undef MULTIPLY
#undef DIVIDE

/* --------------------------- END ratint.c ------------------------- */
