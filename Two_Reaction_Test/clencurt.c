#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 15.9}{Clenshaw-Curtis Quadrature Formulas}
             {Clenshaw-Curtis Quadrature Formulas}*/

/*.BE*/
/* ------------------------ MODULE clencurt.c ----------------------- */

#include <basis.h>
#include <gax.h>
/*.BA*/

int ClenCurt (
/*.IX{ClenCurt}*/
              REAL  func (REAL),
              int   m,
              REAL* t,
              int   n,
              REAL* Tk,
              REAL* Ak,
              REAL* Resultat
             )
/***********************************************************************
* Computes the value of the integral of func (x) over the interval     *
* (a,b) with the partition                                             *
*              t: a = t[0] < t[1] < .. < t[m] = b                      *
* by using the summed  Clenshaw-Curtis formula.                        *
* This program uses precomputed [0..n] weight vectors and Chebyshev    *
* nodes.                                                               *
.BE*)
*                                                                      *
* Parameter:                                                           *
*   REAL   func (REAL  )  Name of integrating function                 *
*   int    m              number of sun=b intervals                    *
*   REAL   t []           partition                                    *
*   int    n              n + 1 = number of nodes, n > 1, n even       *
*                                   (n + 2 = global error order)       *
*   REAL   Tk []          Chebyshev nodes                              *
*   REAL   Ak []          Weights                                      *
*                           Tk and Ak must be made available before    *
*                           calling this function by the procedure     *
*                           ClenCurtStGew for example                  *
*   REAL   *Resultat      Compute integral value                       *
*                                                                      *
* REturn value :                                                       *
*   0:                    o.k.                                         *
*   1:                    improper number of nodes                     *
*   2:                    improper number of sub intervals             *
*                                                                      *
* Author:                 Uli Eggermann, 10.3.1991                     *
.BA*)
***********************************************************************/
/*.BE*/
{
  int    j, k;
  REAL   v, h, sum;

  if (n < 2 || n % 2 != 0) return 1;        /* n positive ? n even  ? */
  if (m < 1)               return 2;        /* partition              */

  for (*Resultat = ZERO, j = 0; j < m; j++) /* loop over intervals    */
  {
    v = 0.5 * (t [j+1] - t [j]);            /* half the interval size */
    h = 0.5 * (t [j+1] + t [j]);            /* Interval center        */

    for (sum = ZERO, k = 0; k <= n; k++)    /* Chebyshev loop         */
      sum += Ak [k] * func (v * Tk [k] + h);

    *Resultat += v * sum;
  }
  return 0;
}

/* ------------------------- END clencurt.c ------------------------- */
