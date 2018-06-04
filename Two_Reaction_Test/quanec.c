#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"



/*.BA*/
/*.FE{C 15.3}{Newton-Cotes Formulas}{Newton-Cotes Formulas}*/

/*.BE*/
/* ------------------------- MODULE quanec.c ------------------------ */

#ifndef QuaNec_C
#define QuaNec_C

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>
#include <gax.h>
/*.BA*/

int QuaNeC (
/*.IX{QuaNeC}*/
            REAL  a,
            REAL  b,
            int   AnzInt,
            int   nrV,
            REAL  f(REAL),
            REAL* Qwert,
            int*  Ordnung,
            REAL* Sch
           );
/***********************************************************************
* Compute the value of the integral of f(x) over the interval [a,b]    *
* using summed Newton-Cotes formulas.                                  *
* AnzInt denotes the number of subintervals.                           *
* If AnzInt > 1, we give an error estimate from a second run with      *
* (AnzInt-1) subintervals.                                             *
.BE*)
*                                                                      *
* Input parameters :                                                   *
*   REAL   a, b         Interval of integration                        *
*   int    AnzInt       number of subintervals                         *
*   int    nrV          Number for method used :                       *
*                         1: Trapezoidal rule                          *
*                         2: Simpson's rule                            *
*                         3: 3/8 formula                               *
*                         4: 4/90 formula                              *
*                         5: 5/288 formula                             *
*                         6: 6/840 formula                             *
*                         7: 7/17280 formula                           *
*   REAL   f (REAL)     Function                                       *
*   REAL   *Qwert       value for integral                             *
*   int    *Ordnung     error order of method used                     *
*   REAL   *Sch         error estimate for Qwert                       *
*                                                                      *
* Return value :                                                       *
*   0:                  o.k.                                           *
*   1:                  nrV improper                                   *
*   2:                  AnzInt improper                                *
*                                                                      *
* Author                Uli Eggermann, 2.2.1991                        *
.BA*)
***********************************************************************/
/*.BE*/
REAL Q_1 [] = { 1./2, 1./2 };
REAL Q_2 [] = { 1./6, 4./6, 1./6 };
REAL Q_3 [] = { 1./8, 3./8, 3./8, 1./8 };
REAL Q_4 [] = { 7./90, 32./90, 12./90, 32./90, 7./90 };
REAL Q_5 [] = { 19./288, 75./288, 50./288, 50./288,
                75./288, 19./288 };
REAL Q_6 [] = { 41./840, 216./840, 27./840, 272./840,
                27./840, 216./840, 41./840 };
REAL Q_7 [] = { 751./17280,  3577./17280, 1323./17280,
                2989./17280, 2989./17280, 1323./17280,
                3577./17280, 751./17280 };

REAL *QuadArr [] = { Q_1, Q_2, Q_3, Q_4, Q_5, Q_6, Q_7 };

int QuaNeC (REAL a,        REAL  b,      int  AnzInt,   int   nrV,
            REAL f(REAL),  REAL* _Qwert, int* _Ordnung, REAL* _Sch)
{
  #define Qwert   (*_Qwert)          /* done for readability:         */
  #define Ordnung (*_Ordnung)        /* Qwert, Ordnung and Sch are    */
  #define Sch     (*_Sch)            /* used as references            */

  int  i, j, k, n;
  REAL xk, q[2], h[2];

  if (nrV < 1 || nrV > 7) return (1);         /* incorrect method #   */
  if (AnzInt < 1)         return (2);     /* # subintervals incorrect */
  Ordnung = (nrV / 2) * 2 + 2;                /* 2, 4, 4, 6, 6, 8, 8  */
  if (a == b) { Qwert = Sch = ZERO; return (0); }          /* trivial */

  for (i = 0; i <= (AnzInt == 1 ? 0 : 1); i++)
  {
    n    = AnzInt - i;
    h[i] = (b - a) / (n * nrV);
    q[i] = ZERO;
    for (k = 0; k < n; k++)
    {
      xk = a + k * nrV * h[i];
      for (j = (k == 0 ? 0 : 1); j <= nrV; j++)
          q[i] += (QuadArr [nrV - 1] [j]
                 * f (xk + j * h[i])
                 * (((k < n - 1) && (j == nrV)) ? 2 : 1));
    }
    q[i] *= h[i] * nrV;
  }
  if (AnzInt > 1)
  {
    Qwert = q [1];
    Sch   = (q[0]-q[1]) / (POW (AnzInt/(AnzInt-ONE),(REAL)Ordnung) - 1);
  }
  else Qwert = q [0];

  return (0);

  #undef Qwert
  #undef Ordnung
  #undef Sch
}                                                           /* QuaNeC */
#endif
/*.BA*/



/*.FE{C 15.7}
     {Gau"s Quadrature Formulas}
     {Gau"s Quadrature Formulas}*/


int OrtogP (
/*.IX{OrtogP}*/
            int n,
            REAL* Integrale,
            REAL* StStelle,
            REAL* Gewicht
           )

/***********************************************************************
* Compute nodes and weights for generalized Gaussian quadrature.       *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   int    n             number of nodes                               *
*   REAL   Integrale []  vector with integral values:                  *
*                          the kth entry is the integral of the        *
*                          weighted function g(x) * x^k:               *
*                                         b       k                    *
*                          Integrale[i] = S g(x) x  , k=0,...,2n-1     *
*                                         a                            *
*   REAL   StStelle []   node vector   [0..n-1]                        *
*   REAL   Gewicht []    weight vector [0..n-1]                        *
*                                                                      *
*  Return value                                                        *
*   0:                   o.k.                                          *
*   1:                   error in Gauss                                *
*   2:                   error in Mueller                              *
*   3:                   lack of available memory                      *
*                                                                      *
* subroutines used :     gauss, mueller                                *
*                                                                      *
* Author                 Uli Eggermann, 8.24.1990                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int  j,  k, *perm;
  REAL f, *a, **gl, *rs, **lumat;
  void *vmblock;

  #define zi  rs                      /* rename often used vector rs  */
  #define eta rs
/***********************************************************************
* allocate buffers, free if unsuccessful                               *
***********************************************************************/
  vmblock = vminit();
  a     = (REAL *)vmalloc(vmblock, VEKTOR,  n + 1, 0);
  rs    = (REAL *)vmalloc(vmblock, VEKTOR,  n,     0);
  perm  = (int *)vmalloc(vmblock, VVEKTOR,  n,     sizeof(*perm));
  gl    = (REAL **)vmalloc(vmblock, MATRIX, n,     n);
  lumat = (REAL **)vmalloc(vmblock, MATRIX, n,     n);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 3;
  }
/**********************************************************************
* Prepare linear system from the integral data and solve using gauss. *
* The solution vector contains the coefficients of the Lagrange       *
* polynomial with leading coefficient 1.                              *
**********************************************************************/
  for (j=0; j<n; j++)
  {
    for (k=0; k<n;k++)
      gl [j] [k] = - Integrale [j + k];
    rs [j] = Integrale [j + n];
  }

  j = gauss (0, n, gl, lumat, perm, rs, a, &k);              /* gauss */

  if (j == 0)
  {
    vmfree(vmblock);
    return 1;
  }
/***********************************************************************
* find roots of polynomial                                             *
***********************************************************************/
  a [n] = ONE;                            /* leading coefficient = 1  */
  k   = mueller (n, a, 1, StStelle, zi);                    /* MULLER */
  if (k != 0)
  {
    vmfree(vmblock);
    return 2;
  }
/***********************************************************************
* Compute weights                                                      *
***********************************************************************/
  for (j=0; j<n; j++)
  {
    f = eta [n-1] = ONE;
    for (k=n-2; k>=0; k--)
    {
      eta [k] = a [k+1] + eta [k+1] * StStelle [j];
      f = eta [k] + f * StStelle [j];
    }
    Gewicht [j] = ZERO;
    for (k=0; k<n; k++)
      Gewicht [j] += eta [k] / f * Integrale [k];
  }
  vmfree(vmblock);
  return 0;
} /* OrtogP */

/* --------------------------- END quanec.c ------------------------- */
