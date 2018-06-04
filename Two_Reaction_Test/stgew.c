#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE stgew.c ------------------------- */

/***********************************************************************
* Nodes and weights for various quadrature formulas:                   *
*                                                                      *
*  - AdaQuaStGew    :  Adaptive Gauss quadrature                       *
*  - ClenCurtStGew  :  Clenshaw-Curtis quadrature formulas             *
***********************************************************************/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>
#include <gax.h>

/***********************************************************************
* Constants and variables:                                             *
***********************************************************************/

static int  ITERMAX = 300;              /* Maximal number iterations  */
static REAL ABSERR  = 0.0;              /* allowable absolute error   */
#define RELERR  (1000.0 * MACH_EPS)     /* relative error allowed     */
#define FKTERR  (100.0 * MACH_EPS)      /* Max. error in function val.*/

/***********************************************************************
* Prototypes:                                                          *
***********************************************************************/

static int pegasus2 (REAL  fkt (REAL, REAL *, int),
/*.IX{pegasus2}*/
                     REAL* _x1,
                     REAL* _x2,
                     REAL* _f2,
                     int*  _it,
                     REAL* LegCoef,
                     int   LegGrad)
/***********************************************************************
* Compute a root of the continuous function fkt(x) via pegasus method. *
* Here the signs of fkt (x1) and fkt (x2) must be different.           *
*                                                                      *
* Input parameters:                                                    *
*   REAL   fkt (REAL)   Name of function                               *
*   REAL   *x1, *x2     inclusion interval with fkt(x1) * fkt(x2) <= 0 *
*                                                                      *
* Output parameters:                                                   *
*   REAL   *x2          approximate root                               *
*   REAL   *f2          functional value at  x2                        *
*   int    *iter        number of iterations                           *
*                                                                      *
* Return value :                                                       *
*   -1:  no inclusion: fkt(x2) * fkt(x1) > 0                           *
*    0:  root has been found with  abs(f2) < FKTERR                    *
*    1:  stopped with | xold - xnew |  <  ABSERR + xnew * RELERR,      *
*         => check function value  f2 !                                *
*    2:  max number of iterations reached                              *
*                                                                      *
* Constants used :                                                     *
*    ABSERR, RELERR, MACH_EPS, EPSROOT, ITERMAX                        *
***********************************************************************/
{
  #define x1 (*_x1) /* serves for readability: avoids use of *s       */
  #define x2 (*_x2)
  #define f2 (*_f2)
  #define it (*_it)

  #define samesign(a,b)  (((a)>=ZERO && (b)>=ZERO) || \
                          ((a)<ZERO && (b)<ZERO))

  REAL f1, x3, f3;
  int  res = 2;

  f1 = fkt (x1, LegCoef, LegGrad);         /* Functional value f (x1) */
  f2 = fkt (x2, LegCoef, LegGrad);         /* Functional value f (x2) */

  if (f1 == 0.0)                                    /* x1 is solution */
  {
    x2 = x1;
    f2 = ZERO;
    return 0;
  }
  if (f2 == ZERO)                                   /* x2 is solution */
    return 0;

  if (samesign (f1, f2))                 /* no inclusion => error    */
    return -1;

  for (it = 0; it < ITERMAX; it++)             /* Pegasus iterations  */
  {
    x3  = x2 - (f2 / (f2 - f1)) * (x2 - x1);          /* new x-value  */
    f3  = fkt (x3, LegCoef, LegGrad);                 /* new f-value  */

    if (samesign (f2, f3))                /* for equal signs:         */
      f1 *= f2 / (f2 + f3);               /*  decrease  f1 modulus    */
    else
    {                           /* otherwise new inclusion interval  */
      x1 = x2;
      f1 = f2;
    }
    x2 = x3;                                   /* re-store            */
    f2 = f3;
    if (FABS (f2) < FKTERR)                    /* root found !        */
    {
      res = 0;
      break;
    }
    if (FABS (x2 - x1) <= FABS (x2) * RELERR + ABSERR)
    {
      res = 1;                           /* Stop: step size too small */
      break;
    }
  }
  if (FABS (f1) < FABS (f2))            /* used smaller modulus value */
  {                                     /* for  f2                    */
    x2 = x1;
    f2 = f1;
  }
  #undef x1
  #undef x2
  #undef f2
  #undef it
  #undef samesign
  return res;
}
/*.BA*/

int LegendreCoeff (REAL* s,
/*.IX{LegendreCoeff}*/
                   int   Grad,
                   REAL  *Leg_r)
/***********************************************************************
* Compute the coefficients of the desired Legendre polynomial of degree*
* Grad. Its coefficients are entered into the vector s.                *
*                                                                      *
* The recursion formula for the coeffcients is given by :              *
*                                                                      *
*            (n+1) P   (x) = (2n+1) x P (x) - n P   (x)                *
*                   n+1                n         n-1                   *
*                                                                      *
* globale variables used : REAL* Leq_r, int GewichtsBerechnung         *
***********************************************************************/
/*.BE*/
{
  int i, n;
  REAL f1, f2, *Leg_q;
  int  Leg_r_NULL;
  void *vmblock;
  Leg_r_NULL = (Leg_r == NULL);
  if (Grad==0)                                   /* Trivial case zero */
  {
    *s = ONE;
    return 0;
  }
  if (Grad==1)                                   /* Special case one */
  {
    *s++ = ZERO;
    *s   = ONE;
    if (! Leg_r_NULL)               /* caller provides storage and    */
      Leg_r[0] = ONE;               /* thus needs Leg_r values?       */
    return 0;
  }

  /*********************************************************************
  * allocate space for the aux vectors Leg_q and Leg_r                 *
  *********************************************************************/
  vmblock = vminit();
  Leg_q = (REAL *)vmalloc(vmblock, VEKTOR, Grad - 1, 0);
  if (Leg_r_NULL)
    Leg_r = (REAL *)vmalloc(vmblock, VEKTOR, Grad,     0);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 1;
  }
  Leg_q[0] = ONE;
  Leg_r[0] = ZERO;
  Leg_r[1] = ONE;

  for (n = 1; n < Grad; n++)                             /* Recursion */
  {
    f1 =    n    / (REAL) (n+1);
    f2 = (2*n+1) / (REAL) (n+1);

    for (i = n+1; i >= 0; i-- ) s[i]  = 0.0;
    for (i = n+1; i >= 1; i-=2) s[i]  = f2 * Leg_r[i-1];
    for (i = n-1; i >= 0; i-=2) s[i] -= f1 * Leg_q[ i ];

    if (n != Grad-1)
    {
      for (i = 0; i <= n;   i++) Leg_q[i] = Leg_r[i];
      for (i = 0; i <= n+1; i++) Leg_r[i] = s[i];
    }
  }
  vmfree(vmblock);
  return 0;
}
/*.BA*/

REAL LegPolWert (REAL x, REAL *LegCoef, int LegGrad)
/*.IX{LegPolWert}*/
/***********************************************************************
* Computes a functional value for a Legendre polynomial.               *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*   REAL x                 value for evaluation                        *
*                                                                      *
* Return value :           value of Legendre polynomial at x           *
*                                                                      *
* global variables:        REAL* LegCoef                               *
*                          int   LegGrad                               *
*                                                                      *
* used in:                 LegendreNullst (for  pegasus2)              *
*                          AdaQuaStGew                                 *
.BA*)
***********************************************************************/
/*.BE*/
{
  REAL f;
  int  n = LegGrad;
  for (f = LegCoef[n--]; n >= 0; n--)
       f = f * x + LegCoef[n];
  return f;
}

#if 0
void PrintCoeff (REAL* C, int Grad)
{
  int n;
  printf ("LegCoeff %d:", Grad);
  for (n = Grad; n >= 0; n--) printf (" %+8.5"LZP"g", C[n]);
  printf ("\n");
}
#endif
/*.BA*/

int LegendreNullst (REAL* Root,
/*.IX{LegendreNullst}*/
                    int   Grad,
                    REAL  *LegCoef,    /* for array of cofficients    */
                    REAL  *Leg_r)
/***********************************************************************
* Computes the roots of a Legendre polynomial.                         *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*   int Grad             degree of Legendre polynomial                 *
*                                                                      *
* Output parameter:                                                    *
*                                                                      *
*   REAL Root []         all roots of the Legendre polynomial          *
*                                                                      *
* global variables:      REAL* LegCoef, int LegGrad                    *
*                                                                      *
* subroutines used:      LegendreCoeff, LegPolWert, pegasus2           *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  REAL z = 3.141592654 / (Grad-.5);
  REAL x1, x2, f2, a = ONE, b;                     /* Pegasus bounds  */
  int  i, n, m = Grad / 2;
  int  LegGrad;                /* for length of coefficienten array   */
  void *vmblock;

  LegGrad = Grad;                                 /* global  variable */

  vmblock = vminit();
  LegCoef = (REAL *)vmalloc(vmblock, VEKTOR, Grad + 1, 0);
  if (! vmcomplete(vmblock))
     return 1;

  if (LegendreCoeff (LegCoef, Grad, Leg_r) != 0)     /* coefficients  */
  {
     vmfree(vmblock);
     return 1;
  }

#if 0                                     /* 1, if coefficients should*/
    PrintCoeff (LegCoef, Grad);           /*    be put out for test   */
#endif                                    /*    purposes, else  0     */

  for (n = 0; n < m; n++)
  {
    x1 = a;
    x2 = b = .5 * (COS ((n + .5) * z) + COS ((n + 1) * z));
    if (pegasus2 (LegPolWert, &x1, &x2, &f2, &i, LegCoef, LegGrad) != 0)
    {
      vmfree(vmblock);
      return 2;
    }
    Root[n] = - x2; a = b;
  }
  if (Grad % 2 != 0) Root [m++] = ZERO;
  for (n = m; n < Grad; n++)
    Root [n] = - Root [Grad - n - 1];

  vmfree(vmblock);
  return 0;
}
/*.BA*/

int AdaQuaStGew (
/*.IX{AdaQuaStGew}*/
                 int Grad,
                 REAL* StStelle,
                 int AuchGewichte,
                 REAL* Gewicht
                )
/***********************************************************************
* Computes the nodes and weights, if desired, of an adaptive quadrature*
* formula of given degree.                                             *
* The interval is [-1, 1].                                             *
.BE*)
*                                                                      *
* Parameter:                                                           *
*   int  Grad           degree of quadrature formula                   *
*   int  AuchGewichte   Flag for computing weights also                *
*   REAL StStelle []    computed nodes                                 *
*   REAL Gewicht  []    computed nodes (if AuchGewichte = 1)           *
*                                                                      *
* Return value :                                                       *
*   0:     all ok                                                      *
*   1:     Legendre polynomial roots not found                         *
*   2:     improper degree                                             *
*                                                                      *
* subroutines used:    LegendreNullst, LegPolWert                      *
*                                                                      *
* global variables:    REAL* Leg_r, int GewichtsBerechnung             *
*                                                                      *
* Author:              Uli Eggermann,  1.13.1991                       *
.BA*)
***********************************************************************/
/*.BE*/
{
   int i; REAL d, x0;
   REAL* LegCoef;              /* for coefficient array               */
   int   LegGrad;              /* for length of coefficient array     */
   REAL* Leg_r;                /* globally used dummy array           */
   void  *vmblock;

   if (Grad < 0) return (2);       /* no roots for  Grad < Null       */
   if (Grad == 0)
   {
     *StStelle = ZERO;
     if (AuchGewichte) *Gewicht = TWO;
     return (0);
   }

   vmblock = vminit();
   LegCoef = (REAL *)vmalloc(vmblock, VEKTOR, Grad + 1, 0);
   if (! vmcomplete(vmblock))
     return 1;

   if (AuchGewichte)
   {
     Leg_r = (REAL *)vmalloc(vmblock, VEKTOR, Grad, 0);
     if (! vmcomplete(vmblock))
     {
       vmfree(vmblock);
       return 1;
     }
   }
   else
     Leg_r = NULL;
   if (LegendreNullst (StStelle, Grad,           /* the nodes are the */
                       LegCoef, Leg_r) != 0)
   {                              /* roots of the Legendre polynomial */
      vmfree(vmblock);
      return 1;
   }

   if (AuchGewichte)                      /* if desired, weights also */
   {
     LegGrad = Grad-1;
     for (i = 0; i < Grad; i++)
     {
       d = LegPolWert (x0 = StStelle[i], Leg_r, LegGrad) * Grad;
       Gewicht[i] = 2 * (1 - x0 * x0) / (d * d);
     }
   }
   vmfree(vmblock);
   return (0);
}
/*.BA*/

int ClenCurtStGew (
/*.IX{ClenCurtStGew}*/
                   int n,
                   REAL* StStelle,
                   REAL* Gewicht
                  )
/***********************************************************************
* Computes the nodes and weights of a Clenshaw-Curtis quadrature       *
* formula of local error order  n+3 for the reference interval [-1,1]. *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   int n                n+1 is the number of nodes and weights:       *
*                        n > 1, n odd                                  *
*   REAL   StStelle []   computed nodes                                *
*   REAL   Gewicht  []   computed weights                              *
*                                                                      *
* Return value :                                                       *
*   0:                   all ok                                        *
*   1:                   n too small or odd                            *
*                                                                      *
* Author:                Uli Eggermann,  9.24.1990                     *
.BA*)
***********************************************************************/
/*.BE*/

{
  int k, j, m;
  REAL p, f, g, h, i = ONE, d;
  if (n < 2 || n % 2 != 0) return 1;              /* n - Test         */
  m = n / 2;                                      /* m = n/2          */
  d = (REAL)(n * n - 1);                          /* d = n*n+1        */
  h = 2. / (n * d);                               /* h = 2/(n(n*n+1)) */
  f = 4. / n;                                     /* f = 4/n          */
  p = f * ATAN (1.);                              /* p = pi/n         */

  Gewicht  [0] =       Gewicht  [n] = 1 / d;       /*  left end point */
  StStelle [0] = -ONE; StStelle [n] = ONE;         /* right end point */

  for (k = 1; k <= m; k++, i = -i)                /* Interval (-1,0]  */
  {
    for (g = ZERO, j = 1; j < m; j++)
         g += COS (2 * j * k * p) / (4 * j * j - 1);
    Gewicht  [k] = h * (d + i) - f * g;
    StStelle [k] = - COS (k * p);
  }
  StStelle [m] = ZERO;                          /* center of interval */

  for (k = m+1; k < n; k++)                       /* Interval (0, 1)  */
  {
    Gewicht  [k] =   Gewicht  [n - k];
    StStelle [k] = - StStelle [n - k];
  }
  return 0;
} /* ClenCurtStGew */

/* --------------------------- END stgew.c -------------------------- */
