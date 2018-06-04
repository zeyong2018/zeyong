#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 3.3.3}{Bauhuber's Method}{Bauhuber's Method}*/

/*.BE*/
/* ------------------------ MODULE fbauhube.c ----------------------- */

#include <basis.h>
#include <u_proto.h>

#define ITERMAX 1000               /* Maximal number of function      */
                                   /* evaluations per root            */
#define EPS \
  (REAL)((REAL)64.0 * MACH_EPS)    /* Accuracy for functional value   */
#define BETA (REAL)((REAL)8.0 * EPS)
#define QR (REAL)0.1               /* Real and imaginary parts of the */
#define QI (REAL)0.9               /* spiralization constant          */


static void scpoly
                   (int      n,      /* length of vector .............*/
                    REAL     ar[],   /* Real part of the vector ......*/
                    REAL     ai[],   /* Imaginary part of the vector .*/
                    REAL   * scal);  /* Scaling factor ...............*/

static int bauroot
                   (int    n,        /* largest degree ...............*/
                    int    iu,       /* lowest degree ................*/
                    REAL   ar[],     /* Real parts of the coefficients*/
                    REAL   ai[],     /* Imaginary parts, coefficients */
                    REAL   * x0r,    /* Real part of the root ........*/
                    REAL   * x0i);   /* Imaginary part of the root ...*/

static void chorner
                    (int      n,     /* highest degree in polynomial .*/
                     int      iu,    /* lowest degree in polynomial ..*/
                     REAL     ar[],  /* Real parts of coefficients ...*/
                     REAL     ai[],  /* Imaginary parts, coefficients */
                     REAL     xr,    /* Real part of x ...............*/
                     REAL     xi,    /* Imaginary part of x ..........*/
                     REAL   * pr,    /* Real part of function value ..*/
                     REAL   * pi,    /* Imaginary part of function v. */
                     REAL   * p1r,   /* Real part of first derivative */
                     REAL   * p1i,   /* Imaginary part, first deriv. .*/
                     REAL   * p2r,   /* Real part, second derivative .*/
                     REAL   * p2i,   /* Imaginary part, second deriv. */
                     REAL   * rf1);  /* Error estimate for 1st deriv. */

static void polydiv
                    (int    n,       /* maximal degree ...............*/
                     int    iu,      /* minimal degree ...............*/
                     REAL   ar[],    /* Real parts of coefficients ...*/
                     REAL   ai[],    /* Imaginary parts, coefficients */
                     REAL   x0r,     /* Real part of x ...............*/
                     REAL   x0i);    /* Imaginary part of x ..........*/

/*.BA*/

int bauhub              /* Bauhuber's method for complex polynomials .*/
/*.IX{bauhub}*/
           (
            int    real,          /* Are the coefficients real ? .....*/
            int    scale,         /* Scaling ? .......................*/
            int     n,            /* degree of polynomial ............*/
            REAL  ar[],           /* Real parts of coefficients ......*/
            REAL  ai[],           /* Imaginary parts, coefficients ...*/
            REAL  rootr[],        /* Real parts of roots .............*/
            REAL  rooti[],        /* Imaginary parts of roots ........*/
            REAL  absf[]          /* Absolute value of function values*/
           )
/*====================================================================*
 *                                                                    *
 *  bauhub uses Bauhuber's Method to find all real or complex roots   *
 *  of a polynomial of degree n :                                     *
 *                                               n-1           n      *
 *      P(x) = a[0] + a[1] * x + ... + a[n-1] * x    + a[n] * x ,     *
 *                                                                    *
 *  where a[i], i=0, ..., n, can be complex.                          *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Application:                                                     *
 *   ===========                                                      *
 *      Find roots of arbitrary polynomials with complex coefficients.*
 *      If the polynomial roots are ill-condi=ditioned, i.e., if small*
 *      changes in the coefficients lead to large changes in the roots*
 *      the polynomial should not be scaled. Otherwise scaling helps  *
 *      with stability and performance.                               *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      real     int real;                                            *
 *        = 0    Polynomial coefficients are complex                  *
 *       != 0    Polynomial coefficients are real                     *
 *      scale    int scale;                                           *
 *        = 0    no scaling                                           *
 *       != 0    Scaling, see polysc()                                *
 *      n        int n;                                               *
 *               degree of the polynomial ( >= 1 )                    *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the polynomial           *
 *               coefficients   ( ar[0], ..., ar[n] )                 *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      rootr    REAL   rootr[];   (Vector of length  n+1 )           *
 *               rootr[0],..,rootr[n-1] are the real parts of the n   *
 *               roots                                                *
 *      rooti    REAL   rooti[];   (Vector of length n+1 )            *
 *               rooti[0],..,rooti[n-1] are the imaginary parts       *
 *               of the roots                                         *
 *      absf     REAL   absf[];                                       *
 *               absf[0],..,absf[n-1] are the magnitudes of the       *
 *               polynomial values at the computed roots              *
 *                                                                    *
 *   Return value:                                                    *
 *   ============                                                     *
 *      = 0      all is ok                                            *
 *      = 1      n < 1 or invalid input parameter                     *
 *      = 2      ar[n] = 0.0 and ai[n] = 0.0                          *
 *      = 3      Iteration maximum ITERMAX exceeded                   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *     bauroot():   determines one root of the polynomial             *
 *     scpoly():    Scales the polynomial                             *
 *     chorner():   Evaluates the polynomial                          *
 *     polydiv():   factors off a root                                *
 *     comabs():    magnitude of a complex number                     *
 *                                                                    *
 *====================================================================*/
{
  int    i, res;
  REAL x0r, x0i, tempr, tempi, t1, t2, t3, t4, t5;
  REAL scalefak = ONE;

  if (n < 1) return (1);

  if (ar == NULL || ai == NULL) return (1);
  if (rootr == NULL || rooti == NULL) return (1);

  if (ar[n] == ZERO && ai[n] == ZERO)  /* Leading coefficient must    */
    return (2);                        /* differ from zero            */

  for (i = 0; i <= n; i++)        /* store given coefficients in root */
  {
    rootr[i] = ar[i];
    rooti[i] = ai[i];
    if (i < n) absf[i] = ZERO;
  }

  scalefak = ONE;
  if (scale)                          /* Scale polynomial, if desired */
    scpoly (n, rootr, rooti, &scalefak);

  x0r = ZERO;                                       /* Starting value */
  x0i = ZERO;

  for (i = 0; i < n; i++)
  {                                           /* compute the ith root */
    res = bauroot (n, i, rootr, rooti, &x0r, &x0i);

    rootr[i] = scalefak * x0r;                          /* store root */
    rooti[i] = scalefak * x0i;

    if (res) return (3);                /* Iteration maximum reached? */

    /* Polynomial value of input polynomial at (rootr[i], rooti[i])   */

    chorner (n, 0, ar, ai, rootr[i], rooti[i],
             &tempr, &tempi, &t1, &t2, &t3, &t4, &t5);

    absf[i] = comabs (tempr, tempi);            /* store error        */

    polydiv (n, i, rootr, rooti, x0r, x0i);     /* reduce degree      */

    if (real)                           /* New starting value         */
      x0i = -x0i;                       /* depending on real x..      */
    else
      x0r = x0i = ZERO;
  }

  return (0);
}


static void scpoly
/*.IX{scpoly}*/
                   (int       n,     /* length of vector .............*/
                    REAL    ar[],    /* Real part of vector ..........*/
                    REAL    ai[],    /* Imaginary part of vector .....*/
                    REAL *  scal     /* Scaling factor ...............*/
                   )
/*====================================================================*
 *                                                                    *
 *  scalpoly skales the polynomial P :                                *
 *                                               n-1           n      *
 *      P(x) = a[0] + a[1] * x + ... + a[n-1] * x    + a[n] * x ,     *
 *                                                                    *
 *  where all a[i], i=0, ..., n, can be complex.                      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               degree of the polynomila ( >= 1 )                    *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients         *
 *               a[0],..,a[n]                                         *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients         *
 *               a[0],..,a[n] of the scaled polynomial.               *
 *      scal     REAL   *scal;                                        *
 *               Scaling factor                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *      comabs():  Modulus of a complex number                        *
 *                                                                    *
 *   From the  C - library : pow()                                    *
 *                                                                    *
 *   Macros:     max                                                  *
 *====================================================================*/
{
  REAL p, pot;
  int    i;

  *scal = ZERO;
                             /* scal =                                */
  p = comabs (ar[n], ai[n]); /*               a[i]  1/(n-i)           */
  for (i = 0; i < n; i++)    /*    max{ cabs( ---- )       ,i=0..n-1} */
                             /*               a[n]                    */
    if (ar[i] != ZERO || ai[i] != ZERO)
    {
      ai[i] /= p;
      ar[i] /= p;

      pot = POW (comabs (ar[i],ai[i]), 1.0/(n-i));
      *scal = max (*scal, pot);
    }

  ar[n] /= p;                /* Absolute value of a[n] = 1            */
  ai[n] /= p;

  if (*scal == ZERO) *scal = ONE;

  for (p = ONE, i = n-1; i >= 0; i--)
  {
    p *= *scal;                /*                    n-i              */
    ar[i] /= p;                /* a[i] = a[i] / (scal    ), i=0..n-1  */
    ai[i] /= p;                /*                                     */
  }
}


static void chorner        /* Complex Horner scheme ..................*/
/*.IX{chorner}*/
                    (
                     int     n,           /* largest pol. degree .....*/
                     int     iu,          /* lowest pol. degree ......*/
                     REAL    ar[],        /* Real parts of coeff. ....*/
                     REAL    ai[],        /* Imaginary parts of coeff.*/
                     REAL    xr,          /* Real part of x ..........*/
                     REAL    xi,          /* Imaginary part of x .....*/
                     REAL *  pr,          /* Real part of function v. */
                     REAL *  pi,          /* Imaginary part, function */
                     REAL *  p1r,         /* Real part first deriv. ..*/
                     REAL *  p1i,         /* Imaginary part first der.*/
                     REAL *  p2r,         /* Real part 2nd derivative */
                     REAL *  p2i,         /* Imaginary part 2nd der. .*/
                     REAL *  rf1          /* Error estimate 1st der. .*/
                    )
/*====================================================================*
 *                                                                    *
 *  Horner scheme for polynomial with complex coefficients.           *
 *  We compute :                                                      *
 *    1. Polynomial value of the polynomial P (complex) of degree     *
 *       n - iu,                                                      *
 *    2. value of first derivative at x,                              *
 *    3. value of 2nd derivative at x,                                *
 *    4. an error estimate for the first derivative.                  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               Maximal degree of the polynomial ( >= 1 )            *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients of the  *
 *               polynomial with the coefficients a[iu], ..., a[n]    *
 *      x0r,x0i  REAL   x0r, x0i;                                     *
 *               Real and imaginary parts of the point of evaluation  *
 *                                                                    *
 *   Ausgabeparameter:                                                *
 *   ================                                                 *
 *      pr, pi   REAL   *pr, *pi;                                     *
 *               Real and imaginary part of the polynomial            *
 *      p1r, p1i REAL   *p1r, *p1i;                                   *
 *               Real and imaginary parts of the 1st derivative there *
 *      p2r, p2i REAL   *p2r, *p2i;                                   *
 *               Real and imaginary parts of the 2nd derivative       *
 *      rf1      REAL   *rf1;                                         *
 *               Error estimate for the first derivative              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *     comabs():    modulus of a complex number                       *
 *     pow ():      Power function                                    *
 *====================================================================*
 *                                                                    *
 *   Constant used:  EPS                                              *
 *   =============                                                    *
 *                                                                    *
 *====================================================================*/
{
  register i, j;
  int      i1;
  REAL   temp, temp1;

  *p2r = ar[n];
  *p2i = ai[n];

  *pr = *p1r = *p2r;
  *pi = *p1i = *p2i;

  *rf1 = comabs (*pr, *pi);
  i1 = n - iu;

  for (j = n - iu, i = n - 1; i >= iu; i--, j--)
  {
    temp = *pr;                        /* Polynomial value (pr,pi)    */
    *pr = *pr * xr - *pi * xi + ar[i];
    *pi = *pi * xr + temp * xi + ai[i];
    if ( i == iu ) break;

    temp = *p1r;                       /* 1st derivative (p1r,p1i)    */
    *p1r = *p1r * xr - *p1i * xi;
    *p1i = *p1i * xr + temp * xi;

    temp = comabs (*p1r, *p1i);        /* Error estimate for the 1st  */
    *p1r += *pr;                       /* derivative of P             */
    *p1i += *pi;

    temp1 = comabs (*pr, *pi);
    temp = max (temp, temp1);
    if ( temp > *rf1 )
    {
      *rf1 = temp;
      i1 = j - 1;
    }

    if (i - iu <= 1) continue;

    temp = *p2r;                       /* 2nd derivative (p2r,p2i)    */
    *p2r = *p2r * xr - *p2i * xi + *p1r;
    *p2i = *p2i * xr + temp * xi + *p1i;
  }

  temp = comabs (xr, xi);
  if ( temp != ZERO )
    *rf1 *= POW (temp, (REAL) i1) * (i1 + 1);
  else
    *rf1 = comabs (*p1r, *p1i);

  *rf1 *= EPS;

  *p2r += *p2r;
  *p2i += *p2i;
}


static int bauroot
/*.IX{bauroot}*/
                   (
                    int     n,        /* highest degree ..............*/
                    int     iu,       /* lowest degree ...............*/
                    REAL    ar[],     /* Real parts of coefficients ..*/
                    REAL    ai[],     /* Imaginary parts, coefficients*/
                    REAL *  x0r,      /* Real part of the root .......*/
                    REAL *  x0i       /* Imaginary part of the root ..*/
                   )
/*====================================================================*
 *                                                                    *
 *  bauroot computes one root of the polynomial P of degree n-iu:     *
 *                                                 n-iu               *
 *      P(x) = a[iu] + a[iu+1] * x + ... + a[n] * x                   *
 *                                                                    *
 *  with complex  coefficients a[i], i=iu, ..., n.                    *
 *  This program uses Newton's method on the function P(x) / P'(x).   *
 *  The iteration is stabilized using spiralization and extrapolation.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               Maximal degree of the polynomial ( >= 1 )            *
 *      iu       int iu;                                              *
 *               Index for the constant term of the polynomial,       *
 *               n-iu is the degree of the polynomial with            *
 *               coefficients a[iu], ..., a[n]                        *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients         *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      x0r,x0i  REAL   *x0r, x0i;                                    *
 *               Real and imaginary part of the computed root         *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all is ok                                            *
 *      = 1      Division by zero                                     *
 *      = 2      Iteration number ITERMAX exceeeded                   *
 *      = 3      Improper input                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *     chorner():   Complex Horner scheme                             *
 *     comabs():    magnitude of a complex number                     *
 *     comdiv():    Complex division                                  *
 *     quadsolv():  solves quadratic equations                        *
 *                                                                    *
 *   Constants used : TRUE, FALSE, ITERMAX,                           *
 *   ===============  QR, QI, MACH_EPS, EPS, EPSROOT, BETA            *
 *                                                                    *
 *====================================================================*/
{
  int  rc, result = 2, endit = FALSE,
       iter = 0, i = 0;

  REAL xoldr = ZERO, xoldi = ZERO,
       xnewr, xnewi, h, h1, h2, h3, h4, dzmax, dzmin,
       dxr = ZERO, dxi = ZERO, tempr, tempi,
       abs_pold, abs_pnew, abs_p1new, temp, ss, u, v, bdze = ZERO,
       pr, pi, p1r, p1i, p2r, p2i;

  if (n < 1) return (3);

  if (n - iu == 1)                          /* Polynomial of degree 1 */
  {
    quadsolv (ZERO, ZERO, ar[n], ai[n], ar[n-1], ai[n-1], x0r, x0i);
    return (0);
  }

  if (n - iu == 2)                          /* Polynomial of degree 2 */
  {
    quadsolv (ar[n],ai[n], ar[n-1],ai[n-1], ar[n-2],ai[n-2], x0r,x0i);
    return (0);
  }

  xnewr = *x0r;  xnewi = *x0i;
  endit = FALSE;

  chorner (n, iu, ar, ai, xnewr, xnewi,     /* Evaluate polynomial   */
           &pr, &pi, &p1r, &p1i, &p2r, &p2i, &ss);
  iter++;

  abs_pnew = comabs (pr, pi);
  if (abs_pnew < EPS) return (0);           /* Starting value is a   */
                                            /* good approximation    */
  abs_pold = abs_pnew;
  dzmin = BETA * (EPSROOT + comabs (xnewr, xnewi));

  while ( iter < ITERMAX )   /* Bauhuber-Iteration */
  {
    abs_p1new = comabs (p1r, p1i);

    if (abs_pnew > abs_pold)                /* Spiralization step     */
    {
      i = 0;                                /* dx = dx * q            */
      iter++;
      temp = dxr;

      dxr = QR * dxr - QI * dxi;
      dxi = QR * dxi + QI * temp;
    }
    else
    {
      dzmax = ONE + comabs (xnewr, xnewi);
      h1 = p1r * p1r - p1i * p1i - pr * p2r + pi * p2i;
      h2 = TWO * p1r * p1i - pr * p2i - pi * p2r;
      if (    abs_p1new > (REAL)10.0 * ss
           && comabs (h1, h2) > (REAL)100.0 * ss * ss)
      {
        i++;
        if ( i > 2 ) i = 2;
        tempr = pr * p1r - pi * p1i;
        tempi = pr * p1i + pi * p1r;

        rc = comdiv (-tempr, -tempi, h1, h2, &dxr, &dxi);
        if (rc != 0) return (1);

        if ( comabs (dxr, dxi) > dzmax )
        {
          temp = dzmax / comabs (dxr, dxi);       /* Newton step     */
          dxr *= temp;
          dxi *= temp;
          i = 0;
        }

        if (   i == 2
            && comabs (dxr, dxi) < dzmin / EPSROOT
            && comabs (dxr, dxi) > ZERO)
        {
          i = 0;                            /* Extrapolation step    */
          rc = comdiv (xnewr - xoldr, xnewi - xoldi,
                                   dxr, dxi, &h3, &h4);
          if (rc != 0) return (1);

          h3 += ONE;
          h1 = h3 * h3 - h4 * h4;
          h2 = TWO * h3 * h4;
          rc = comdiv (dxr, dxi, h1, h2, &h3, &h4);
          if (rc != 0) return (1);

          if ( comabs (h3, h4) < (REAL)50.0 * dzmin )
          {
            dxr += h3;
            dxi += h4;
          }
        }

        xoldr = xnewr;
        xoldi = xnewi;
        abs_pold = abs_pnew;
      }
      else
      {
        i = 0;                             /* Close to a saddle point */
        h = dzmax / abs_pnew;
        dxr = h * pr;
        dxi = h * pi;

        xoldr = xnewr;
        xoldi = xnewi;
        abs_pold = abs_pnew;

        do
        {
          chorner (n, iu, ar, ai, xnewr+dxr, xnewi+dxi,
                         &u, &v, &h, &h1, &h2, &h3, &h4);
          iter++;

          dxr += dxr;
          dxi += dxi;                        /* dx = dx * 2.0         */
        }
        while ( ABS (comabs (u,v) / abs_pnew - ONE) < EPSROOT );
      }
    }

    if (endit)
    {
      if ( comabs (dxr, dxi) < (REAL)0.1 * bdze )
      {
        xnewr += dxr;  xnewi += dxi;
      }

      result = 0;
      break;                                     /* stop iteration   */
    }
    else
    {
      xnewr = xoldr + dxr;
      xnewi = xoldi + dxi;
      dzmin = BETA * (EPSROOT + comabs (xnewr, xnewi));
      chorner (n, iu, ar, ai, xnewr, xnewi,
               &pr, &pi, &p1r, &p1i, &p2r, &p2i, &ss);
      iter++;
      abs_pnew = comabs ( pr, pi);

      if (abs_pnew == ZERO)
      {
        result = 0;
        break;
      }

      if (comabs (dxr, dxi) < dzmin || abs_pnew < EPS)
      {
        endit = TRUE;
        bdze = comabs (dxr, dxi);
      }
    }
  } /* End Bauhuber iteration */

  *x0r = xnewr;
  *x0i = xnewi;

  return (result);
}


static void polydiv
/*.IX{polydiv}*/
                    (
                     int   n,        /* highest degree ...............*/
                     int   iu,       /* lowest degree ................*/
                     REAL  ar[],     /* Real parts of coefficients ...*/
                     REAL  ai[],     /* Imaginary parts, coefficients */
                     REAL  x0r,      /* Real part of x ...............*/
                     REAL  x0i       /* Imaginary part of x ..........*/
                    )
/*====================================================================*
 *                                                                    *
 *  polydiv computes the coefficients of the polynomial Q, with       *
 *  P(x) = Q(x) * (x - x0), where x0 is a computed root of P.         *
 *  Both P and Q and x0 may be complex.                               *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               Highest degree of the polynomial ( >= 1 )            *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficienten of P   *
 *               of degree n-iu, with a[iu], ..., a[n]                *
 *      x0r,x0i  REAL   x0r, x0i;                                     *
 *               Real and imaginary parts of the root x0              *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients         *
 *               ar[iu+1],..,ar[n] of the remainder polynomial Q      *
 *                                                                    *
 *====================================================================*/
{
  register i;
  REAL     temp;

  for (i = n - 1; i > iu; i--)
  {
    temp = ar[i+1];
    ar[i] += temp * x0r - ai[i+1] * x0i;
    ai[i] += ai[i+1] * x0r + temp * x0i;
  }
}

/* ------------------------- END fbauhube.c ------------------------- */
