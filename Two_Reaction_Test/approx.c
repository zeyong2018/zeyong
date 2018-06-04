#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE approx.c ------------------------ */

/***********************************************************************
*                                                                      *
* Linear approximation                                                 *
* --------------------                                                 *
*                                                                      *
* exported functions:                                                  *
*   - gfq():        find the coefficients of a approximating           *
*                   polynomial using the discrete least square method  *
*                   of Gauss.                                          *
*   - pol_appr():   discrete linear least squares via orthogonal       *
*                   polynomials                                        *
*   - opolwert():   evaluate the polynomial from  pol_appr()           *
*   - opolkoeff():  compute the standard coefficients from the         *
*                   orthogonal coefficients returned from pol_appr()   *
*   - lin_happr():  linear least squares via Householder               *
*                   transformation                                     *
*   - lin_hwert():  evaluate function from lin_happr()                 *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Juergen Dietel, Computer center, RWTH Aachen   *
* Sources:              existing C, Pascal and QuickBASIC programs     *
* Date:                 9.24.1992-2.27.1997                            *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  POW, NULL, MACH_EPS, sqr, PI, REAL,   */
                       /*       ansatzfnk, approxfnk, ableitfnk, ONE, */
                       /*       ZERO, TWO, HALF, THREE, FABS, SQRT,   */
                       /*       SIN, COS                              */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR, MATRIX                        */
#include <u_proto.h>   /*  for  choly, fhouse                         */
#include <approx.h>    /*  for  gfq, pol_appr, opolwert, opolkoeff,   */
                       /*       lin_happr, lin_hwert                  */



/*.BA*/
/*.FE{C 8.1.3.1}
     {Normal Equations for Discrete Linear Least Squares}
     {Normal Equations for Discrete Linear Least Squares}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int gfq            /* Normal equations for least square approximation */
/*.IX{gfq}*/
    (
     int  n,                 /* degree of approximation polynomial ...*/
     int  m,                 /* number of nodes - 1 ..................*/
     REAL x[],               /* nodes: x-values ......................*/
     REAL y[],               /*        y-values ......................*/
     REAL w[],               /* weights ..............................*/
     REAL c[]                /* coeffic. of approximating polynomial  */
    )                        /* error code ...........................*/

/***********************************************************************
* Find the coefficients of a linear least square approximating         *
* polynomial of degree n using the discrete least square method of     *
* Gauss.                                                               *
.BE*)
* The algorithm uses the Gaussian normal equations, which are solved   *
* using the Cholesky method: The coefficienten c[j] of the vector      *
*                                                                      *
*        ( w[0] * (c[0] + c[1] * x[0]^1 + ... + c[n] * x[0]^n) )       *
*        ( w[1] * (c[0] + c[1] * x[1]^1 + ... + c[n] * x[1]^n) )       *
*        ( ...                            ...                  )       *
*        ( w[m] * (c[0] + c[1] * x[m]^1 + ... + c[n] * x[m]^n) )       *
*                                                                      *
* are to be computed so that it approximates the given right hand side *
*                                                                      *
*                        ( w[0] * y[0] )                               *
*                        ( w[1] * y[1] )                               *
*                        ( ...         )                               *
*                        ( w[m] * y[m] )                               *
*                                                                      *
* optimally in the euclidean norm.                                     *
* c is the solution of the linear system                               *
*                        a * c  =  b              (Normal equations)   *
* where                                                                *
*      a[i][j] = w[0] * x[0]^(j+i)    + ... + w[m] * x[m]^(j+i),       *
*      b[i]    = w[0] * y[0] * x[0]^i + ... + w[m] * y[m] * x[m]^i     *
* for i,j=0, ..., n.                                                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* m: Index of terminal node                                            *
* n: degree of algebraic polynomial                                    *
* x: [0..m] vector of x-values                                         *
* y: [0..m] vector of y-values                                         *
* w: [0..m] vector of weights                                          *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* c: [0..n] vector of polynomial coefficients                          *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: error in input: n or m too small                                  *
* 2: error in  choly()                                                 *
* 3: lack of memory                                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* vminit, vmalloc, vmcomplete, vmfree, VEKTOR, MATRIX, REAL, choly,    *
* ZERO                                                                 *
.BA*)
***********************************************************************/
/*.BE*/

{
  void *vmblock;      /* List of dynamically allocated vectors and    */
                      /* matrices                                     */
  REAL **a,           /* [0..n,0..n] matrix of the linear system      */
       *b,            /* [0..n] vector of the right hand side         */
       summand;       /* aux variable                                 */
  int  i, j;          /* Loop variables                               */


  if (n < 0 || m < n)                           /* n or m too small ? */
    return 1;


  /* ------- allocate space for matrix a and right hand side b ------ */

  vmblock = vminit();                             /* allocate storage */
  a = (REAL **)vmalloc(vmblock, MATRIX, n + 1, n + 1);
  b = (REAL *) vmalloc(vmblock, VEKTOR, n + 1, 0);
  if (! vmcomplete(vmblock))                      /* any problems ?   */
    return 3;


  for (i = 0; i <= n; i++)            /* compute first column, last   */
    a[i][0] = a[n][i] = b[i] = ZERO;  /* row and the right hand side  */
                                      /* of the system                */
  for (j = 0; j <= m; j++)
  {
    summand = w[j];
    for (i = 0; i < n; i++)
      a[i][0] += summand,
      b[i]    += summand * y[j],
      summand *= x[j];
    a[n][0] += summand;
    b[n]    += summand * y[j];
    for (i = 1; i <= n; i++)
      summand *= x[j],
      a[n][i] += summand;
  }


  for (i = n - 1; i >= 0; i--)    /* complete matrix by replicating   */
    for (j = 0; j < n; j++)       /* row i in row i + 1 in shifted    */
      a[i][j + 1] = a[i + 1][j];  /* form                             */


  i = choly(0, n + 1, a, b, c);                /* solve linear system */


  vmfree(vmblock);                            /* free storage buffers */

  if (i != 0)                     /* error solving linear system ?    */
    return 2;
  else                            /* no error ?                       */
    return 0;
}
/*.BA*/



/*.FE{C 8.1.3.2}
     {Discrete Least Squares via Orthogonal Polynomials}
     {Discrete Least Squares via Algebraic Polynomials and
      Orthogonal Polynomials}*/

/*.BE*/
/* ------------------------------------------------------------------ */

static REAL qwert
/*.IX{qwert}*/
    (
     int  k,
     REAL x,
     REAL b[],
     REAL d[]
    )

/***********************************************************************
* Evaluate the orthogonal polynomial Q of degree k at x                *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* k: degree of polynomial                                              *
* x: x-value                                                           *
* b: [1..k] aux vectors \  for computing the orthogonal polynomials    *
* d: [2..k]-Hilfsvektor /  up to degree k                              *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* Qk(x)                                                                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ONE                                                            *
***********************************************************************/

{
  int  i;            /* Loop index                                    */
  REAL qi = ZERO,    /* Function value of Qi at x                     */
       qi1,          /* Function value of Q(i-1) at x                 */
       qi2;          /* Function value of Q(i-2) at x                 */


  switch (k)
  {
    case 0:  return ONE;

    case 1:  return x - b[1];

    default: qi2 = ONE;
             qi1 = x - b[1];
             for (i = 2; i <= k; i++, qi2 = qi1, qi1 = qi)
               qi = (x - b[i]) * qi1 - d[i] * qi2;
             return qi;
  }
}



#define  QQ  1     /* Flag for computing the inner products (Qk,Qk)   */
#define  YQ  2     /* Flag for computing the inner products (y,Qk)    */
#define XQQ  3     /* Flag for computing the inner products (x*Qk,Qk) */


/* ------------------------------------------------------------------ */

static REAL skal
/*.IX{skal}*/
    (
     int  k,
     int  m,
     int  qq_yq_xqq,
     REAL x[],
     REAL y[],
     REAL w[],
     REAL b[],
     REAL d[]
    )

/***********************************************************************
* Compute the inner products for  pol_appr(). Depending on the flag    *
* qq_yq_xqq one of the inner products  (Qk,Qk), (f,Qk) or is computed. *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* k:         degree of the orthogonal polynomial used in the inner     *
*            product                                                   *
* m:         number of nodes - 1                                       *
* qq_yq_xqq: Flag designating type of inner product to be computed     *
*            = QQ:  (Qk,Qk)                                            *
*            = YQ:  (y,Qk)                                             *
*            = XQQ: (x*Qk,Qk)                                          *
* x:         [0..m] vector of x-values                                 *
* y:         [0..m] vector of y-values                                 *
* w:         [0..m] weight vector                                      *
* b: [1..k] aux vectors \  for orthogonal polynomials of lower degree  *
* d: [2..k]             /                                              *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* inner product                                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* qwert, QQ, YQ, XQQ, REAL, sqr, ZERO                                  *
***********************************************************************/

{
  int  i;                        /* Loop index                        */
  REAL summe;                    /* aux sum                           */


  summe = ZERO;

  switch (qq_yq_xqq)
  {
    case QQ:
      for (i = 0; i <= m; i++)
        summe += w[i] * sqr(qwert(k, x[i], b, d));
      break;
    case YQ:
      for (i = 0; i <= m; i++)
        summe += w[i] * y[i] * qwert(k, x[i], b, d);
      break;
    case XQQ:
      for (i = 0; i <= m; i++)
        summe += w[i] * x[i] * sqr(qwert(k, x[i], b, d));
  }

  return summe;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int pol_appr       /* discr. lin. least squares via orthog. polynom.  */
/*.IX{pol\unt appr}*/
    (
     int  n,                 /* degree of least square appr. polynom. */
     int  m,                 /* number of nodes - 1 ..................*/
     REAL x[],               /* nodes: x-values ......................*/
     REAL y[],               /*        y-values ......................*/
     REAL w[],               /* weights ..............................*/
     REAL c[],               /* coefficients of optimal polynomial ...*/
     REAL b[],               /* aux variables for lower degree .......*/
     REAL d[]                /* orthogonal polynomials ...............*/
    )                        /* error code ...........................*/

/***********************************************************************
* Compute the coefficients c, b and d for a least squares approximating*
* polynomial P of degree n for m+1 nodes (x[i],y[i]), i=0, ..., m and  *
* m+1 weights w[i], using discrete orthogonal polynomials.             *
.BE*)
*                                                                      *
* The optimal solution P has the form:                                 *
*                                                                      *
*     P(x) = c[0] * Q (x) + c[1] * Q (x) + .. + c[n] * Q (x),          *
*                    0              1                   n              *
*                                                                      *
* where                                                                *
*                                                                      *
*    Q (x) = 1;    Q (x) = x - b[1];                                   *
*     0             1                                                  *
*                                                                      *
*    Q (x) = (x - b[j]) * Q   (x) - d[j] * Q   (x),       j=2,..,n     *
*     j                    j-1              j-2                        *
*                                                                      *
* and                                                                  *
*    b[j] = (x * Q   (x), Q   (x)) / ( Q   (x), Q   (x)), j=1,..,n     *
*                 j-1      j-1          j-1      j-1                   *
*                                                                      *
*    d[j] = (Q   (x), Q   (x)) / (Q   (x), Q   (x)),      j=2,..,n     *
*             j-1      j-1         j-2      j-2                        *
*                                                                      *
*    c[j] = (y, Q (x)) / (Q (x), Q (x)),                  j=0,..,n .   *
*                j         j      j                                    *
*                                                                      *
* Here (u, v) denotes the inner product of u and v, weighted by w, i.e.
*       (u, v)  =  w[0] * u[0] * v[0] + ... + w[m] * u[m] * v[m].      *
*                                                                      *
* This polynomial is evaluated as :                                    *
*                                                                      *
*                            P(x) = s                                  *
*                                    0                                 *
*                                                                      *
* where                                                                *
*                                                                      *
*     s  = c[n],         s    = c[n-1] + s  * (x - b[n]),              *
*      n                  n-1             n                            *
*                                                                      *
*                                                                      *
*     s  = c[j] + s    * (x - b[j+1]) - s    * d[j+2],  j = n-2,..,0   *
*      j           j+1                   j+2                           *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n: degree of polynomial                                              *
* m: number of nodes - 1                                               *
* x: [0..m] vector of x-values                                         *
* y: [0..m] vector of y-values                                         *
* w: [0..m] weight vector                                              *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* c: [0..n] coefficient vector for the orthogonal polynomial system    *
* b: [1..n] aux vectors \  to evaluate the polynomial Qk we need       *
* d: [2..n]             /  b[1],...,b[k] and d[2],...,d[k].            *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: n > m or m < 1 or n < 0                                           *
* 2: m = n (Interpolation), and two x-values are identical             *
* 3: Not all weights are positive                                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* skal, QQ, YQ, XQQ, REAL, ZERO                                        *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  k, j;                             /* Loop indeces              */
  REAL qk_qk,                            /* inner product (Qk,Qk)     */
       qk_qk_m1,                         /* inner product (Qk-1,Qk-1) */
       qk_qk_m2,                         /* inner product (Qk-2,Qk-2) */
       summe1,                           /* inner product (xQ0,Q0)    */
       summe2;                           /* inner product (Q0,Q0)     */


  if (m < n || m < 1 || n < 0)         /* too few nodes or polynomial */
    return 1;                          /* degree erroneous ?          */

  if (m == n)                 /* For  m = n the x[k] must be distinct */
    for (k = 1; k <= m; k++)
      for (j = 0; j < k; j++)
        if (x[k] == x[j])
          return 2;

  for (k = 0; k <= n; k++)              /* one weight not positive ?  */
    if (w[k] <= ZERO)
      return 3;


  for (j = 0, summe1 = ZERO; j <= m; j++)
    summe1 += w[j] * x[j];
  for (j = 0, summe2 = ZERO; j <= m; j++)
    summe2 += w[j];
  b[1] = summe1 / summe2;

  qk_qk_m2 = skal(0, m, QQ, x, y, w, b, d);
  qk_qk_m1 = skal(1, m, QQ, x, y, w, b, d);

  c[0]     = skal(0, m, YQ, x, y, w, b, d) / qk_qk_m2;
  c[1]     = skal(1, m, YQ, x, y, w, b, d) / qk_qk_m1;

  for (k = 2; k <= n; k++, qk_qk_m2 = qk_qk_m1, qk_qk_m1 = qk_qk)
    b[k]  = skal(k - 1, m, XQQ, x, y, w, b, d) / qk_qk_m1,
    d[k]  = qk_qk_m1 / qk_qk_m2,
    qk_qk = skal(k,     m, QQ,  x, y, w, b, d),
    c[k]  = skal(k,     m, YQ,  x, y, w, b, d) / qk_qk;


  return 0;
}


#undef  QQ
#undef  YQ
#undef XQQ



/* ------------------------------------------------------------------ */
/*.BA*/

REAL opolwert      /* Evaluate the polynomial from  pol_appr() .......*/
/*.IX{opolwert}*/
    (
     int  n,                 /* degree of polynomial .................*/
     REAL x,                 /* x-value ..............................*/
     REAL b[],               /* coeffici. for orthogonal polynomials  */
     REAL d[],
     REAL c[]                /* coefficients of optimal polynomial ...*/
    )                        /* value of polynomial at x .............*/

/***********************************************************************
* Evaluate the polynomial P from  pol_appr() at x.                     *
.BE*)
* Due to the two step recursion this code has the same problem as      *
* qwert(). Hence we again do not use recursion here.                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n: degree of P                                                       *
* x: x-value                                                           *
* b: [1..n] aux vectors \  to evaluate the orth. polynomials. for Qk   *
* d: [2..n]-Hilfsvektor /  we need b[1],...,b[k] and  d[2],...,d[k]    *
*                          (see pol_appr()).                           *
* c: [0..n] coefficient vector for the expansion of P wrt. the         *
*    orthogonal system of polynomials                                  *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* P(x)                                                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  k;            /* Loop index                                    */
  REAL sk = ZERO,    /* current value s[k] of the horner like scheme  */
       sk1,          /* s[k+1]                                        */
       sk2;          /* s[k+2]                                        */


  switch (n)
  {
    case 0:  return c[0];

    case 1:  return c[0] + c[1] * (x - b[1]);

    default: sk2 = c[n];
             sk1 = c[n - 1] + c[n] * (x - b[n]);
             for (k = n - 2; k >= 0; k--, sk2 = sk1, sk1 = sk)
               sk = c[k] + sk1 * (x - b[k + 1]) - sk2 * d[k + 2];
             return sk;
  }
}



/* ------------------------------------------------------------------ */

/*.BA*/

int opolkoeff      /* orthogonal coeffic.  -->  standard coeffic. ....*/
/*.IX{opolkoeff}*/
    (
     int  n,                 /* degree of polynomial .................*/
     REAL b[],               /* coefficients for computation .........*/
     REAL d[],               /* of orthogonal polynomials ............*/
     REAL c[],               /* coeff. for c[0]*Q_0 + ... + c[n]*Q_n  */
     REAL a[]                /* coeff. for a[0]*x^0 + ... + a[n]*x^n  */
    )                        /* error code ...........................*/

/***********************************************************************
* Diskrete linear approximation using orthogonal polynomials:          *
* Calculation of standard coefficients from orthog. coefficients.      *
.BE*)
* This function calculates the standard coefficients in vector a from  *
* the coefficients of the approximating polynomial given in the three  *
* vectors b, d and c produced by pol_appr():                           *
*                                                                      *
*                                   1                    n             *
*         P(x)  =  a[0]  +  a[1] * x   +  ... +  a[n] * x .            *
*                                                                      *
* While c contains the coefficients for the dev. of the discrete orth. *
* polynomials, the vectors b and d are used to calculate the orth.     *
* polynomials Q , Q , ..., Q  recursively (of course here we do NOT    *
*              0   1        n                                          *
* work it out recursively, but in normal iterations, because this is   *
* more efficient).                                                     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n  degree of polynomial P                                            *
* b  [1..n]-help vector \  for calculation of the orth. polyn., for    *
* d  [2..n]-help vector /  Q_i we use b[1],...,b[i] and d[2],...,d[i]  *
*                          (see pol_appr()).                           *
* c  [0..n]-vector with the coefficients of the development of         *
*    polynomial P for discrete orthogonal polynomials.                 *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* a  [0..n]-vector with the polynomial coefficients in the form        *
*                                    1                   n             *
*            P(x) = a[0]  +  a[1] * x   + ... +  a[n] * x              *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* Errorcode:                                                           *
*   0: no error                                                        *
*   1: n < 0                                                           *
*   3: not enough memory                                               *
*                                                                      *
* benutzte globale Namen:                                              *
* =======================                                              *
* REAL, vmalloc, VEKTOR, vmcomplete, vmfree                            *
.BA*)
***********************************************************************/
/*.BE*/

{

#define ZYKLISCH_TAUSCHEN(v0, v1, v2)                                  \
            /* exchange the three REAL-pointers v0,v1,v2 cyclically */ \
  {                                                                    \
    REAL *tmp;                                      /* help pointer */ \
    tmp = v0;                                                          \
    v0  = v1;                                                          \
    v1  = v2;                                                          \
    v2  = tmp;                                                         \
  }

#define RETURN(fehlercode)                                             \
              /* stop function including releasing allocated blocks */ \
  {                                                                    \
    vmfree (vmblock);                                                  \
    return (fehlercode);                                               \
  }

  REAL *ai2;         /* [0..i-2]-coefficient-vector of Q_{i-2}        */
                     /*          in standard form                     */
  REAL *ai1;         /* [0..i-1]-coefficient-vector of Q_{i-1}        */
                     /*          in standard form                     */
  REAL *ai;          /* [0..i]-coefficient-vector of Q_i              */
                     /*          in standard form                     */
  REAL bi;           /* temporary variable containing the value b[i]  */
  REAL di;           /* temporary variable containing the value d[i]  */
  REAL ci;           /* temporary variable containing the value c[i]  */
  int  i;            /* index variable for all orthogonal polynomials */
  int  j;            /* index variable for the a-coefficients of an   */
                     /* orthogonal polynomial                         */
  void *vmblock;     /* list of all dynamically allocated vectors     */


                               /* check input data                    */
  if (n < 0)                   /* illegal polynomial degree?          */
    return 1;

                               /* allocate memory for ai2, ai1, ai    */
  vmblock = vminit();          /* initialize memory block             */
  ai2 = (REAL *) vmalloc(vmblock, VEKTOR, n + 1, 0);
  ai1 = (REAL *) vmalloc(vmblock, VEKTOR, n + 1, 0);
  ai  = (REAL *) vmalloc(vmblock, VEKTOR, n + 1, 0);
  if (! vmcomplete(vmblock))   /* return, if one of the allocations   */
    RETURN(3);                 /*         has not been o.k.           */

  ai2[0] = 1.0;                /* set coefficients of Q_0             */
  a[0]   = c[0];               /* set a[0]                            */

  if (n != 0)                  /* polynomial degree not zero?         */
  {                            /* at least one more step              */

    ai1[1] = 1.0;              /* set coefficients of Q_1             */
    ai1[0] = -b[1];

    a[1]  = c[1];              /* set a[1], increase a[0]             */
    a[0] += c[1] * ai1[0];
  }

  /* calculate the a-coeff. of Q_i from the a-coeff. of Q_{i-1} and   */
  /* Q_{i-2} by using the recursive formula:                          */
  /*    Q (x) = (x - b[i]) * Q   (x) - d[i] * Q   (x),   i=2,..,n     */
  /*     i                    i-1              i-2                    */
  /* and finally set a[i] and increase the coefficients               */
  /* a[i-1],...,a[0] by the new values                                */

  for (i = 2; i <= n; i++)
  {
    /* calculate ai from ai1 and ai2, i.e.
    /* ai[i]...ai[0] from ai1[i-1]...ai1[0] and ai2[i-2]...ai2[0]     */
    bi        = b[i];
    di        = d[i];
    ai[i]     = 1.0;                    /* 1*a1[i-1] with a1[i-1] = 1 */
    ai[i - 1] = ai1[i - 2] - bi * ai1[i - 1];
    for (j = i - 2; j >= 1; j--)
      ai[j] = ai1[j - 1] - bi * ai1[j] - di * ai2[j];
    ai[0] = -bi * ai1[0] - di * ai2[0];

    /* set a[i], increase a[i-1],...,a[0] by new values               */
    a[i] = ci = c[i];
    for (j = i - 1; j >= 0; j--)
      a[j] += ci * ai[j];

    /* exchange the three REAL-pointers ai2, ai1 and ai cyclically    */
    /* for loop-continuation                                          */
    ZYKLISCH_TAUSCHEN (ai2, ai1, ai);
  }

  /* Now ai[0],...,ai[n] contain the standard coefficients of Q_n     */
  /* and a[0],...,a[n] contain the standard coefficients of the       */
  /* approximation polynomial.                                        */


  vmfree(vmblock);

  return 0;
#undef ZYKLISCH_TAUSCHEN
}
/*.BA*/
/*.BA*/



/*.FE{C 8.1.3.4}
     {Householder Transformation for Linear Least Squares}
     {Solving Linear Least Squares Problems using Householder
      Transformations}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int lin_happr      /* lin. least squares via Householder transform. ..*/
/*.IX{lin\unt happr}*/
    (
     int       m,            /* number of nodes ......................*/
     int       n,            /* number of functions -1 ...............*/
     REAL      x[],          /* nodes : x-values .....................*/
     REAL      y[],          /*         y-values .....................*/
     REAL      w[],          /* positive weights .....................*/
     ansatzfnk phi,          /* model functions ......................*/
     REAL      c[],          /* optimal coefficients .................*/
     REAL      *r            /* error of the opt. solution ...........*/
    )                        /* error code ...........................*/

/***********************************************************************
* Compute the coefficients c of the linear least square approximation  *
*                                                                      *
*          PHI(x) = c[0] * phi(0,x) + ... + c[n] * phi(n,x)            *
*                                                                      *
* for the given model functions  phi(0,.), ..., phi(n,.).              *
.BE*)
*                                                                      *
* The plot of the real valued function PHI shall approximate the m     *
* nodes  (x[i],y[i]), i=0, ..., m-1, with m >= n) optimally wrt. the   *
* mean root square error, i.e., the expression                         *
*    (*)    (y[0] - PHI(x[0]))^2 + ... + (y[m-1] - PHI(x[m-1]))^2      *
* is to be minimized.                                                  *
*                                                                      *
* The optimal coefficcients c[..] are determined via Householder       *
* transform. This avoids the generally ill-conditioned normal equations*
* for this overdetermined problem.                                     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* m:   number of nodes                                                 *
* n:   number of model functions - 1                                   *
* x:   [0..m-1] x-values                                               *
* y:   [0..m-1] y-values                                               *
* w:   [0..m-1] weight vector                                          *
* phi: points to a function which supplies the value of any of the n   *
*      model functions at x. phi has the form :                        *
*               REAL phi(int i, REAL x)                                *
*               {                                                      *
*                 return <value of the ith model function at x>        *
*               }                                                      *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* c: [0..n] coefficient vector for optimal solution                    *
* r: mean root square error of the approximation                       *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: all is ok                                                         *
* 1: m <= n or n < 1                                                   *
* 2: lack of memory                                                    *
* 3: The model functions phi(0.),...,phi(n,.) seem to be lin. dependent*
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* sqr, vminit, vmalloc, vmcomplete, vmfree, MATRIX, house, ansatzfnk,  *
* SQRT, REAL, ZERO                                                     *
.BA*)
***********************************************************************/
/*.BE*/

{
  void *vmblock;  /* List of dyn. allocated vectors and matrices      */
  int  i, j,      /* Loop indices                                     */
       res;       /* error code of  house()                           */
  REAL **mat,     /* [0..m-1,0..n] overdetermined system matrix       */
       *b,        /* [0..m-1] right hand side vector                  */
       we,        /* square root of w[i]                              */
       msqe;      /* mean square error                                */


  if (m <= n || n < 1)                /* too few nodes or functions ? */
    return 1;


  vmblock = vminit();                            /* allocate storage  */
  mat = (REAL **)vmalloc(vmblock,
                         MATRIX, m, n + 1);
  b   = (REAL *)vmalloc(vmblock, VEKTOR, m, 0);
  if (! vmcomplete(vmblock))                 /* memory insufficient ? */
    return 2;


  for (i = 0; i < m; i++)                    /*  form linear system : */
  {
    we = SQRT(w[i]);
    b[i] = we * y[i];             /* form right hand side and ...     */
    for (j = 0; j <= n; j++)      /* ... the system matrix, in which  */
      mat[i][j] = we * (*phi)(j, x[i]);   /* each row is multiplied by*/
  }                                       /* the respective weight    */


  res = house(m, n + 1, mat, b);  /* solve system via Householder     */
  for (i = 0; i <= n; i++)        /* copy solution in b to c          */
    c[i] = b[i];

  if (res != 0)                        /* Linear system numerically   */
  {                                    /* singular ?                  */
    vmfree(vmblock);
    return 3;
  }


  for (msqe = ZERO, i = n + 1; i < m; i++)/* compute mean square error*/
    msqe += sqr(b[i]);

  *r = SQRT(msqe);


  vmfree(vmblock);

  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

REAL lin_hwert     /* Evaluate function from lin_happr() .............*/
/*.IX{lin\unt hwert}*/
    (
     REAL      x0,           /* x-value ..............................*/
     int       n,            /* number of functions - 1 ..............*/
     ansatzfnk phi,          /* model functions ......................*/
     REAL      c[]           /* optimal coefficients .................*/
    )                        /* return value .........................*/

/***********************************************************************
* Evaluate the approximating function PHI from lin_happr() at x0:      *
*                                                                      *
*         PHI(x0) = c[0] * phi(0,x0) + ... + c[n] * phi(n,0).          *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0:  x-value                                                         *
* n:   number of functions - 1                                         *
* phi: points to a function which supplies the value of any of the n   *
*      model functions at x0. phi has the form :                       *
*               REAL phi(int i, REAL x)                                *
*               {                                                      *
*                 return <value of the ith model function at x>        *
*               }                                                      *
* c:   optimal coefficient vector from lin_happr()                     *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* PHI(x0)                                                              *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ansatzfnk, ZERO                                                *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  i;
  REAL funktionswert;

  for (funktionswert = ZERO, i = 0; i <= n; i++)
    funktionswert += c[i] * (*phi)(i, x0);

  return funktionswert;
}

/* --------------------------- END approx.c ------------------------- */
