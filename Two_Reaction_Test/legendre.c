#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE legendre.c ------------------------ */

/***********************************************************************
*                                                                      *
* Compute roots of a Legendre polynomial                               *
* --------------------------------------                               *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Hermann-Josef Rheinbach (FORTRAN)              *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               equivalent C, Pascal, QuickBASIC and FORTRAN   *
*                       codes                                          *
* Date:                 1.13.1992                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  PI, MACH_EPS, horner, FABS, REAL,    */
                        /*       COS, ZERO, ONE, HALF                 */
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vmfree, vminit, */
                        /*       VEKTOR                               */
#include <legendre.h>   /*  for  legendre                             */



static REAL gxpega    /* compute one root of a polynomial ............*/
/*.IX{gxpega}*/
                  (
                   REAL a,
                   REAL b,
                   REAL *c,
                   int  n
                  )

/***********************************************************************
* Compute the roots of a polynomial of degree n, with coefficients in  *
* the vector c. The root lies inside the interval [a,b]. The method is *
* an adapted pegasus method.                                           *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* a  \ end point of inclusion interval for root                        *
* b  /                                                                 *
* c  [0..n] coefficient vector for the polynomial:                     *
*    c[i] is the coefficient of x^i.                                   *
* n  degree of poynomial                                               *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* polynomial root in [a,b]                                             *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, horner, MACH_EPS, FABS, ZERO                                   *
***********************************************************************/

{
  REAL xdiff,
       x1,
       x2,
       x3,
       f1,
       f2,
       f3,
       s12;


  x1 = a;
  x2 = b;
  f1 = horner(n, c, x1);
  f2 = horner(n, c, x2);
  xdiff = x2 - x1;

  for ( ; ; )
  {
    s12 = xdiff / (f2 - f1);
    x3  = x2 - f2 * s12;
    f3  = horner(n, c, x3);
    if (f2 * f3 < ZERO)
      x1 = x2,
      f1 = f2;
    else
      f1 *= f2 / (f2 + f3);

    x2 = x3;
    f2 = f3;
    if (FABS(f2) < MACH_EPS * (REAL)100.0)
      if (FABS(f1) < FABS(f2))
        return x1;
      else
        return x2;
    xdiff = x2 - x1;
    if (FABS(xdiff) < MACH_EPS * (REAL)100.0)
      if (FABS(f1) < FABS(f2))
        return x1;
      else
        return x2;
  }
}

/* ------------------------------------------------------------------ */
/*.BA*/

int legendre    /* Compute all roots of a Legendre polynomial ........*/
/*.IX{legendre}*/
            (
             int  grad,                      /* polynomial degree ....*/
             REAL *alpha                     /* roots ................*/
            )                                /* error code ...........*/

/***********************************************************************
* Compute all roots of a Legendre poynomial of degree grad.            *
* All lie inside the interval  [-1,1].                                 *
.BE*)
*                                                                      *
* Input parameter:                                                     *
* ================                                                     *
* grad  degree of Legendre polynomial                                  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* alpha  [0..grad-1] root vector                                       *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok                                                          *
* = 3: lack of memory                                                  *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* gxpega, REAL, PI, COS, vminit, vmalloc, vmcomplete, vmfree, VEKTOR,  *
* ZERO, ONE, HALF                                                      *
.BA*)
***********************************************************************/
/*.BE*/

{
  float zw;
  REAL  *p_old,          /* [0..grad] vectors ........................*/
        *p_midd,         /*                                           */
        *c,              /*                                           */
        xk,
        xk_plus,
        xk_inv,
        xfa,
        grenza,
        grenzb;
  int   i,               /* Loop variables ...........................*/
        k,               /*                                           */
        gradhalb;        /* greatest integer in  grad / 2 ............*/
  void  *vmblock;        /* List of dynamic allocations ..............*/


  if (grad == 1)
  {
    alpha[0] = ZERO;
    return 0;
  }

  vmblock = vminit();                 /* initialize storage           */

  p_old  = (REAL *)vmalloc(vmblock, VEKTOR, grad + 1, 0);
  p_midd = (REAL *)vmalloc(vmblock, VEKTOR, grad + 1, 0);
  c      = (REAL *)vmalloc(vmblock, VEKTOR, grad + 1, 0);

  if (! vmcomplete(vmblock))                     /* lack of storage ? */
  {
    vmfree(vmblock);         /* free storage and report error         */
    return 3;
  }

  p_old[0]  = ONE;                     /* compute coefficients of     */
  p_midd[0] = ZERO;                    /* Legendre polynomials        */
  p_midd[1] = ONE;
  xk        = ZERO;

  for (k = 1; k < grad; k++)
  {
    xk++;
    xk_plus = xk + ONE;
    xk_inv  = ONE / xk_plus;
    xfa     = (xk + xk_plus) * xk_inv;
    for (i = 0; i <= k; i++)
      c[i + 1] = p_midd[i] * xfa;
    c[0] = ZERO;
    xfa = xk * xk_inv;
    for (i = 0; i <= k - 1; i++)
    {
      c[i]     -= p_old[i] * xfa;
      p_old[i] =  p_midd[i];
    }
    p_old[k] = p_midd[k];
    for (i = 0; i <= k + 1; i++)
      p_midd[i] = c[i];
  }


  grenza   = ONE;                          /* compute roots symmetric */
  zw       = (float)(PI / (grad - HALF));  /* to origin               */
  gradhalb = grad / 2;

  for (i = 1; i <= gradhalb; i++)        /* compute negative roots    */
  {
    grenzb  = HALF * (COS(((REAL)i - HALF) * zw) + COS((REAL)i * zw));
    alpha[i - 1] = -gxpega(grenza, grenzb, c, grad);
    grenza  = grenzb;
  }

  for (i = 0, k = grad - 1; i < gradhalb; i++, k--) /* find positive  */
    alpha[k] = -alpha[i];                           /* roots by       */
                                                    /* reflection     */
  if (grad % 2)                             /* odd degree polynomial? */
    alpha[k] = ZERO;                                /* zero is a root */

  vmfree(vmblock);

  return 0;
}

/* -------------------------- END legendre.c ------------------------ */
