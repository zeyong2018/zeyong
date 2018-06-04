#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE rk_fehl.c ------------------------ */

/***********************************************************************
*                                                                      *
* Solve a first order ordinary differential equation sysatem using the *
* -------------------------------------------------------------------- *
* Runge-Kutta-Fehlberg method [O(h^6)] with estimates for the local    *
* -----------------------------------------------------------------    *
* error and step size control                                          *
* ---------------------------                                          *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Richard Reuter (FORTRAN)                       *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 1.10.1992                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  MACH_EPS, FABS, max, SQRT, REAL, ONE, */
                       /*       dglsysfnk, min, boolean, FALSE, TRUE, */
                       /*       ZERO, TWO, THREE, EIGHT, HALF, NINE   */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR                                */
#include <rk_fehl.h>   /*  for  rk_fehl                               */



/* ------------------------------------------------------------------ */

#define MAX_FCT  10000 /* maximal number of evaluations for the right */
                       /* right hand side of the differential equation*/



/* ------------------------------------------------------------------ */
/*.BA*/

int rk_fehl   /* Runge-Kutta-Fehlberg method for 1st order DE systems */
/*.IX{rk\unt fehl}*/
           (
            REAL      *x,        /* initial / final x-value ..........*/
            REAL      xend,      /* desired x-value ..................*/
            int       n,         /* number of DEs ....................*/
            REAL      y[],       /* initial/final y-value ............*/
            dglsysfnk dgl,       /* righ thand side for the system ...*/
            REAL      *h,        /* initial/final step size ..........*/
            REAL      hmax,      /* maximal step size ................*/
            REAL      epsabs,    /* absolute error bound .............*/
            REAL      epsrel     /* relative error bound .............*/
           )                     /* error code .......................*/

/***********************************************************************
* Solve a system of first order ordinary differential equations using  *
* the Runge-Kutta-Fehlberg method of order O(h^6) with estimation of   *
* the local error and step size control.                               *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x       initial x-value                                              *
* xend    desired x-value  (xend < x is allowed)                       *
* n       number of DEs in system                                      *
* y       intial y-value at x                                          *
* dgl     pointer to function that evaluates the right hand side of    *
*         the system  y' = f(x,y)                                      *
* h       starting step size; if h is put in unreasonably, we adjust h *
*         internally; h may be negative if xend < x.                   *
* hmax    upper limit for the magnitude of used step sizes (hmax > 0.0)*
* epsabs  \ local error bounds for relative/absolute errors wrt.       *
* epsrel  / current step size. If in each component of the computed    *
*           solution  y[i]                                             *
*               |estimate for local errro|  <=                         *
*                          |h| * (epsrel * |y[i]| + epsabs),           *
*           then we accept the current step.                           *
*           If epsabs = 0  =>  Test for relative error only;           *
*           if epsrel = 0  =>  Test for absolute error only.           *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x       final x-value reached (usually x = xend)                     *
* y       final y-value of solution at x                               *
* h       final step size used                                         *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok  (solution found at  xend)                               *
* = 1: after MAX_FCT calls of dgl() we stopped without reaching xend.  *
*      One can call rk_fehl() again with unaltered parameters.         *
* = 2: improper input, i. e.                                           *
*          epsabs < 0             or                                   *
*          epsrel < 0             or                                   *
*          epsabs + epsrel = 0    or                                   *
*          hmax <= 0.                                                  *
* = 3: the optimal step size cannot be represented in the computer     *
* = 4: The test of the local error can't be executed because epsabs    *
*      and a component of the computed approximate solution are zero   *
*      at the same time. Therefore a pure test of the local error is   *
*      impossible for this component.                                  *
* = 5: lack of available memory                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* MAX_FCT, REAL, MACH_EPS, boolean, FALSE, TRUE, min, max, dglsysfnk,  *
* FABS, vminit, vmalloc, vmcomplete, vmfree, VEKTOR, SQRT, ONE, ZERO,  *
* TWO, THREE, EIGHT, HALF, NINE                                        *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL    *yt,       /* aux vectors for Runge-Kutta step              */
          *t,
          *r,        /* Vector with estimates for the local error     */
          *k1,       /* \                                             */
          *k2,       /*  \                                            */
          *k3,       /*   \  aux vectors for                          */
          *k4,       /*   /  Runge-Kutta step                         */
          *k5,       /*  /                                            */
          *k6,       /* /                                             */
          da,        /* Lenght of integration interval                */
          xt,        /* aux again for R-K                             */
          hmx,       /* maximal step size, at most |da|               */
          hf = ZERO, /* aux variable for step size                    */
          quot,      /* Measure for accuracy of integration step      */
          et,        /* epsrel * |yt[i]| + epsabs                     */
          tr;        /* aux variable for QUOT                         */
  int     i,         /* Loop variable                                 */
          aufrufe;   /* number of calls of  dgl                       */
  boolean repeat,    /* Flag indicating step with same h              */
          ende;      /* Flag indicating method can be stopped         */
  void    *vmblock;  /* List of dynamic allocations                   */


  aufrufe = 0;
  repeat  = FALSE;
  ende    = FALSE;
  da      = xend - *x;


  if (epsrel < ZERO || epsabs < ZERO ||           /* check input      */
      epsrel + epsabs == ZERO || hmax <= ZERO)
    return 2;
  if (FABS(da) <= (REAL)13.0 * MACH_EPS * max(FABS(*x), FABS(xend)))
    return 3;


  vmblock = vminit();                 /* initialize storage           */
  yt = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  t  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  r  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  k1 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  k2 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  k3 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  k4 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0),
  k5 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0),
  k6 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))                      /* lack of memory ? */
  {
    vmfree(vmblock);         /* free allocations and report error     */
    return 4;
  }


  hmx = min(hmax, FABS(da));
  if (FABS(*h) <= (REAL)13.0 * MACH_EPS * FABS(*x))
    *h = hmx;

  for ( ; ; )
  {
    if (! repeat)                              /* limit h to hmx      */
    {
      *h = min(FABS(*h), hmx);        /* and choose h so as to reach  */
      *h = (da > ZERO) ? *h : -*h;    /* the final point xend, if     */
      if (FABS(xend - *x) <= (REAL)1.25 * FABS(*h))       /* possible */
      {
        hf   = *h;            /* If ende has beenset and h = xend - x */
        ende = TRUE;          /* is admissable we can stop after the  */
        *h   = xend - *x;     /* next integration                     */
      }

      (*dgl)(*x, y, k1);                  /* perform one integration  */
      aufrufe++;
    }

    xt = *h * (REAL)0.25;
    for (i = 0; i < n; i++)
      yt[i] = y[i] + xt * k1[i];
    xt += *x;
    (*dgl)(xt, yt, k2);

    for (i = 0; i < n; i++)
      yt[i] = y[i] + (*h) * (k1[i] * (THREE / (REAL)32.0) +
                             k2[i] * (NINE / (REAL)32.0));

    xt = *x + *h * (REAL)0.375;
    (*dgl)(xt, yt, k3);
    for (i = 0; i < n; i++)
      yt[i] = y[i] + *h * (k1[i] * ((REAL)1932.0 / (REAL)2197.0) -
                           k2[i] * ((REAL)7200.0 / (REAL)2197.0) +
                           k3[i] * ((REAL)7296.0 / (REAL)2197.0)
                          );

    xt = *x + *h * ((REAL)12.0 / (REAL)13.0);
    (*dgl)(xt, yt, k4);
    for (i = 0; i < n; i++)
      yt[i] = y[i] + (*h) * (k1[i] * ((REAL)439.0 / (REAL)216.0) -
                             k2[i] * EIGHT +
                             k3[i] * ((REAL)3680.0 / (REAL)513.0) -
                             k4[i] * ((REAL)845.0 / (REAL)4104.0)
                            );

    xt = *x + *h;
    (*dgl)(xt, yt, k5);
    for (i = 0; i < n; i++)
      yt[i] = y[i] + *h * (-k1[i] * (EIGHT / (REAL)27.0) +
                            k2[i] * TWO -
                            k3[i] * ((REAL)3544.0 / (REAL)2565.0) +
                            k4[i] * ((REAL)1859.0/ (REAL)4104.0) -
                            k5[i] * ((REAL)11.0 / (REAL)40.0)
                          );

    xt = *x + HALF * *h;
    (*dgl)(xt, yt, k6);
    for (i = 0; i < n; i++)
    {
      t[i] = k1[i] * ((REAL)25.0 / (REAL)216.0) +
             k3[i] * ((REAL)1408.0 / (REAL)2565.0) +
             k4[i] * ((REAL)2197.0 / (REAL)4104.0) -
             k5[i] * (REAL)0.2;
      yt[i] = y[i] + *h * t[i];
    }

    /* yt is the preliminary result of the integration step.          */
    /* Next we compute r, the estimate for the local error, wrt. the  */
    /* current step size.                                             */

    for (i = 0; i < n; i++)
      r[i] = k1[i] / (REAL)360.0 -
             k3[i] * ((REAL)128.0 / (REAL)4275.0) -
             k4[i] * ((REAL)2197.0 / (REAL)75240.0) +
             k5[i] / (REAL)50.0 +
             k6[i] * (TWO / (REAL)55.0);

    for (i = 0, quot = ZERO; i < n; i++)          /* check accuracy   */
    {
      et = epsrel * FABS(yt[i]) + epsabs;
      if (et == ZERO)                         /* If the approximate   */
        return 5;                             /* solution vanishes a  */
      tr   = FABS(r[i]) / et;                 /* pure check of the    */
      quot = max(quot, tr);                   /* relative error is    */
    }                                         /* not useful.          */

    if (quot <= ONE)                          /* result acceptable ?  */
    {
      for (i = 0; i < n; i++)
        y[i] = yt[i];
      *x += *h;
      if (ende)                     /* stop at xend, if reached       */
      {
        *h = hf;
        vmfree(vmblock);
        return 0;
      }
      quot = max(quot, (REAL)0.00065336);        /* prepare next step */
    }

    quot = min(quot, (REAL)4096.0);      /* increase h maximally by a */
    *h = (REAL)0.8 * *h /                /* factor of 5, or decrease  */
         SQRT(SQRT(quot));               /* h at most by the factor 10*/

    if (FABS(*h) <= (REAL)13.0 * MACH_EPS * FABS(*x))
    {
      vmfree(vmblock);
      return 3;
    }
    aufrufe += 5;

    if (aufrufe >= MAX_FCT)
    {
      vmfree(vmblock);
      return 1;
    }

    if (quot > ONE)                  /* repeat step with smaller h ? */
    {
      repeat = TRUE;
      ende   = FALSE;
    }
    else                                       /* step acceptable ?  */
      repeat = FALSE;
  }
}

/* ------------------------- END rk_fehl.c -------------------------- */
