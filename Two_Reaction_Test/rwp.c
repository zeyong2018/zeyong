#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* --------------------------- MODULE rwp.c ------------------------- */

/***********************************************************************
*                                                                      *
* Solve a two point boundary problem of first order with the shooting  *
* -------------------------------------------------------------------  *
* method                                                               *
* ------                                                               *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Klaus Niederdrenk (FORTRAN)                    *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 6.2.1992; 10.30.1995                           *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  max, MACH_EPS, FABS, POW, ONE,       */
                        /*       copy_vector, REAL, dglsysfnk, HALF,  */
                        /*       rndbedfnk, ZERO                      */
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vmfree, vminit, */
                        /*       VEKTOR, VVEKTOR, MATRIX              */
#include <u_proto.h>    /*  for  gauss                                */
#include <awp.h>        /*  for  awp                                  */
#include <ab_mou.h>     /*  for  prae_korr                            */
#include <einb_rk.h>    /*  for  einb_rk                              */
#include <bulirsch.h>   /*  for  bul_stoe                             */
#include <implruku.h>   /*  for  implruku                             */
#include <rwp.h>        /*  for  rwp                                  */


/* produce a vector of length `n' (and its name)                      */
#define zeig(v, n)                 \
  {                                \
    int i;                         \
    printf("%-9s", #v": ");        \
    for (i = 0; i < n; i++)        \
      printf("%16"LZP"g", v[i]);   \
    printf("\n");                  \
  }


/* ------------------------------------------------------------------ */
/*.BA*/

int rwp    /* Shooting method for boundary value problem of 1st order */
/*.IX{rwp}*/
       (
        REAL      a,           /* left end point .....................*/
        REAL      b,           /* right end point ....................*/
        REAL      h,           /* starting step size .................*/
        REAL      y_start[],   /* initial approximation or solution ..*/
                               /* initial value problem  y(a) ........*/
        int       n,           /* number of differntial equations ....*/
        dglsysfnk dgl,         /* right hand side for the system .....*/
        rndbedfnk rand,        /* Function for the boundary conditions*/
        int       awpnumm,     /* Number of desired IVP solver .......*/
        REAL      epsawp,      /* error bound for initial value       */
                               /* problem ............................*/
        REAL      epsrb,       /* error bound for boundary value      */
                               /* problem ............................*/
        long      fmax,        /* maximal number of calls of dgl() ...*/
        int       itmax,       /* maximal number of Newton iterations */
        int       *act_iter    /* actual number of Newton steps ......*/
       )                       /* error code .........................*/

/***********************************************************************
* This function solves the first order boundary value problem          *
*                                                                      *
*     y' = F(x,y),    a <= x <= b,    R(y(a),y(b)) = 0.                *
*                                                                      *
* It uses the shooting method starting with an approximate y_start for *
* the associated initial value y(a), from which it constructs an       *
* approximate solution using as initial value solver the function awp()*
* The nonlinear system arising in the shooting method is solved using  *
* Newton's method.                                                     *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* a        left end point                                              *
* b        right end point  (b > a)                                    *
* h        suitable step size for solving the associated initial value *
*          problem for the shooting method approximately               *
* y_start  [0..n-1] vector, initial approximation for y(a)             *
* n        number of differential equations                            *
* dgl      pointer to the function which evaluates the right hand side *
*          of the differential equation  y' = f(x,y)                   *
* rand     pointer to the function that evaluates the boundary cond.   *
* awpnumm  label for the desired IVP solver in the shooting method:    *
*          = 1: Runge-Kutta embedding formula of 4/5th order; England  *
*               formula from awp().                                    *
*          = 2: Predictor-corrector method of order 4 by Adams-        *
*               Bashforth-Moulton (from prae_korr())                   *
*          = 3: Runge-Kutta embedding formula of 7/8th order (from     *
*               einb_rk())                                             *
*          = 4: Extrapolation method of  Bulirsch-Stoer (from          *
*               bul_stoe())                                            *
*          = 5: implicit Runge-Kutta-Gauss method (from im[lruku())    *
* epsawp   desired accuracy for the solution of the associated initial *
*          value problem                                               *
* epsrb    desired accuracy for which the approximation y_start for    *
*          the initial value  y(a) of a solution of the boundary value *
*          problem should satisfy the boundary condition R.            *
* fmax     upper bound for number of calls of dgl() in the system of   *
*          differential equations when solving the associated initial  *
*          value problem                                               *
* itmax    upper bound for number of Newton iterations when solving    *
*          the nonlinear systems in the shooting method                *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y_start   [0..n-1] approximation for the initial value y(a) of a     *
*           solution  y of the boundary problem                        *
* act_iter  number of Newton iterations performed                      *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* error code                                                           *
* = 0: all is ok                                                       *
* = 1: at least one of the accuracy bounds  epsawp, epsrb too small    *
* = 2: b <= a                                                          *
* = 3: h <= 0                                                          *
* = 4: n <= 0                                                          *
* = 5: improper input for awpnumm                                      *
* = 6: maximal number of allowed function evaluations exceeded         *
* = 7: act_iter > itmax. number of allowd Newton steps exceeded without*
*      finding a suitable value for  y_start.                          *
* = 8: The Jacobi matrix for Newton iterations is singular.            *
* = 9: lack of memory space                                            *
* >  9: error in one of the IVP solvers at first node                  *
* > 19: error in one of the IVP solvers at second node                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, dglsysfnk, rndbedfnk, max, MACH_EPS, copy_vector, FABS,        *
* vminit, vmalloc, vmcomplete, vmfree, VEKTOR, VVEKTOR, MATRIX, awp,   *
* gauss, POW, ZERO, ONE, HALF, prae_korr, einb_rk, bul_stoe, implruku  *
.BA*)
***********************************************************************/
/*.BE*/

{
  #define RETURN(rc)      /* for shortness:                        */  \
    {                                                                  \
      vmfree(vmblock);    /* free already assigned memory          */  \
      return (rc);        /* and report error                      */  \
    }

  void        *vmblock;   /* List of dynamically allocated vectors and*/
                          /* matrices                                 */
  REAL        *yk,        /* [0..n-1] vector with value of solution   */
                          /* at right end point                       */
              *yaj,       /* [0..n-1] vector of modified y_start used */
                          /* for forming the Jacobi matrix            */
              *r,         /* [0..n-1] vector with left boundary       */
                          /* condition, which serves as the right hand*/
                          /* side for the linear system in the Newton */
                          /* method                                   */
              *rj,        /* [0..n-1] vector with left end boundary   */
                          /* condition of the modified initial value  */
                          /* problem                                  */
              *d,         /* [0..n-1] vector, the solution of the     */
                          /* linear system                            */
              **amat,     /* [0..n-1,0..n-1] array of the Jacobi      */
                          /* matrix for the  Newton step              */
              *yk2,       /* [0..n-1] vector with the values of the   */
                          /* solution of the DE at the left endpoint  */
                          /* for implruku()                           */
              *g,         /* [0..n-1] vector of weights for implruku()*/
              xk,         /* left end point for  awp()                */
              hk,         /* desired step size for awp()              */
              epsabs,     /* absolute accuracy desired in  awp()      */
              epsrel,     /* relative accuracy in awp()               */
              epsrel2,    /* relative error bound for implruku()      */
              epsrel3,    /* same as `epsrel2', but possibly altered  */
                          /* during last call of  implruku()          */
              rmax,       /* Maximum norm of the left boundary        */
                          /* condition                                */
              delta,      /* Aux variable for  Jacobi matrix          */
              mach1,      /* accuracy bounds depending on computer    */
              mach2;      /*                                          */
  long        aufrufe;    /* number of function evaluations in awp()  */
  int         *pivot,     /* Permutation vector for Gaussian elimin.  */
              i,          /* loop counter                             */
              jacobi,     /* counter for columns of Jacobi matrix     */
              mark,       /* error code from gauss()                  */
              sign_det,   /* sign of the permutation in gauss()       */
              fehler;     /* error code in  awp()                     */

#define SINGU   1         /* Name for the return value of gauss(),    */
                          /* which indicates that the matrix is       */
                          /* singular                                 */
#define ENGL45  6         /* awp() shall use the England embedding    */
                          /* formula of 4th and 5th order             */

#define FORMEL11      11  /* Macros for calling einb_rk(); to improve */
#define MIT_HULL       1  /* readability                              */
#define VON_VORNE   TRUE
#define NOT_SAVE   FALSE

#define MMAX          12  /* for calling implruku()                   */
#define STUETZ  "stuetz"
#define AUS           ""
#define PROT          ""


  mach1  = POW(MACH_EPS, (REAL)0.75);
  mach2  = (REAL)100.0 * MACH_EPS;


  if (epsrb < mach2)               /* check input  */
    return 1;
  if (b <= a)
    return 2;
  if (h <= mach2 * FABS(a))
    return 3;
  if (n <= 0)
    return 4;
  if (awpnumm < 1 || awpnumm > 5)
    return 5;


  /* -------- allocate storage for aux vectors ---------------------- */

  #define MYALLOC(n)  (REAL *)vmalloc(vmblock, VEKTOR, n, 0)
  vmblock = vminit();             /* initialize storage               */
  yk    = MYALLOC(n);             /* allocate buffer for local aux    */
  yaj   = MYALLOC(n);             /* vectors                          */
  r     = MYALLOC(n);
  rj    = MYALLOC(n);
  d     = MYALLOC(n);
  yk2   = MYALLOC(n);
  g     = MYALLOC(n);
  pivot = (int *)vmalloc(vmblock, VVEKTOR, n, sizeof(int));
  amat  = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  #undef MYALLOC
  if (! vmcomplete(vmblock))   /* insufficient memory?                */
    RETURN(9);

  *act_iter = 0;
  epsabs    = HALF * epsawp;
  epsrel    = epsabs;

  for (i = 0; i < n; i++)       /* weight vector for implruku(), set  */
    g[i] = ONE;                 /* up as all ones                     */
  /* `epsrel2' is assigned the same value as `epsrel', but at least   */
  /*  1e-10, because otherwise implruku() may not work correctly      */
  /*  depending on type of compiler and floating point type chosen    */
  /* (This affects mainly Borland C++ for OS/2 and `long double').    */
  epsrel2 = (epsrel < (REAL)1.0e-10) ? (REAL)1.0e-10 : epsrel;

  /*********************************************************************
  * If  y_start proves to be a sufficiently good approximation for y(a)*
  * in the following loop, we report the success and check the         *
  * remaining input data in awp().                                     *
  *********************************************************************/

  for ( ; ; )
  {
#ifdef DEBUG
    zeig(y_start, n);
#endif
    copy_vector(yk, y_start, n);
    xk = a;
    hk = h;

    switch (awpnumm)
    {
      case 1:  fehler = awp(&xk, b, n, dgl, yk, epsabs, epsrel, &hk,
                            ENGL45, fmax, &aufrufe); break;
      case 2:  fehler = prae_korr(&xk, yk, n, dgl, b, &hk, epsabs,
                                  epsrel, fmax, &aufrufe, ONE,
                                  VON_VORNE); break;
      case 3:  fehler = einb_rk(&xk, b, n, dgl, yk, epsabs, epsrel,
                                FORMEL11, MIT_HULL, VON_VORNE, NOT_SAVE,
                                fmax, &aufrufe); break;
      case 4:  fehler = bul_stoe(&xk, b, n, dgl, yk, epsabs, epsrel,
                                 &hk, ONE, VON_VORNE, fmax, &aufrufe);
                                 break;
      case 5:  epsrel3 = epsrel2;
               copy_vector(yk2, yk, n);                   /* yk2 = yk */
               fehler = implruku(dgl, n, MMAX, STUETZ, AUS, PROT,
                                 - ONE, &epsrel3, fmax, &aufrufe, g,
                                 &xk, b, yk2, yk); break;
      default: RETURN(5);
    }
    if (fehler != 0)                         /* error in IVP solver ? */
      RETURN(10 + fehler);

#ifdef DEBUG
    zeig(yk, n);
#endif
    (*rand)(y_start, yk, r);

    for (i = 0, rmax = ZERO; i < n; i++)
      rmax = max(rmax, FABS(r[i]));
    if (rmax < epsrb)               /* boundary condition satisfied ? */
      RETURN(0);                           /* report success          */

    if (++*act_iter > itmax)    /* approximation not good enough      */
    {                           /* after  itmax Newton steps ?        */
      RETURN(6);                /* report failure                     */
    }

    /* find a better approximation  y_start for y(a) by using the     */
    /* Newton method with a Jacobi matrix that is approximated by one-*/
    /* sides difference quotients                                     */

    for (jacobi = 0; jacobi < n; jacobi++)
    {
      copy_vector(yk,  y_start, n);
      copy_vector(yaj, y_start, n);

      if (FABS(yk[jacobi]) < mach2)
      {
        yk[jacobi] += mach1;
        delta      =  ONE / mach1;
      }
      else
      {
        yk[jacobi] *= ONE + mach1;
        delta      =  ONE / (mach1 * yk[jacobi]);
      }
      yaj[jacobi] = yk[jacobi];
      xk = a;
      hk = h;
#ifdef DEBUG
      zeig(yaj, n);
#endif
      switch (awpnumm)
      {
        case 1:  fehler = awp(&xk, b, n, dgl, yk, epsabs, epsrel, &hk,
                              ENGL45, fmax, &aufrufe); break;
        case 2:  fehler = prae_korr(&xk, yk, n, dgl, b, &hk, epsabs,
                                    epsrel, fmax, &aufrufe, ONE,
                                    VON_VORNE); break;
        case 3:  fehler = einb_rk(&xk, b, n, dgl, yk, epsabs, epsrel,
                                  FORMEL11, MIT_HULL, VON_VORNE,
                                  NOT_SAVE, fmax, &aufrufe); break;
        case 4:  fehler = bul_stoe(&xk, b, n, dgl, yk, epsabs, epsrel,
                                   &hk, ONE, VON_VORNE, fmax, &aufrufe);
                 break;
        case 5:  epsrel3 = epsrel2;
                 copy_vector(yk2, yk, n);                 /* yk2 = yk */
                 fehler = implruku(dgl, n, MMAX, STUETZ, AUS, PROT,
                                   -ONE, &epsrel3, fmax, &aufrufe, g,
                                   &xk, b, yk2, yk); break;
        default: RETURN(5);
      }
      if (fehler != 0)                       /* error in IVP solver ? */
        RETURN(20 + fehler);

#ifdef DEBUG
      zeig(yk, n);
#endif
      (*rand)(yaj, yk, rj);

      for (i = 0; i < n; i++)
        amat[i][jacobi] = (rj[i] - r[i]) * delta;
    }

    mark = gauss(0, n, amat, amat, pivot, r, d, &sign_det);
    if (mark == SINGU)                     /* Jacobi matrix singular? */
      RETURN(8);                           /* return error            */

    for (i = 0; i < n; i++)                        /* correct y_start */
      y_start[i] -= d[i];
  }
}

/* ---------------------------- END rwp.c --------------------------- */
