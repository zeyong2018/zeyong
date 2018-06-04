#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE ab_mou.c ------------------------- */

/***********************************************************************
*                                                                      *
* Solve an ordinary system of differential equations of first order    *
* -------------------------------------------------------------------- *
* using the predictor-corrector method of Adams-Bashforth-Moulton      *
* ----------------------------------------------------------------     *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Jobst Hoffmann (FORTRAN)                       *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 7.9.1992                                       *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  MACH_EPS, FABS, max, dglsysfnk, min,  */
                       /*       boolean, FALSE, TRUE, norm_max,       */
                       /*       copy_vector, REAL, ZERO, ONE, THREE,  */
                       /*       FIVE, SIX, EIGHT, TEN, NINE           */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR                                */
#include <ab_mou.h>    /*  for  prae_korr                             */




/* ------------------------------------------------------------------ */

typedef struct     /* Type definitions for aux vectors ...............*/
/*.IX{hilfstyp}*/
{
  REAL *f[5],      /* buffder for the starting values:  f[1],...,f[4] */
                   /* contain the newest entries from rk_start() and  */
                   /* abm_schritt(), f[0] is also needed              */
       *tmp,       /* aux buffer for the vectors k[i] from one Runge- */
                   /* Kutta step; used in abm_schritt() for a linear  */
                   /* combination of nodes in correktor step          */
       *ki_sum,    /* Linear combination A[0]*k[0]+..+A[3]*k[3] in the*/
                   /* Runge-Kutta method and end result of            */
                   /* rk_schritt()                                    */
       *y1,        /* in rk_start() solution for step size  3*h, or   */
                   /* in abm_schritt() solution of the predictor step */
                   /* (also used fro aux purposes in rk_schritt() )   */
       *y2,        /* newest approximation                            */
       *diff;      /* eror estimate for computed solution             */
} hilfstyp;



/* ------------------------------------------------------------------ */

static void init0_vector
/*.IX{init0\unt vector}*/
                        (
                         REAL vektor[],
                         int  n
                        )

/***********************************************************************
* initialize the [0..n-1] vector vektor to be the zero vector          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO                                                           *
***********************************************************************/

{
  for (n--; n >= 0; n--)
    *vektor++ = ZERO;
}



/* ------------------------------------------------------------------ */

static void inc_vector
/*.IX{inc\unt vector}*/
                      (
                       REAL ziel[],
                       REAL quelle[],
                       REAL faktor,
                       int  n
                      )

/***********************************************************************
* add factor * quelle to the vector ziel                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  for (n--; n >= 0; n--)
    *ziel++ += faktor * *quelle++;
}



/* ------------------------------------------------------------------ */

static void add_vector
/*.IX{add\unt vector}*/
                      (
                       REAL summe[],
                       REAL summand1[],
                       REAL summand2[],
                       REAL faktor,
                       int  n
                      )

/***********************************************************************
* add factor * summand2 to the [0..n-1] vector summand1, store result  *
* summe.                                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  for (n--; n >= 0; n--)
    *summe++ = *summand1++ + faktor * *summand2++;
}



/* ------------------------------------------------------------------ */

static void rk_schritt
/*.IX{rk\unt schritt}*/
                      (
                       REAL      x0,
                       REAL      y0[],
                       int       n,
                       dglsysfnk dgl,
                       hilfstyp  *hilf,
                       REAL      h
                      )

/***********************************************************************
* perform one Runge-Kutta integration                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0    x-value from wchich to integrate                               *
* y0    initial value of solution at  x0                               *
* n     number of DEs                                                  *
* dgl   pointer to function that evaluates the right hand side of the  *
*       the DE system  y' = f(x,y)                                     *
* h     step size for current Runge-Kutta integration                  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* hilf  stores the result in ki_sum. The entries in y1 and tmp have no *
*       meaning.                                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hilfstyp, inc_vector, add_vector, REAL, dglsysfnk, ZERO, ONE, SIX,   *
* THREE                                                                *
***********************************************************************/

{
  /* --- Coefficients of the classical Runge-Kutta method ----------- */
  static REAL A[4] = { ONE / SIX, ONE / THREE, ONE / THREE, ONE / SIX },
              a[4] = { ZERO, HALF, HALF, ONE };  /* use also as the   */
                                                 /* b[j][s]           */

  int         i, j;                              /* loop counter      */


  (*dgl)(x0, y0, hilf->tmp);    /* k[0] <- right hand side at (x0,y0) */

  for (i = 0; i < n; i++)                    /* initialize ki_sum with*/
    hilf->ki_sum[i] = A[0] * hilf->tmp[i];   /* first term            */

  for (j = 1; j < 4; j++)      /* compute k[1]..k[3], add to ki_sum   */
  {                            /* multiplied by factor                */
    add_vector(hilf->y1, y0,                   /* y1 <- y0 +          */
               hilf->tmp, a[j] * h, n);        /* a[j] * h * k[j - 1] */
    (*dgl)(x0 + a[j] * h, hilf->y1, hilf->tmp);     /* compute k[j]   */

    inc_vector(hilf->ki_sum,
               hilf->tmp, A[j], n);          /* ki_sum += A[j] * k[j] */
  }
}



/* ------------------------------------------------------------------ */

static int rk_start
/*.IX{rk\unt start}*/
                   (
                    REAL      x,
                    REAL      *x0,
                    REAL      y[],
                    int       n,
                    dglsysfnk dgl,
                    REAL      xend,
                    REAL      *h,
                    REAL      hmax,
                    int       new_step,
                    int       *methode,
                    long      *aufrufe,
                    hilfstyp  *hilf
                   )

/***********************************************************************
* Find the starting values using the Runge-Kutta method as needed in   *
* prae_korr() for the Adams-Bashforth-Moulton method                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x         x-value for start of integration                           *
* y         Initial value for the solution at x                        *
* n         number of DE equations                                     *
* dgl       pointer to a function that evaluates the right hand side   *
*           of the system of DEs   y' = f(x,y)                         *
* xend      x-value where solution is wanted; xend may be less than x  *
* h         Step size                                                  *
* hmax      maximal step size, must be positive                        *
* new_step  != 0: check new step size for properness                   *
*            = 0: do not check new step size                           *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x0       x-value until which we have integrated                      *
* h        step size for next integration                              *
* methode  value 0; hence the error estimate in prae_korr() is per-    *
*          formed for the factor for Runge-Kutta values                *
* aufrufe  actual number of calls of dgl()                             *
* hilf     Structure, which contains info of the solution:             *
*          y1:         approximate solution needed for error estimation*
*                      (Step size  3*h)                                *
*          y2:         actaul approximate solution                     *
*          f[2]..f[4]: starting values                                 *
*          The entries in  ki_sum and tmp represent aux values.        *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok                                                          *
* = 1: new step size too small relative to machine constant            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hilfstyp, rk_schritt, inc_vector, add_vector, REAL, MACH_EPS, min,   *
* copy_vector, FABS, THREE, EIGHT, dglsysfnk                           *
***********************************************************************/

{
  int j;                                          /* loop counter     */


  if (new_step)                 /* check new step size for properness */
  {
    *h = min(*h, hmax);
    *h = min(FABS(*h), FABS(xend - x) / THREE);
    *h = (xend > x) ? *h : -*h;
    if (FABS(*h) <= EIGHT * MACH_EPS * FABS(x))
      return 1;
  }

  *x0 = x;                                      /* save initial value */


  copy_vector(hilf->y2, y, n);                             /* y2 <- y */

  for (j = 2; j < 5; j++)                 /* three steps with steps h */
  {
    rk_schritt(x, hilf->y2, n, dgl, hilf, *h);   /* accumulate ki_sum */
    x += *h;
    inc_vector(hilf->y2, hilf->ki_sum, *h, n);    /* y2 += h*ki_sum   */
    (*dgl)(x, hilf->y2, hilf->f[j]);    /* copy remainder of starting */
  }                                     /* values in  f[2]..f[4]      */


  /* after three steps with size h, we now perform one step with 3*h  */
  /* result put into y1.                                              */

  rk_schritt(*x0, y, n, dgl, hilf, THREE * *h);   /* compute ki_sum   */
  add_vector(hilf->y1, y, hilf->ki_sum, THREE * *h, n);


  *x0      += THREE * *h;
  *methode =  0;
  *aufrufe += 13;       /* 13 calls of dgl() for the starting values  */

  return 0;
}



/* ------------------------------------------------------------------ */

static void abm_schritt
/*.IX{abm\unt schritt}*/
                       (
                        REAL      *x0,
                        int       n,
                        dglsysfnk dgl,
                        REAL      *h,
                        int       *methode,
                        long      *aufrufe,
                        hilfstyp  *hilf
                       )

/***********************************************************************
* Perform one step of the Adams-Bashforth-Moulton method               *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0    initial x-value                                                *
* n     number of DEs in system                                        *
* dgl   pointer to the function with the right hand side of the DE     *
*       system   y' = f(x,y)                                           *
* h     step size                                                      *
* hilf  Structure with the vectors f[1]..f[4] for the starting values  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* x0       final x-value of the integration                            *
* methode  equal to 1; hencethe error estimation in prae_korr() is     *
*          performed with the factor for Adams-Bashforth-Moulton.      *
* aufrufe  current number of calls of dgl()                            *
* hilf     Structure with the following output :                       *
*          y1:   approximate solution used for the error estimation    *
*                (from predictor step)                                 *
*          y2:   approximate solution from corrector step              *
*          f[4]: new node for starting values                          *
*          f[0],...,f[3] are alterd and tmp is used as aux storage     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hilfstyp, init0_vector, inc_vector, add_vector, REAL, dglsysfnk,     *
* ONE, FIVE, NINE                                                      *
***********************************************************************/

{
  /* --- the coefficients for  Adams-Bashforth-Moulton -------------- */
  static REAL prae[4] = { -NINE, (REAL)37.0,
                          (REAL)-59.0, (REAL)55.0 };    /* Predictor  */
  static REAL korr[4] = { ONE, -FIVE,
                          (REAL)19.0, NINE };           /* corrector  */

  REAL        *tmp;     /* aux vector for cyclic exchanges of adresses*/
                        /* of starting vectors                        */
  int         j;        /* loop counter                               */


  tmp = hilf->f[0];                /* the starting values are expected*/
  for (j = 0; j < 4; j++)          /* reside in f[0], ..., f[3] , but */
    hilf->f[j] = hilf->f[j + 1];   /* are stored in f[1]..f[4], we    */
  hilf->f[4] = tmp;                /* rotate their pointer            */

  init0_vector(hilf->y1, n);               /* one predictor step      */
  for (j = 0; j < 4; j++)
    inc_vector(hilf->y1, hilf->f[j], prae[j], n);
  add_vector(hilf->y1, hilf->y2, hilf->y1, *h / (REAL)24.0, n);

  *x0 += *h;                            /* move on by h; compute rhs  */
  (*dgl)(*x0, hilf->y1, hilf->f[4]);    /* of DE at (x0,y1) : this    */
                                        /* yields the start for the   */
                                        /* corrector                  */

  init0_vector(hilf->tmp, n);               /* one corrector step     */
  for (j = 0; j < 4; j++)
    inc_vector(hilf->tmp, hilf->f[j + 1], korr[j], n);
  inc_vector(hilf->y2, hilf->tmp, *h / (REAL)24.0, n); /* y2: new appr*/
                                                       /* solution    */
  (*dgl)(*x0, hilf->y2, hilf->f[4]);     /* f[4] comtains the node for*/
                                         /* the starting values       */
                                        /* derived from the corrector */

  *methode =  1;
  *aufrufe += 2;              /* 2 calls of dgl() for one  A-B-M step */
}



/* ------------------------------------------------------------------ */

static void *init_praeko
/*.IX{init\unt praeko}*/
                        (
                         hilfstyp *hilf,
                         int      n
                        )

/***********************************************************************
* allocate storage for the dynamically allocated aux vectors, which    *
* reside inside the structure  hilf                                    *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* hilf  Structure, giving valid addresses for all pointers             *
*                                                                      *
* Return value :                                                       *
* =============                                                        *
* leading address for storage, in which all dynamically allocated aux  *
* vectors are managed. Zero address if lack of memory is experienced   *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hilfstyp, vminit, vmalloc, vmcomplete, vmfree, VEKTOR, NULL, REAL    *
***********************************************************************/

{
  void *vmblock;                  /* List of dynamically allocated    */
                                  /* vectors and matrices             */
  int  i;                         /* loop counter                     */


  #define MYALLOC(n)  (REAL *)vmalloc(vmblock, VEKTOR, n, 0)

  vmblock = vminit();                 /* initialize storage           */
  for (i = 0; i < 5; i++)
    hilf->f[i] = MYALLOC(n);
  hilf->tmp    = MYALLOC(n);
  hilf->ki_sum = MYALLOC(n);
  hilf->y1     = MYALLOC(n);
  hilf->y2     = MYALLOC(n);
  hilf->diff   = MYALLOC(n);

  #undef MYALLOC

  if (! vmcomplete(vmblock))                     /* lack of storage ? */
  {
    vmfree(vmblock);                    /* free buffers, report error */
    return NULL;
  }

  return vmblock;                             /* return valid address */
}



/* ------------------------------------------------------------------ */
/*.BA*/

int prae_korr   /* Predictor-corrector meth. for 1st order DE systems */
/*.IX{prae\unt korr}*/
             (
              REAL      *x,          /* initial/final x-value ........*/
              REAL      y[],         /* initial value/ solution ......*/
              int       n,           /* number of DEs ................*/
              dglsysfnk dgl,         /* right hand side for DEs ......*/
              REAL      xend,        /* desired final x-value ........*/
              REAL      *h,          /* starting/final step size .....*/
              REAL      epsabs,      /* absolute error bound .........*/
              REAL      epsrel,      /* relative error bound .........*/
              long      fmax,        /* maximal # of calls for dgl() .*/
              long      *aufrufe,    /* actual # of calls of dgl() ...*/
              REAL      hmax,        /* maximal step size ............*/
              boolean   neu          /* delete old data ? ............*/
             )                       /* error code ...................*/

/***********************************************************************
* Solve a system of ordinary differential equations                    *
*                                                                      *
*      y' = f(x,y)        with the initial condition y(x0) = y0        *
*                                                                      *
* using the predictor-corrector method of Adams-Bashforth-Moulton.     *
* The needed starting values are obtained by an initial run of Runge-  *
* Kutta with the same order as the  A-B-M method.                      *
* We work with automatic step size control with doubling/halfing of    *
* the step size according to the subsequent error estimates.           *
.BE*)
*                                                                      *
* REMARK :                                                             *
*   This procedure frees certain buffers only at the beginning of a    *
*   new run or when memory is not sufficient, or for improper input.   *
*   Hence we advise to terminate a connected run of this function with *
*   one call with improper input such as n < 0 to clear all memory.    *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        x-value where a corresponding y-value is known              *
* y        [0..n-1] vector with the initial y-data                     *
* dgl      pointer to a function f which evaluates the right hand side *
*          of the DE system   y' = f(x,y)                              *
* n        number of DEs                                               *
* xend     x-value for which we want to find the solution; may be less *
*          than x                                                      *
* h        step size for next step; this is usually determined inside  *
*          this function                                               *
* epsabs\  error bounds for absolute and relative errors; each >= 0    *
* epsrel/  We apply a mixed test :                                     *
*              |local error|  <=  |y| * epsrel + epsabs.               *
*          For  epsrel = 0, we test for absolute error only.           *
*          For  epsabs = 0, we test for relative error.                *
*          Each error bound is internally set to 10 * machine constant *
*          if it should be prescribed to be too small.                 *
* fmax     upper bound for the allowed number of evaluations of the    *
*          right hand side  dgl() od the DE system                     *
* hmax     maximal  step size; must be positive                        *
* neu      If TRUE, we start off using Runge-Kutta, even if the earlier*
*          integration counld be continued with another A-B-M step     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value of integration (if run was successful: x=xend)*
* y        [0..n-1] approximate solution vector at x                   *
* h        final step size; should be used for start of next call.     *
*          If changed, please reset flag as well.                      *
* aufrufe  counter for calls of  dgl()                                 *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* error code :                                                         *
*   = 0: no error, xend was reached                                    *
*   = 1: after fmax calls of dgl() we have not reached xend:           *
*        a new call with unchanged parameters might be successful      *
*        (or try raising the error bounds)                             *
*   = 2: the step size is less than 8 * machine constant.              *
*        Before subsequent calls increase h and the error bounds.      *
*   = 3: epsabs or epsrel negative, or both equal to zero.             *
*   = 4: xend = x                                                      *
*   = 5: hmax negative                                                 *
*   = 6: n <= 0                                                        *
*   = 7: lack of available memory                                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* rk_start, abm_schritt, init_praeko, MACH_EPS, max, boolean, EIGHT,   *
* norm_max, copy_vector, FABS, vmfree, NULL, REAL, ZERO, ONE, TEN,     *
* dglsysfnk, FALSE, TRUE                                               *
.BA*)
***********************************************************************/
/*.BE*/

{
#define CHECK  TRUE           /* check new step size in rk_start()    */

  /* Vector with factors for error estimation :                       */
  /* guess[0] for Runge-Kutta, guess[1] for Adams-Bashforth-Moulton   */
  static
    REAL     guess[2] = { ONE / (REAL)80.0,
                          (REAL)-19.0 / (REAL)270.0 },

             x0,              /* aux storage for       x              */
             h_save;          /* aux storage for       h              */
  static
    int      methode;         /* Number of most recently used method  */
                              /* (Runge-Kutta or Adams-Bashforth-M.)  */
  static
    hilfstyp hilf;            /* local aux variables                  */

  static
    void     *vmblock = NULL; /* List of dynamic allocations          */
  REAL       ynorm,           /* Maximum norm from hilf.y2            */
             diffnorm;        /* Maximum norm from hilf.diff          */
  int        i;               /* loop counter                         */
  boolean    folgeaufruf = FALSE;  /* Flag indicating we can continue */
                              /* an earlier call with A-B-M method    */


  /* ------------------- Check input data --------------------------- */

  if (epsabs < ZERO || epsrel < ZERO || epsabs + epsrel <= ZERO)
  {
    vmfree(vmblock);
    vmblock     = NULL;
    folgeaufruf = FALSE;
    return 3;
  }
  if (xend == *x)
  {
    vmfree(vmblock);
    vmblock     = NULL;
    folgeaufruf = FALSE;
    return 4;
  }
  if (hmax <= ZERO)
  {
    vmfree(vmblock);
    vmblock     = NULL;
    folgeaufruf = FALSE;
    return 5;
  }
  if (n <= 0)
  {
    vmfree(vmblock);
    vmblock     = NULL;
    folgeaufruf = FALSE;
    return 6;
  }


  /* ------------- Prepare integration loop ------------------------- */

  epsrel = max(epsrel, TEN * MACH_EPS);
  epsabs = max(epsabs, TEN * MACH_EPS);

  *aufrufe = 0l;

  if (neu)                           /*  Force start with R-K step ?  */
    folgeaufruf = FALSE;


  if (folgeaufruf)                   /* Use values of previous call ? */
    abm_schritt(&x0, n, dgl, h,      /* Perform A-B-M step            */
                &methode, aufrufe, &hilf);

  else                               /* very first call ?             */
  {
    vmfree(vmblock);                  /* clear storage from previous  */
                                      /* calls                        */
    if ((vmblock = init_praeko(&hilf, n)) == NULL) /* allocate aux    */
      return 7;                                    /* vectors         */

    (*dgl)(*x, y, hilf.f[1]);      /* store beginning of starting     */
    ++*aufrufe;                    /* values in hilf.f[1]             */
    h_save = *h;                   /* save starting step size         */
    if (rk_start(*x, &x0, y, n, dgl,     /* new step size too small ? */
                 xend, h, hmax, CHECK,
                 &methode, aufrufe, &hilf))
    {
      *h          = h_save;
      folgeaufruf = FALSE;
      return 0;
    }
    /* --------- starting values available in                -------- */
    /* --------- hilf.f[1], hilf.f[2], hilf.f[3], hilf.f[4]. -------- */
  }


  /* --------------------- Integration  loop ------------------------ */

  for ( ; ; )
  {
    if (*aufrufe > fmax)                /* excessive function calls ? */
    {
      *x = x0;
      copy_vector(y, hilf.y2, n);   /* newest approximation for  y    */
      folgeaufruf = TRUE;
      return 1;                     /* stop and report excessive calls*/
    }

    for (i = 0; i < n; i++)                   /* errro estimation     */
      hilf.diff[i] = guess[methode] *
                     (hilf.y2[i] - hilf.y1[i]);
    diffnorm = norm_max(hilf.diff, n);
    ynorm    = norm_max(hilf.y2,   n);


    if (diffnorm >= epsrel * ynorm + epsabs)   /* error too large ?   */
    {
      *h *= HALF;                              /* halve the step size */
                                               /* and  repeat         */
      if (FABS(*h) <= EIGHT *                  /* step size too small?*/
                      MACH_EPS * FABS(x0))
      {
        vmfree(vmblock);
        vmblock     = NULL;
        folgeaufruf = FALSE;
        return 2;
      }

      rk_start(*x, &x0, y, n, dgl,          /* compute new starting   */
               xend, h, hmax, ! CHECK,      /* values with  Runge-    */
               &methode, aufrufe, &hilf);   /* Kutta                  */
    }

    else                             /* error not excessive ?         */
    {                                /* step was successful, continue */
      *x = x0;                       /* on with the previous step size*/
      for (i = 0; i < n; i++)        /* add estimated error onto new  */
        hilf.y2[i] += hilf.diff[i];  /* approximation                 */
      copy_vector(y, hilf.y2, n);    /* newest approximation for y    */
      (*dgl)(x0, y, hilf.f[4]);      /* correct last node for the next*/
                                     /* A-B-M step                    */
      ++*aufrufe;

      if (diffnorm <= (REAL)0.02 *             /* accuracy exessive ? */
                      (epsrel * ynorm + epsabs))
      {
        *h     += *h;                             /* double step size */
        h_save =  max(*h, h_save);

        if ((*h > ZERO && x0 >= xend) ||           /* reached xend ? */
            (*h < ZERO && x0 <= xend))
        {
          folgeaufruf = FALSE;
          return 0;                                 /* all done !     */
        }

        /* Continue integration with doubled step size.               */
        /* First find a new set of starting values via Runge-Kutta.   */

        copy_vector(hilf.f[1],  /* use final value of starting value  */
                    hilf.f[4],  /* as the first entry for the new set */
                    n);
        if (rk_start(*x, &x0, y, n, dgl, /* new step size too small ? */
                     xend, h, hmax, CHECK,
                     &methode, aufrufe, &hilf))
        {
          *h          = h_save;
          folgeaufruf = FALSE;
          return 0;
        }
      }

      else                                  /* last step successful ? */
      {
        if ((*h > ZERO && x0 >= xend) ||    /* has xend been reached ?*/
            (*h < ZERO && x0 <= xend))
        {
          folgeaufruf = FALSE;
          return 0;
        }

        if ((*h > ZERO && x0 + *h >= xend) ||    /* sufficintly close */
            (*h < ZERO && x0 + *h <= xend))      /* to  xend ?        */
        {
          /* The next step would lead beyond xend. Hence we noe use   */
          /* xend - x0 as the step size, and start anew using R-K.    */

          h_save = *h = xend - x0;
          copy_vector(hilf.f[1],   /* assign the final entry of the   */
                      hilf.f[4],   /* old starting values to be the   */
                      n);          /* first value of the new one      */
          if (rk_start(*x, &x0, y, n, dgl, /* new step size too small?*/
                       xend, h, hmax, CHECK,
                       &methode, aufrufe, &hilf))
          {
            *h          = h_save;
            folgeaufruf = FALSE;
            return 0;
          }
        }

        else
          abm_schritt(&x0, n, dgl, h,       /* perform one A-B-M step */
                      &methode,
                      aufrufe, &hilf);
      }
    }
  }
}

/* -------------------------- END ab_mou.c -------------------------- */
