#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE bulirsch.c ------------------------ */

/***********************************************************************
*                                                                      *
* Solve a system of first degree ordinary differential equations using *
* -------------------------------------------------------------------- *
* the extrapolation method of Bulirsch-Stoer-Gragg                     *
* -------------------------------------------------                    *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Jobst Hoffmann (FORTRAN)                       *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 3.11.1992                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  MACH_EPS, FABS, max, dglsysfnk, min,  */
                       /*       boolean, FALSE, TRUE, LOG, POW, FOUR, */
                       /*       copy_vector, REAL, ZERO, HALF, ONE,   */
                       /*       TEN                                   */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR                                */
#include <bulirsch.h>  /*  for  bul_stoe                              */



/* ------------------------------------------------------------------ */

static int ext_max;             /* maximal order of the extrapolation */
static int bufol[12] = {  2,  4,  6,  8, 12,  16,    /* the Bulirsch  */
                         24, 32, 48, 64, 96, 128 };  /* sequence      */



/* ------------------------------------------------------------------ */

static void extrapol
/*.IX{extrapol}*/
    (
     int       *row,
     REAL      *fhilf[],
     REAL      y[],
     int       n,
     REAL      epsabs,
     REAL      epsrel,
     REAL      *x,
     REAL      *h,
     REAL      h0,
     boolean   *ahead,
     int       *index
    )

/***********************************************************************
* Perform one extrapolation step for bul_stoe()                        *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* row     pointer to the arrays  bufol and fhilf                       *
* fhilf   Matrix for extrapolation values                              *
* y       y-value of solution at x                                     *
* n       Anzahl der Differentialgleichungen                           *
* epsabs  absolute error bound                                         *
* epsrel  relative error bound                                         *
* x       x-value                                                      *
* h       local step sizes                                             *
* h0                                                                   *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ahead  *ahead  = 1:  Step is acceptable                              *
*        *ahead != 1:  Step not acceptable                             *
* index  pointer to the arrays bufol and fhilf                         *
* x      final x-value                                                 *
* h      final step size                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ext_max, bufol, REAL, TRUE, min, max, boolean, copy_vector, FABS,    *
* POW, ZERO, ONE                                                       *
***********************************************************************/

{
  int  column,      /* column of extrapolation scheme                 */
       i;           /* loop counter                                   */
  REAL diff_max,    /* maximal difference of two columns in scheme    */
       fhilf_i1,    /* column in extrapolation scheme                 */
       fhilf_i2,    /* column to the right of fhilf_i1                */
       bufol_r,     /* element in Bulirsch sequence                   */
       bufol_i,     /* ditto                                          */
       y_max,       /* maximal value in a column                      */
       help;        /* aux variable                                   */


  for (column = 2, y_max = ZERO;
       column <= min(*row + 1, ext_max) && (! *ahead);
       column++)
  {
    *index = min(11 - column, *row - column + 3);
    for (i = 0, diff_max = ZERO; i < n; i++)
    {
      fhilf_i1 = fhilf[*index - 1][i];
      fhilf_i2 = fhilf[*index - 2][i];
      bufol_r  = (REAL)bufol[*row];
      bufol_i  = (REAL)bufol[*index - 2];
      fhilf[*index - 2][i] = fhilf_i1 + (fhilf_i1 - fhilf_i2) /
                             ((bufol_r / bufol_i) * (bufol_r / bufol_i)
                              - ONE);
      fhilf_i1 = fhilf[*index - 1][i];
      fhilf_i2 = fhilf[*index - 2][i];
      y_max    = max(y_max, FABS(fhilf_i2));
      diff_max = max(diff_max, FABS(fhilf_i2 - fhilf_i1));
    }
    if (diff_max < epsrel * y_max + epsabs)    /* Step acceptable ?   */
    {
      *x += h0;
      copy_vector(y, fhilf[*index - 2], n);
      help   = (REAL)(column - ext_max);
      *h     = (REAL)0.9 * *h *            /* Step size for next step */
               POW((REAL)0.6, help);
      *row   = -1;
      *ahead = TRUE;
    }
  }
}



/* ------------------------------------------------------------------ */
/*.BA*/

int bul_stoe       /* Extrapolation method for 1st order DE systems ..*/
/*.IX{bul\unt stoe}*/
    (
     REAL      *x,             /* initial x-value/ final x-value .....*/
     REAL      xend,           /* desired end point ..................*/
     int       n,              /* number of DEs ......................*/
     dglsysfnk dgl,            /* right hand side of DE system .......*/
     REAL      y[],            /* initial y-value/ final y-value .....*/
     REAL      epsabs,         /* absolute error bound ...............*/
     REAL      epsrel,         /* relative error bound ...............*/
     REAL      *h,             /* initial/final step size ............*/
     REAL      hmax,           /* maximal step size ..................*/
     int       neu,            /* use an outside given x ? ...........*/
     long      fmax,           /* maximal # of calls of  dgl() .......*/
     long      *aufrufe        /* actual # of calls of dgl() .........*/
    )                          /* error code .........................*/

/***********************************************************************
* Compute the solution of an intial value problem for a first order    *
* system of ordinary differential equations                            *
*                                                                      *
*      y' = f(x,y)        with initial value y(x0) = y0                *
*                                                                      *
* by using the extrapolation method of Bulirsch-Stoer-Gragg.           *
* The maximal extrapolation order is set internally depending on the   *
* machine constant; the step size control follows the ideas in :       *
*     HALL, G.; WATT, J.M: Modern Numerical Methods for Ordinary       *
*                          Differential Equations, Oxford 1976,        *
*                          Clarendon Press,  [HALL76]                  *
.BE*)
* REMARK :                                                             *
* The dynamic allocations from this function are only set free for a   *
* new run or with improper input (fehler = TRUE or Return value >= 2). *
* Hence we advise the user to end work with this function by one call  *
* with an improper value such as n = -2 in order to free up all storage*
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        x-value                                                     *
* y        [0..n-1] vector with initial y-value at x                   *
* dgl      pointer to a function which evaluates the right hand side   *
*          odf the DE system   y' = f(x,y)                             *
* n        number of equations in the system                           *
* xend     desired final value for x, xend < x is allowed              *
* h        step size for next step                                     *
* hmax     maximal step size; hmax > 0                                 *
* epsabs\  error bounds for absolute and relative errors, both non     *
* epsrel/  negative.                                                   *
*          We test for the mixed error:                                *
*              |lokaler Fehler|  <=  |y| * epsrel + epsabs.            *
*          If epsrel = 0, we test the absolute error; if epsabs = 0 we *
*          test the relative error. epsrel will be set to 10 * machine *
*          constant if it has been chosen smaller than this.           *
* neu      Flag, indicating whether this function is called to compute *
*          on after one partial run that was cut short due to too many *
*          evaluations:                                                *
*          neu = FALSE: continue with old values                       *
*          neu = TRUE:  start new                                      *
* fmax     upper bound for the number of allowed function evaluations  *
*          of the right hand side via dgl()                            *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value of integration; usually x = xend              *
* y        [0..n-1] y-value vector at x                                *
* h        final step size; should be used unchanged in the next call  *
*          if neu = TRUE; should be redefined for neu = TRUE.          *
* aufrufe  number of calls of dgl()                                    *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok; after a new setting for xend, the function can be used  *
*      again for the same problem, or all inputs may be changed        *
* = 1: method has not reached xend after  MAX_FCT right hand side      *
*      evaluations; try another call with same parameters or change    *
*      error bounds.                                                   *
* = 2: step size less than 4 * machine constant.                       *
*      For a future call, increase step size and error bounds.         *
* = 3: epsabs or epsrel negative, or both zero.                        *
* = 4: xend = x                                                        *
* = 5: hmax negative.                                                  *
* = 6: n <= 0.                                                         *
* = 7: lack of memory                                                  *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* extrapol, ext_max, bufol, REAL, TRUE, FALSE, MACH_EPS, min,          *
* max, boolean, copy_vector, FABS, POW, LOG, vminit, vmalloc,          *
* vmcomplete, vmfree, VEKTOR, dglsysfnk, ZERO, HALF, TEN               *
.BA*)
***********************************************************************/
/*.BE*/

{
  static void *vmblock = NULL;    /* List of dynamic allocations      */
  static REAL *fhilf[12],     /* Index 0..7: Matrix for the           */
                              /*             extrapolation scheme     */
                              /* Index 8..11: aux vectors             */
              x0;             /* x-value on call of dgl()             */
  static int  row;            /* row of extrapolation scheme          */
  REAL        absh,           /* magnitude of step size               */
              absh0,          /* magnitude of aux step size  h0       */
              h0 = ZERO,      /* aux step size                        */
              hilf;           /* aux variable                         */
  int         i, j,           /* loop counters                        */
              count,          /* Loop variable for mid-point rule     */
              index;          /* Index for the vectors fhilf and bufol*/
  boolean     ahead;          /* Flag indicating acceptance of step   */


  if (epsabs < ZERO ||                                 /* check input */
      epsrel < ZERO ||
      epsabs + epsrel <= ZERO)
  {
    vmfree(vmblock);
    vmblock = NULL;
    return 3;
  }
  else if (xend == *x)
  {
    vmfree(vmblock);
    vmblock = NULL;
    return 4;
  }
  else if (hmax <= ZERO)
  {
    vmfree(vmblock);
    vmblock = NULL;
    return 5;
  }
  if (n <= 0)
  {
    vmfree(vmblock);
    vmblock = NULL;
    return 6;
  }


  if (! neu) /* new call or repeat call due to excessive evaluations? */
  {
    *x = fhilf[9][0];
    copy_vector(y, fhilf[8], n);
  }
  else
  {                                               /* allocate storage */
    vmfree(vmblock);
    #define MYALLOC(n)  (REAL *)vmalloc(vmblock, VEKTOR, n, 0)
    vmblock = vminit();
    for (i = 0; i < 12; i++)
      fhilf[i] = MYALLOC(n);
    #undef MYALLOC
    if (! vmcomplete(vmblock))                    /* lack of memory ? */
    {
      vmfree(vmblock);               /* free storage and report error */
      return 7;
    }
    row = -1;
  }

  *aufrufe = 0l;
  ahead    = TRUE;
  ext_max  = (int)(-LOG(MACH_EPS) / LOG(TWO) / (REAL)7.0 + HALF);

  epsrel   = max(epsrel, TEN * MACH_EPS);


  for ( ; ; )
  {
    if (neu)
    {
      if (ahead)                                       /* new step ? */
      {
        absh = FABS(*h);
        absh = min(absh, hmax);
        h0 = min(absh, FABS(xend - *x));
        h0 = (xend > *x) ? h0 : -h0;
        absh0 = FABS(h0);
        if (absh0 <= FOUR * MACH_EPS * FABS(*x))
          return 0;
        ahead = FALSE;
      }
      do
      {
        row++;                    /* find step size for extrapolation */
        *h = h0 / bufol[row];
        x0 = *x;                  /* Euler step; save initial values  */
        copy_vector(fhilf[8], y, n);
        (*dgl)(x0, fhilf[8], fhilf[11]);
        for (i = 0; i < n; i++)
          fhilf[9][i] = fhilf[8][i] + *h * fhilf[11][i];
        x0 += *h;
        (*dgl)(x0, fhilf[9], fhilf[11]);

                                               /* use mid-point rule  */

        for (count = 1; count <= bufol[row] - 1; count++)
        {
          for (i = 0; i < n; i++)
            fhilf[10][i] = fhilf[8][i] + TWO * *h * fhilf[11][i];
          x0 += *h;
          (*dgl)(x0, fhilf[10], fhilf[11]);

          for (j = 8; j < 11; j++)             /* store for next step */
            copy_vector(fhilf[j], fhilf[j + 1], n);
        }

        (*dgl)(x0, fhilf[9], fhilf[11]);           /* stabilize with  */
        for (i = 0; i < n; i++)                    /* trapezoidal rule*/
          fhilf[row][i] = HALF * (fhilf[9][i] + fhilf[8][i] +
                                  *h * fhilf[11][i]);
        *aufrufe += bufol[row] + 2;
      }
      while (row == 0);

                                                     /* Extrapolation */

      extrapol(&row, fhilf, y, n, epsabs, epsrel, x, h, h0,
               &ahead, &index);
      if (*aufrufe >= fmax)                        /* too many calls ?*/
      {                                            /* => stop         */
        fhilf[9][0] = *x;
        copy_vector(fhilf[8], y, n);
        *x = x0;
        copy_vector(y, fhilf[index - 2], n);
        return 1;
      }
    }

    if (! ahead || ! neu)              /* do we need to repeat step ? */
    {
      neu = TRUE;                          /* set flag for a new call */
      if (row >= min(7, ext_max - 1))

      /* store differently since the extrapolation scheme has at most */
      /* 11 rows, but we allow only 8 extrapolations                  */

        for (j = 0; j < 7; j++)
          copy_vector(fhilf[j], fhilf[j + 1], n);

      /* accuracy could not be reached from the whole extrapolation   */
      /* table; repeat step with smaller step size                    */

      if (row >= ext_max + 2)
      {
        hilf = (REAL)(row - ext_max + 1);
        h0   = (REAL)0.9 * *h * POW((REAL)0.6, hilf);
        if (FABS(h0) <= FOUR * MACH_EPS * FABS(x0))
        {
          vmfree(vmblock);
          vmblock = NULL;
          return 2;
        }
        row = -1;
      }
    }
  }
}

/* -------------------------- END bulirsch.c ------------------------ */
