#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* -------------------------- MODULE glsp.c ------------------------- */

#include <basis.h>    /*  for  REAL, NULL, SQRT, FABS, ACOS, PI      */
#include <vmblock.h>
#include <u_proto.h>  /*  for  diag5pd, diag5                        */
#include <glsp.h>     /*  for  glspnp, glsppe, glsppa, glsptr,       */
                      /*       fzyfsy, fzyfsl                        */

static
int glsp1a  (int   n,     REAL* x,      REAL* f,
             REAL* w,     REAL  marg_0, REAL  marg_n,
             int   rep,   REAL* A,      REAL* B,
             REAL* C,     REAL* D,      REAL* h,
             REAL* h1,    REAL* h2,     REAL* md,
             REAL* ud1,   REAL* ud2);

static
int glsp2a  (int   n,     REAL* x,      REAL* f,
             REAL* w,     REAL  marg_0, REAL  marg_n,
             int   rep,   REAL* A,      REAL* B,
             REAL* C,     REAL* D,      REAL* h,
             REAL* h1,    REAL* h2,     REAL* md,
             REAL* ud1,   REAL* ud2);

static
int glsp3a  (int   n,     REAL* x,      REAL* f,
             REAL* w,     REAL  marg_0, REAL  marg_n,
             REAL* A,     REAL* B,      REAL* C,
             REAL* D,     REAL* ld1,    REAL* ld2,
             REAL* ud1,   REAL* ud2,
             REAL* h1,    REAL* h);



/*.BA*/
/*.FE{C 11.3}
     {Non-Parametric Cubic Fitting Splines}
     {Non-Parametric Cubic Fitting Splines}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int glspnp (int   n,
/*.IX{glspnp}*/
            REAL* xn, REAL* fn, REAL* w,
            int   marg_cond,
            REAL  marg_0, REAL marg_n,
            REAL* a, REAL* b, REAL* c, REAL* d
           )
/***********************************************************************
*  Compute the coefficients of a nonparametric cubic fitting spline.   *
*  The type of end point derivative is designated via  marg_cond.      *
*                                                                      *
*  The spline has the form :                                           *
*                                                                      *
*  s(x) = a[i] + b[i]*(x-xn[i]) + c[i]*(x-xn[i])^2 + d[i]*(x-xn[i])^3  *
*                                                                      *
*  for x in  [ xn[i], xn[i+1] ] ,   i = 0, ..., n-1 .                  *
.BE*)
*                                                                      *
*  Input parameters:                                                   *
*                                                                      *
*   int   n           Index of last node                               *
*                     we need :  n > 4 for marg_cond = 1,2,3 (s. below)*
*                                n > 5 for marg_cond = 4     (s. below)*
*   REAL  xn [n+1]    monotonically increasing x-values of nodes       *
*   REAL  fn [n+1]    y-values of nodes; if marg_cond=4 : fn[0] = fn[n]*
*   REAL  w  [n+1]    positive weights; for marg_cond=4 : w[0] = w[n]  *
*   int   marg_cond   type of end point condition :                    *
*                       = 1 : 1st derivatives at ends                  *
*                       = 2 : 2nd derivatives                          *
*                       = 3 : 3rd derivatives                          *
*                       = 4 : periodic spline                          *
*   REAL  marg_0    : end point derivatives at xn[0] and xn[n]         *
*   REAL  marg_n    : (not used for  marg_cond = 4)                    *
*                                                                      *
*                                                                      *
*  Output parameters:                                                  *
*                                                                      *
*   REAL  a [n+1]   : coefficients of spline in positions 0 to n-1     *
*   REAL  b [n+1]   : (the nth entry is used for aux purposes)         *
*   REAL  c [n+1]   :                                                  *
*   REAL  d [n+1]   :                                                  *
*                                                                      *
*                                                                      *
*  Return value :                                                      *
*                                                                      *
*   0 : no error                                                       *
*   1 : error in LU factorization in fdiasy, fdiag or fzyfsy           *
*   2 : n < 5 or n < 6                                                 *
*   3 : x-values not monotonic                                         *
*   4 : fn [0] differs from fn [n] or w [0] different from w [n]       *
*   5 : improper weight                                                *
*   6 : improper value for  marg_cond                                  *
*   7 : lack of memory                                                 *
*                                                                      *
*  subprograms used :        glsp1a, glsp2a, glsp3a, glsppe            *
*                                                                      *
*  Constant used :           NULL                                      *
*                                                                      *
*  REMARK : A natural fitting spline is computed when  marg_cond = 2   *
*           and marg_0 = marg_n = 0.0                                  *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
 int  error = 7, i;
 REAL *h   = NULL, *h1  = NULL, *h2 = NULL, *md  = NULL,
      *ud1 = NULL, *ud2 = NULL, *rs = NULL, *hup = NULL;
 void *vmblock;

 /* check input data */

 if (n < 5)                  return 2;
 for (i=0; i<=n-1; ++i)
   if (xn [i] >= xn [i+1])   return 3;
 for (i=0; i<=n; ++i)
   if (w [i] <= 0.0)         return 5;

 /* set up aux data  */

#define allo(x, n)  (x) = (REAL *)vmalloc(vmblock, VEKTOR, (n), 0)

 vmblock = vminit();

 switch (marg_cond)
 {
   case 1:
   case 2:
   case 3:  allo (h,  n);
            allo (h1, n);
            allo (h2, n);
            allo (md, n);
            allo (ud1,n);
            allo (ud2,n);
            break;
   case 4:  allo (h,   n+1);
            allo (h1,  n+1);
            allo (h2,  n+1);
            allo (md,  n+1);
            allo (rs,  n+1);
            allo (hup, 9*n-11);
            break;
   default: return 6;
 }
 if (! vmcomplete(vmblock))
 {
   vmfree(vmblock);
   return 7;
 }
#undef allo

 /* call subroutines, depending on type of end point conditions  */

 switch (marg_cond)
 {
   case 1: error = glsp1a (n, xn, fn, w, marg_0, marg_n, 0,
                           a, b, c, d, h, h1, h2, md, ud1, ud2);
           break;
   case 2: error = glsp2a (n, xn, fn, w, marg_0, marg_n, 0,
                           a, b, c, d, h, h1, h2, md, ud1, ud2);
           break;
   case 3: error = glsp3a (n, xn, fn, w, marg_0, marg_n, a, b, c, d,
                           h2, md, ud1, ud2, h, h1);
           break;
   case 4: if (n < 6)
             error = -1;
           else
             error = glsppe (n, xn, fn, w, 0, a, b, c, d,
                             h, h1, h2, md, rs, hup);
           break;
 }
 vmfree(vmblock);
 return error;
}

static REAL hf[10];   /* aux vector that may not be altered after     */
                      /* the first call of a run                      */

static
int glsp1a (int   n,
/*.IX{glsp1a}*/
            REAL* xn, REAL* fn, REAL* w,
            REAL  marg_0, REAL marg_n,
            int   rep,
            REAL* a, REAL* b, REAL* c, REAL* d,
            REAL* h, REAL* h1, REAL* h2,
            REAL* md, REAL* ud1, REAL* ud2
           )
/***********************************************************************
*  Compute the coefficients of a cubic fitting spline for given first  *
*  derivatives at the end points.                                      *
*                                                                      *
*  The spline has the form :                                           *
*                                                                      *
*  s(x) = a[i] + b[i]*(x-xn[i]) + c[i]*(x-xn[i])^2 + d[i]*(x-xn[i])^3  *
*                                                                      *
*  for x in  [ xn[i], xn[i+1] ] ,   i = 0, ..., n-1 .                  *
*                                                                      *
*  Input parameters:                                                   *
*                                                                      *
*   int   n           Index of final node : n > 4                      *
*   REAL  xn [n+1]    x-values for nodes, strictly increasing          *
*   REAL  fn [n+1]    y-values of nodes                                *
*   REAL  w  [n+1]    positive weights                                 *
*   REAL  marg_0      1st derivative at  xn[0]                         *
*   REAL  marg_n      1st derivative at  xn[n]                         *
*   int   rep         indicates repeated call:                         *
*                       0 : first call, form matrix for the c[i] and   *
*                           use fdiasy to factor                       *
*                       1 : repeated call, form new right hand side    *
*                           only; here the vectors md,ud1, ud2, h, h1, *
*                           h2 with the factorization data must not be *
*                           altered.                                   *
*                                                                      *
*  Output parameters:                                                  *
*                                                                      *
*   REAL  a [n+1]   : coefficients of the spline                       *
*   REAL  b [n+1]   : in positions  0 to n-1;                          *
*   REAL  c [n+1]   : (position n is used for aux purposes only)       *
*   REAL  d [n+1]   :                                                  *
*                                                                      *
*  aux vectors:                                                        *
*                                                                      *
*   REAL  h  [n]    :                                                  *
*   REAL  h1 [n]    :                                                  *
*   REAL  h2 [n]    :                                                  *
*   REAL  md [n]    : (position zero in  md, ud1 and ud2 not used)     *
*   REAL  ud1[n]    :                                                  *
*   REAL  ud2[n]    :                                                  *
*                                                                      *
*                                                                      *
*  Return value:                                                       *
*                                                                      *
*   0 : no error                                                       *
*   1 : error in factorization  diag5pd                                *
*   2 : n < 5                                                          *
*   3 : improper value for  rep                                        *
*                                                                      *
*  subprograms used:        diag5pd                                    *
*                                                                      *
*  REMARKS :   (i) glsp1a  should only be called from a routine that   *
*                  checks the input, such as  glspnp or glsppa;        *
*             (ii) For parametric splines with non constant weights    *
*                  rep must equal 0.                                   *
*                                                                      *
***********************************************************************/
{
  int i, k, error;
  REAL   h_var_1, h_var_2;

  if (rep != 0 && rep != 1)   return (3);
  if (!rep)
 /*
     First call: determine aux values and LU factorization
 */
  {
    for (i=0; i<=n-1; ++i)
    {
      h [i]  = xn [i+1] - xn [i];
      h1 [i] = 1. / h [i];
      c [i]  = h1 [i] * h1 [i];
      b [i]  = 6. / w [i];
    }
    b [n] = 6. / w [n];
    for (i=0; i<=n-2; ++i)
      h2 [i] = h1 [i] + h1 [i+1];
 /*
    second  co-diagonal
 */
    for (i=1; i<=n-3; ++i)
      ud2 [i] = b [i+1] * h1 [i] * h1 [i+1];
 /*
    first co-diagonal
 */
    for (i=1; i<=n-2; ++i)
      ud1 [i] = h [i] - b [i] * h1 [i] * h2 [i-1]
                      - b [i+1] * h1 [i] * h2 [i];
 /*
    main diagonal
 */
    for (i=1; i<=n-1; ++i)
    {
      k = i - 1;
      md [i] = 2. * (h [k] + h [i]) + b [k] * c [k]
             + b [i] * h2 [k] * h2 [k] + b [i+1] * c [i];
    }
 /*
    The global aux vector hf is used to alter the corners in the system
    matrix. Its 10 elements must be available for repeated calls.
 */
    hf [ 0 ] = h [ 0 ] - b [ 0 ] * c [ 0 ]  - b [ 1 ] * h2[ 0] * h1[ 0];
    hf [ 1 ] = b [ 1 ] * h1[ 0 ] * h1[ 1 ];
    hf [ 2 ] = b [n-1] * h1[n-2] * h1[n-1];
    hf [ 3 ] = h [n-1] - b [n-1] * h2[n-2]  * h1[n-1] - b[ n ] * c[n-1];
    hf [ 4 ] = h1[ 0 ] * (b[ 1 ] + b [ 0 ]) +   2.    * h[ 0 ] * h[ 0 ];
    hf [ 5 ] = h1[n-1] * (b[ n ] + b [n-1]) +   2.    * h[n-1] * h[n-1];
    hf [ 6 ] = b [ 1 ] * h2[ 0 ] + b [ 0 ]  * h1[ 0 ] - h[ 0 ] * h[ 0 ];
    hf [ 7 ] = b [n-1] * h2[n-2] + b [ n ]  * h1[n-1] - h[n-1] * h[n-1];
    hf [ 8 ] = b [ 1 ] * h1[ 1 ];
    hf [ 9 ] = b [n-1] * h1[n-2];

    md  [ 1 ] += hf [0] / hf [4] * hf [6];
    ud1 [ 1 ] -= hf [0] / hf [4] * hf [8];
    md  [ 2 ] -= hf [1] / hf [4] * hf [8];
    md  [n-2] -= hf [2] / hf [5] * hf [9];
    ud1 [n-2] += hf [2] / hf [5] * hf [7];
    md  [n-1] += hf [3] / hf [5] * hf [7];
 }
 /*
     Compute right hand side
 */
  h_var_1 = (fn[1] - fn[0]) * h1[0];
  for (i = 1; i <= n-1; ++i, h_var_1 = h_var_2)
    {
      h_var_2 = (fn[i+1] - fn[i]) * h1[i];
      c [i]   = 3. * (h_var_2 - h_var_1);
    }

  h_var_1 = 3. * ((fn[1] - fn[0]) - marg_0 * h[0]);
  h_var_2 = 3. * (marg_n * h[n-1] - (fn[n] - fn[n-1]));

  c [ 1 ] -= hf [0] / hf [4] * h_var_1;
  c [ 2 ] -= hf [1] / hf [4] * h_var_1;
  c [n-2] -= hf [2] / hf [5] * h_var_2;
  c [n-1] -= hf [3] / hf [5] * h_var_2;
 /*
     Compute the  c[i], i=1, ..., n-1
     by solving the linear system via diag5pd
 */
  if (! rep)
  {
    error = diag5pd(0, n - 1, md + 1, ud1 + 1,      /* for first call */
                    ud2 + 1, c + 1);
    if (error != 0)
      if (error == 1)
        return 2;
      else
        return 1;
  }
  else
    diag5pd(2, n - 1, md + 1, ud1 + 1,          /* for repeated call  */
            ud2 + 1, c + 1);

 /*
     compute remaining spline coefficients
 */
  c [0] = (h_var_1 + c[ 1 ] * hf[6] - c[ 2 ] * hf[8]) / hf[4];
  c [n] = (h_var_2 + c[n-1] * hf[7] - c[n-2] * hf[9]) / hf[5];

  a [0] = fn[0] + 2. / w[0] * h1[0] * (c[0] - c[1]);
  for (i=1; i<=n-1; ++i)
  {
    k = i - 1;
    a [i] = fn[i] - 2. / w[i] * (c[k] * h1[k] - h2[k] * c[i]
            + c[i+1] * h1[i]);
  }
  a [n] = fn [n] - 2. / w[n] * h1[n-1] * (c[n-1] - c[n]);

  for (i=0; i<=n-1; ++i)
  {
    k = i + 1;
    b [i] = h1[i] * (a[k] - a[i]) - h[i] / 3. * (c[k] + 2. * c[i]);
    d [i] = h1[i] / 3. * (c[k] - c[i]);
  }

  return 0;
}

static
int glsp2a (int   n,
/*.IX{glsp2a}*/
            REAL* xn, REAL* fn, REAL* w,
            REAL  marg_0, REAL marg_n,
            int   rep,
            REAL* a, REAL* b, REAL* c, REAL* d,
            REAL* h, REAL* h1, REAL* h2,
            REAL* md, REAL* ud1, REAL* ud2
           )
/***********************************************************************
*  Compute the coefficients of a cubic fitting spline for given second *
*  derivatives at the end points.                                      *
*                                                                      *
*  The spline has the form:                                            *
*                                                                      *
*  s(x) = a[i] + b[i]*(x-xn[i]) + c[i]*(x-xn[i])^2 + d[i]*(x-xn[i])^3  *
*                                                                      *
*  for x in  [ xn[i], xn[i+1] ] ,   i = 0, ..., n-1 .                  *
*                                                                      *
*  Input parameters:                                                   *
*                                                                      *
*   int   n           Index of final node : n > 4                      *
*   REAL  xn [n+1]    x-values for nodes, strictly increasing          *
*   REAL  fn [n+1]    y-values of nodes                                *
*   REAL  w  [n+1]    positive weights                                 *
*   REAL  marg_0      2nd derivative at  xn[0]                         *
*   REAL  marg_n      2nd derivative at  xn[n]                         *
*   int   rep         indicates repeated call:                         *
*                       0 : first call, form matrix for the c[i] and   *
*                           use fdiasy to factor                       *
*                       1 : repeated call, form new right hand side    *
*                           only; here the vectors md, ud1, ud2, h, h1,*
*                           h2 with the factorization data must not be *
*                           altered.                                   *
*                                                                      *
*  Output parameters:                                                  *
*                                                                      *
*   REAL  a [n+1]   : coefficients of the spline                       *
*   REAL  b [n+1]   : in positions  0 to n-1;                          *
*   REAL  c [n+1]   : (position n is only used for aux purposes)       *
*   REAL  d [n+1]   :                                                  *
*                                                                      *
*  aux vectors :                                                       *
*                                                                      *
*   REAL  h  [n]    :                                                  *
*   REAL  h1 [n]    :                                                  *
*   REAL  h2 [n]    :                                                  *
*   REAL  md [n]    : (position zero in  md, ud1 and ud2 not used)     *
*   REAL  ud1[n]    :                                                  *
*   REAL  ud2[n]    :                                                  *
*                                                                      *
*                                                                      *
*  Return value :                                                      *
*                                                                      *
*   0 : no error                                                       *
*   1 : error in  diag5pd                                              *
*   2 : n < 5                                                          *
*   3 : improper value for  rep                                        *
*                                                                      *
*  subprograms used:         diag5pd                                   *
*                                                                      *
*  REMARKS :   (i) glsp2a  should only be called from a routine that   *
*                  checks the input, such as  glspnp or glsppa;        *
*             (ii) For parametric splines with non constant weights    *
*                  rep must equal 0.                                   *
*            (iii) Fa natural fitting spline set marg_0 = marg_n = 0.0.*
*                                                                      *
***********************************************************************/
{
  int i, k, error;
  REAL   h_var_1, h_var_2;

  if (rep != 0 && rep != 1)   return (3);
  if (!rep)
 /*  Compute aux values and three diagonals of linear system          */
  {
    for (i=0; i<=n-1; ++i)
    {
      h [i]  = xn [i+1] - xn [i];
      h1 [i] = 1. / h [i];
      c [i]  = h1 [i] * h1 [i];
      b [i]  = 6. / w [i];
    }
    b [n] = 6. / w [n];

    for (i=0; i<=n-2; ++i)
      h2 [i] = h1 [i] + h1 [i+1];
 /*
    second co-diagonal
 */
    for (i=1; i<=n-3; ++i)
      ud2 [i] = b [i+1] * h1 [i] * h1 [i+1];
 /*
    first co-diagonal
 */
    for (i=1; i<=n-2; ++i)
      ud1 [i] = h [i] - b [i] * h1 [i] * h2 [i-1]
                      - b [i+1] * h1 [i] * h2 [i];
 /*
    main diagonal
 */
    for (i=1; i<=n-1; ++i)
    {
      k = i - 1;
      md [i] = 2. * (h [k] + h [i]) + b [k] * c [k]
             + b [i] * h2 [k] * h2 [k] + b [i+1] * c [i];
    }
  }
 /*
     Compute right hand side
 */
  c [0] = 0.5 * marg_0;
  c [n] = 0.5 * marg_n;

  h_var_2 = (fn[2] - fn[1]) * h1[1];
  h_var_1 = (fn[3] - fn[2]) * h1[2];
  c[1] = 3. * (h_var_2 - (fn[1] - fn[0]) * h1[0]) - c[0] *
         (h[0] - 6. / w[0] * h1[0] * h1[0] - 6. / w[1] * h1[0] * h2[0]);
  c[2] = 3. * (h_var_1 - h_var_2) - c[0] * (6. / w[1]) * h1[0] * h1[1];
  for (i=3; i<=n-3; ++i, h_var_1=h_var_2)
  {
    h_var_2 = (fn[i+1] - fn[i]) * h1[i];
    c [i]   = 3. * (h_var_2 - h_var_1);
  }
  h_var_2 = (fn[n-1] - fn[n-2]) * h1[n-2];
  c [n-2]  = 3. * (h_var_2 - h_var_1)
             - c[n] * 6. / w[n-1] * h1[n-2] * h1[n-1];
  c [n-1]  = 3. * ((fn[n] - fn[n-1]) * h1[n-1] - h_var_2)
             - c[n] * (h[n-1] - 6. / w[n-1] * h1[n-1] * h2[n-2]
                              - 6. / w[n] * c[n-1]);
 /*
     compute coefficients  c[i], i=1, ..., n-1
     via LU factorization in  diag5pd
 */
  if (! rep)
  {
    error = diag5pd(0, n - 1, md + 1, ud1 + 1,    /*  first call      */
                    ud2 + 1, c + 1);
    if (error != 0)
      if (error == 1)
        return 2;
      else
        return 1;
  }
  else
    diag5pd(2, n - 1, md + 1, ud1 + 1,            /*  repeated call   */
            ud2 + 1, c + 1);

 /*
     compute remaining spline coefficients
 */
  a [0] = fn[0] + 2. / w[0] * h1[0] * (c[0] - c[1]);
  for (i=1; i<=n-1; ++i)
  {
    k = i - 1;
    a [i] = fn[i] - 2. / w[i] * (c[k] * h1[k] - h2[k] * c[i]
          + c[i+1] * h1[i]);
  }
  a [n] = fn [n] - 2. / w[n] * h1[n-1] * (c[n-1] - c[n]);

  for (i=0; i<=n-1; ++i)
  {
    k = i + 1;
    b [i] = h1[i] * (a[k] - a[i]) - h[i] / 3. * (c[k] + 2. * c[i]);
    d [i] = h1[i] / 3. * (c[k] - c[i]);
  }

  return 0;
}

static
int glsp3a (int   n,
/*.IX{glsp3a}*/
            REAL* xn, REAL* fn, REAL* w,
            REAL  marg_0, REAL marg_n,
            REAL* a, REAL* b, REAL* c, REAL* d,
            REAL* ld1, REAL* ld2, REAL* ud1, REAL* ud2,
            REAL* h1, REAL* h
           )
/***********************************************************************
*                                                                      *
*  Compute the coefficients of a cubic fitting spline for given third  *
*  derivatives at the end points.                                      *
*                                                                      *
*  The spline has the form:                                            *
*                                                                      *
*  s(x) = a[i] + b[i]*(x-xn[i]) + c[i]*(x-xn[i])^2 + d[i]*(x-xn[i])^3  *
*                                                                      *
*  for x in  [ xn[i], xn[i+1] ] ,   i = 0, ..., n-1 .                  *
*                                                                      *
*  Input parameters:                                                   *
*                                                                      *
*   int   n           Index of final node : n > 4                      *
*   REAL  xn [n+1]    x-values for nodes, strictly increasing          *
*   REAL  fn [n+1]    y-values of nodes                                *
*   REAL  w  [n+1]    positive weights                                 *
*   REAL  marg_0      3rd derivative at  xn[0]                         *
*   REAL  marg_n      3rd derivative at  xn[n]                         *
*                                                                      *
*  Output parameters:                                                  *
*                                                                      *
*   REAL  a [n+1]   : coefficients of the spline                       *
*   REAL  b [n+1]   : in positions  0 to n-1;                          *
*   REAL  c [n+1]   : (position n is only used for aux purposes)       *
*   REAL  d [n+1]   :                                                  *
*                                                                      *
*  aux vectors:                                                        *
*                                                                      *
*   REAL  h  [n]    :                                                  *
*   REAL  h1 [n]    :                                                  *
*   REAL  ld1[n]    :                                                  *
*   REAL  ld2[n]    : (position zero in ld1, ld2, ud1, ud2 not used)   *
*   REAL  ud1[n]    :                                                  *
*   REAL  ud2[n]    :                                                  *
*                                                                      *
*  Return value:                                                       *
*                                                                      *
*    0 : no error                                                      *
*    1 : error in  diag5                                               *
*    2 : n < 5                                                         *
*                                                                      *
*  subprograms used:         diag5                                     *
*                                                                      *
*  REMARK :  glsp3a  should only be called from a routine that checks  *
*            the input, such as  glspnp.                               *
*                                                                      *
***********************************************************************/
{
  int i, k, error;
  REAL   h_var_1, h_var_2;
 /*
     compute aux values
 */
  for (i=0; i<=n-1; ++i)
  {
    h [i]  = xn [i+1] - xn [i];
    h1 [i] = 1. / h [i];
    c [i]  = h1 [i] * h1 [i];
    b [i]  = 6. / w [i];
  }
  b [n] = 6. / w [n];

  for (i=0; i<=n-2; ++i)
    d [i] = h1 [i] + h1 [i+1];
 /*
     Compute five diagonal system matrix A and the right hand side RS
     of the system  A * C = RS
 */
  for (i=3; i<=n-1; ++i)
  {
    ld2 [i]   = b [i-1] * h1 [i-2] * h1 [i-1];  /* 2nd co-diagonals   */
    ud2 [i-2] = ld2 [i];
  }

  h_var_1 = h[1] - b[2] * h1 [1] * d[1];        /* first co-diagonals */
  ld1 [2] = h_var_1 - b[1] * c[1];
  ud1 [1] = h_var_1 - b[1] * h1[1] * d[0];
  for (i=3; i<=n-2; ++i)
    {k = i - 1;
     ld1 [i] = h [k] - b [k] * h1[k] * d[k-1]
                     - b [i] * h1 [k] * d [k];
     ud1 [k] = ld1 [i];
    }
  h_var_1   = h[n-2] - b[n-2] * h1 [n-2] * d[n-3];
  ld1 [n-1] = h_var_1 - b[n-1] * h1[n-2] * d[n-2];
  ud1 [n-2] = h_var_1 - b[n-1] * c[n-2];

  a [1] = 3. * h[0] + 2. * h[1] + b[1] * h1[1]*d[0] + b[2]*c[1];
  for (i=2; i<=n-2; ++i)
  {
    k = i - 1;                                      /*  main diagonal */
    a [i] = 2. * (h [k] + h [i]) + b [k] * c [k]
          + b [i] * d [k] * d [k] + b [i+1] * c [i];
  }
  a [n-1] = 3. * h[n-1] + 2. * h[n-2] + b[n-1] * h1[n-2]*d[n-2]
               + b[n-2] * c[n-2];
  c [0] = 0.5 * marg_0;                           /* right hand side  */
  c [n] = 0.5 * marg_n;

  h_var_2 = (fn[2] - fn[1]) * h1[1];
  h_var_1 = (fn[3] - fn[2]) * h1[2];
  c [1]  = 3. * (h_var_2 - (fn[1] - fn[0]) * h1[0])
           + c[0] * (h[0]*h[0] - b[0] * h1[0] - b[1] * d[0]);
  c [2]  = 3. * (h_var_1 - h_var_2) + c[0] * b[1] * h1[1];
  for (i=3; i<=n-3; ++i, h_var_1=h_var_2)
  {
    h_var_2 = (fn[i+1] - fn[i]) * h1[i];
    c [i]   = 3. * (h_var_2 - h_var_1);
  }
  h_var_2 = (fn[n-1] - fn[n-2]) * h1[n-2];
  c [n-2] = 3. * (h_var_2 - h_var_1)
            - c[n] * b[n-1] * h1[n-2];
  c [n-1] = 3. * ((fn[n] - fn[n-1]) * h1[n-1] - h_var_2)
            - c[n] * (h[n-1]*h[n-1] - b[n-1] * d[n-2] - b[n] * h1[n-1]);
 /*
     Compute coefficients c[i], i=1, ..., n-1
     by LU factorization in diag5
 */
  error = diag5(0, n - 1, ld2 + 1, ld1 + 1, a + 1, ud1 + 1, ud2 + 1,
                c + 1);
  if (error != 0)
    if (error == 2)
      return 1;
    else
      return 2;

  c [0] = c[1] - c[0] * h[0];
  c [n] = c[n-1] + c[n] * h[n-1];
 /*
     Compute remaining spline coefficients
 */
  a [0] = fn[0] + b[0] / 3. * h1[0] * (c[0] - c[1]);
  for (i=1; i<=n-1; ++i)
  {
    k = i - 1;
    a [i] = fn[i] - b[i] / 3. * (c[k] * h1[k] - d[k] * c[i]
          + c[i+1] * h1[i]);
  }
  a [n] = fn [n] - b[n] / 3. * h1[n-1] * (c[n-1] - c[n]);

  for (i=0; i<=n-1; ++i)
  {
    k = i + 1;
    b [i] = h1[i] * (a[k] - a[i]) - h[i] / 3. * (c[k] + 2. * c[i]);
    d [i] = h1[i] / 3. * (c[k] - c[i]);
  }

  return 0;
}
/*.BA*/



/*.FE{C 11.4}
     {Parametric Cubic Fitting Splines}
     {Parametric Cubic Fitting Splines}*/

/*.BE*/
/* ------------------------------------------------------------------ */
/*.BA*/

int glsppa (int   n,
/*.IX{glsppa}*/
            REAL* xn, REAL* fn,
            REAL* wx, REAL* wf,
            REAL* t,
            int   marke_t,
            int   rand,
            REAL* alpha, REAL* beta,
            int   marke_w,
            REAL* ax, REAL* bx, REAL* cx, REAL* dx,
            REAL* ay, REAL* by, REAL* cy, REAL* dy,
            REAL* help
           )
/***********************************************************************
* Compute the coefficients                                             *
*      ax[i], bx[i], cx[i] and dx[i]  as well as                       *
*      ay[i], by[i], cy[i] and dy[i], i = 0, ..., n-1                  *
* for a parametric cubic fitting spline for varying end point          *
* conditions determined in rand.                                       *
*                                                                      *
* The parametric spline with parameter  t[i], i=0, ..., n is composed  *
* of two components  sx and sy in the following form with              *
* u := (t - t[i]) :                                                    *
*                                                                      *
*  sx := sx(t) = ax[i] + bx[i] * u + cx[i] * u^2 + dx[i] * u^3         *
*  sy := sy(t) = ay[i] + by[i] * u + cy[i] * u^2 + dy[i] * u^3         *
*                                                                      *
* for t in the interval [t[i],t[i+1]], i=0, ..., n-1.                  *
.BE*)
*                                                                      *
* sx and sy are non-parametric cubic splines and are computed via the  *
* algorithms of chapter  11.3 .                                        *
*                                                                      *
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*   int  n             index of last node                              *
*   REAL xn[n+1]       nodes : x-values                                *
*   REAL fn[n+1]               y-values at  xn[]                       *
*   REAL wx[n+1]       weights for xn[i]                               *
*   REAL wf[n+1]       weights for fn[i]                               *
*   REAL t [n+1]       parameter values for (xn[i], fn[i]) if marke_t=1*
*   int  marke_t       Control for curve parameter t input:            *
*                        1 : user supplies t values                    *
*                        0 : no values for t[i] supplied; they shall   *
*                            be set up using the formula :             *
*                              t[ 0 ] = 0                              *
*                              t[i+1] = t[i] + sqrt((xn[i+1]-xn[i])^2  *
*                                                +  (fn[i+1]-fn[i])^2) *
*   int  rand          type of end point condition :                   *
*                        1 : 1st end point derivative wrt t given (n>4)*
*                        2 : 2nd derivative                       (n>4)*
*                        3 : derivative dy/dx given at end points (n>4)*
*                        4 : periodic spline; n > 5 and  xn[0]=xn[n],  *
*                            wx[0]=wx[n], fn[0]=fn[n] and wf[0]=wf[n]. *
*   REAL   alpha []  : if rand=1 or rand=2:                            *
*   REAK   beta  []  :   rand-th derivative wrt. t                     *
*                          alpha[1] = sx[rand]*t[0]                    *
*                          alpha[2] = sy[rand]*t[0]                    *
*                          beta [1] = sx[rand]*t[n]                    *
*                          beta [2] = sy[rand]*t[n]                    *
*                      if rand=3:                                      *
*                        1st end point derivatives dy/dx               *
*                          alpha[1] = dy/dx*xn[0]                      *
*                          alpha[2] = not used                         *
*                          beta [1] = dy/dx*xn[n]                      *
*                          beta [2] = not used                         *
*                        If the magnitude of alpha[1] or beta[1]       *
*                        exceeds  1E10), the tangent vector is instead *
*                        computed as:                                  *
*                                . .                                   *
*                               (x,y) = (0, sign (1, fn[1]-fn[0]))     *
*                                          for left end point          *
*                                . .                                   *
*                               (x,y) = (0, sign (1, fn[n]-fn[n-1]))   *
*                                         for right end point          *
*                      rand=4: alpha, beta vectors not used            *
*                                                                      *
*                      A natural parametric fitting spline results from*
*                      rand=2 and alpha[1]=alpha[2]=beta[1]=beta[2]=0. *
*                                                                      *
*     int marke_w      Control for type of weights wx, wf :            *
*                        0 : use weights  wx[i] = wf[i] = 1            *
*                        1 : user supplies weights  wf; program sets   *
*                            wx[i] = wf[i]                             *
*                        2 : user supplies both wx[i] and wf[i]        *
*     REAL help [14*n-9]                                               *
*                      aux vector                                      *
*                                                                      *
* Output parameters:                                                   *
*                                                                      *
*   REAL ax [n+1]   spline coefficients for sx                         *
*   REAL bx [n+1]                                                      *
*   REAL cx [n+1]                                                      *
*   REAL dx [n+1]                                                      *
*   REAL ay [n+1]   ditto for sy                                       *
*   REAL by [n+1]                                                      *
*   REAL cy [n+1]                                                      *
*   REAL dy [n+1]                                                      *
*   REAL  t [n+1]   output for marke_t = 0: the parameter values from  *
*                   glsppa                                             *
*                                                                      *
*  Return value :                                                      *
*                                                                      *
*     0 : all is ok                                                    *
*    -1 : n too small                                                  *
*    -2 : improper values for rand or marke_w                          *
*    -3 : improper weights wx or wf                                    *
*    -4 : parameter values not monotonically increasing (t[i]>=t[i+1]) *
*    -5 : for rand=4:         xn[0] != xn[n]                           *
*                      or     fn[0] != fn[n]                           *
*                      or     wx[0] != wx[n]                           *
*                      or     wf[0] != wf[n]                           *
*     1 : error in  diag5pd or fzyfsy                                  *
*                                                                      *
* subroutines used :                                                   *
*                                                                      *
*    glsp1a : cubic fitting splines with first end point derivatives   *
*    glsp2a : ditto for second derivatives                             *
*    glsppe : periodic cubic fitting splines                           *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{ int i, mrp, error;
  REAL og, alphax, alphay, betax, betay, delta, deltx, delty, root;

 if (n < 5) return -1;                            /* check input data */
 if (rand<1  ||  rand>4 || marke_w<0 || marke_w>2)  return (-2);

 if (marke_w == 2)                                   /* check weights */
 {
   for (i=0; i<n+1; i++)
     if (wx[i] <= 0.0  ||  wf[i] <= 0.0)  return -3;
   mrp = 0;
 }
 else
 {
   if (marke_w == 1)
   {
     for (i=0; i<n+1; i++)
       if (wf[i] <= 0.0)  return -3;
   }
   else                                  /* assign weights = 1 for wf */
     for (i=0; i<n+1; i++)
       wf[i] = 1.0;
   for (i=0; i<n+1; i++)
     wx[i] = wf[i];
   mrp = 1;
 }

 if (marke_t == 0)
   for (i=1,t[0]=0.0; i<n+1; i++)   /* compute parameter values for t */
   {
     deltx = xn[i] - xn[i-1];
     delty = fn[i] - fn[i-1];
     delta = deltx*deltx + delty*delty;
     if (delta <= 0.0)  return -4;
     t[i] = t[i-1] + SQRT(delta);
   }
 else
   for (i=0; i<n; i++)                      /* check parameter values */
     if (t[i+1] <= t[i])  return -4;


                                       /* Compute spline coefficients */
 if (rand == 1)              /* 1st end point derivatives wrt t given */
 {
   error = glsp1a (n,t,xn,wx,alpha[1],beta[1],0,ax,bx,cx,dx,&help[1],
                   &help[n+1],&help[2*n+1],&help[3*n],&help[4*n-1],
                   &help[5*n-2]);
   if (error != 0)  return (error);
   error = glsp1a (n,t,fn,wf,alpha[2],beta[2],mrp,ay,by,cy,dy,&help[1],
                   &help[n+1],&help[2*n+1],&help[3*n],&help[4*n-1],
                   &help[5*n-2]);
 }
 else if (rand == 2)         /* 2nd end point derivatives wrt t given */
 {
   error = glsp2a (n,t,xn,wx,alpha[1],beta[1],0,ax,bx,cx,dx,&help[1],
                   &help[n+1],&help[2*n+1],&help[3*n],&help[4*n-1],
                   &help[5*n-2]);
   if (error != 0)  return (error);
   error = glsp2a (n,t,fn,wf,alpha[2],beta[2],mrp,ay,by,cy,dy,&help[1],
           &help[n+1],&help[2*n+1],&help[3*n],&help[4*n-1],
           &help[5*n-2]);
 }

 else if (rand == 3)          /* 1st derivatives dy/dx given at ends */
 {
   og = 1.E10;
   if ( FABS(alpha[1]) >= og )
   {
     alphax = 0.0;
     alphay = (fn[1]-fn[0] >= 0.0)? 1.0:-1.0;
   }
   else
   {
     root = SQRT ( 1.0/(1.0+alpha[1]*alpha[1]) );
     alphax = (xn[1]-xn[0] >= 0.0)?  root : -root;
     alphay = alphax*alpha[1];
   }
   if ( FABS(beta[1]) >= og )
   {
     betax = 0.0;
     betay = (fn[n]-fn[n-1] >= 0.0)?  1.0 : -1.0;
   }
   else
   {
     root = SQRT ( 1.0/(1.0+beta[1]*beta[1]) );
     betax = (xn[n]-xn[n-1] >= 0.0)?  root : -root;
     betay = betax*beta[1];
   }
   error = glsp1a (n,t,xn,wx,alphax,betax,0,ax,bx,cx,dx,&help[1],
                   &help[n+1],&help[2*n+1],&help[3*n],&help[4*n-1],
                   &help[5*n-2]);
   if (error != 0)  return error;
   error = glsp1a (n,t,fn,wf,alphay,betay,mrp,ay,by,cy,dy,&help[1],
                   &help[n+1],&help[2*n+1],&help[3*n],&help[4*n-1],
                   &help[5*n-2]);
 }

 else                                              /* periodic spline */
 {
   if (n < 6)  return -1;
   error = glsppe (n,t,xn,wx,0,ax,bx,cx,dx,&help[1],&help[n+2],
                 &help[2*n+3],&help[3*n+4],&help[4*n+4],&help[5*n+4]);
   if (error != 0)  return error;
   error = glsppe (n,t,fn,wf,mrp,ay,by,cy,dy,&help[1],&help[n+2],
                 &help[2*n+3],&help[3*n+4],&help[4*n+4],&help[5*n+4]);
 }
 return error;
}

/* ------------------------------------------------------------------ */
/*.BA*/

int glsppe (int   n,
/*.IX{glsppe}*/
            REAL* xn, REAL* fn, REAL* w,
            int   rep,
            REAL* a, REAL* b, REAL* c, REAL* d,
            REAL* h, REAL* h1, REAL* h2, REAL* h3,
            REAL* rs, REAL* hup
           )
/***********************************************************************
*  Computes the coefficients of a periodic cubic fitting spline.       *
*                                                                      *
*  The spline has the representation :                                 *
*                                                                      *
*  S := S(X) = A(i) + B(i)*(X-XN(i)) + C(i)*(X-XN(i))**2               *
*                         + D(i)*(X-XN(i))**3                          *
*                                                                      *
*  for X in the interval  [ XN(i) , XN(i+1) ] ,   i = 0, ..., n-1 .    *
.BE*)
*                                                                      *
* ==================================================================== *
*                                                                      *
*  Input parameters:                                                   *
*  -----------------                                                   *
*                                                                      *
*   Name  Type/size    Meaning                                         *
*  ------------------------------------------------------------------- *
*   int  n             Number of last node (n > 5)                     *
*   REAL xn[n+1]       monotonically increasing x-values of nodes      *
*   REAL fn[n+1]       y-values with  fn[0] = fn[n]                    *
*   REAL w [n+1]       positive weights with  w[0] = w[n]              *
*   int  rep           ontrol parameter for repeated calls             *
*                        0 : first call for which the system matrix    *
*                            is formed and the system solved for the   *
*                            c(i) using  fzyfsy .                      *
*                        1 : reated call; form new right hand side only*
*                            The vectors h, h1, h2, h3 and hup are used*
*                            unaltered from the first call.            *
*                                                                      *
*                                                                      *
*  Output parameters:                                                  *
*  ------------------                                                  *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*   a       REAL  /[n+1]        }  coefficients ot the spline in       *
*   b       REAL  /[n+1]        }  positions   0 to n-1                *
*   c       REAL  /[n+1]        }  ( the last element is used for aux  *
*   d       REAL  /[n+1]        }    purposes only )                   *
*                                                                      *
*                                                                      *
*  Aux variables:                                                      *
*  --------------                                                      *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*   h       REAL  /[n+1]        }                                      *
*   h1      REAL  /[n+1]        }  aux vectors                         *
*   h2      REAL  /[n+1]        }                                      *
*   h3      REAL  /[n+1]        }  (the 0th elements in rs and hup are *
*   rs      REAL  /[n+1]        }   not used)                          *
*   hup     REAL  /[9*n-11]     }                                      *
*                                                                      *
*                                                                      *
*  Return value :                                                      *
*  --------------                                                      *
*                                                                      *
*   = 0 : no error  r                                                  *
*   = 1 : error in  fzyfsy                                             *
*   = 2 : n < 6                                                        *
*   = 3 : improper value for rep                                       *
*   = 4 : fn[0] different from fn[n] or w[0] different from w[n]       *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   subprograms used :        fzyfsy, fzyfsz, fzyfsl                   *
*   ------------------                                                 *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   REMARK : (i)  glsppe  shpuld only be called from a program that    *
*   -------       checks the assumptions for the input data such as    *
*                 glspnp or glsppa;                                    *
*           (ii)  For parametric splines with non constant weights, we *
*                 must have  rep = 0 .                                 *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{ int i, k, error;
  REAL   h_var_1, h_var_2;

  if (rep != 0 && rep != 1)              return 3;
 /*
     Check periodicity
 */
  if (fn[n] != fn[0]  ||  w[n] != w[0])  return 4;

  if (!rep)
 /*
     First call : i.e., we must determine the aux variables and the
     entries of the symmetric almost cyclic five-diagonal system matrix.
 */
  {
    for (i=0; i<=n-1; ++i)                         /*  aux variables  */
    {
      h [i]  = xn [i+1] - xn [i];
      h1 [i] = 1. / h [i];
      c [i]  = h1 [i] * h1 [i];
      h2 [i]  = 6. / w [i];
    }
    h [n]  = h [0];  h1 [n] = h1 [0];  c [n]  = c [0];  h2 [n] = h2 [0];
    for (i=0; i<=n-1; ++i)
      h3 [i] = h1[i] + h1[i+1];

    for (i=1; i<=n-1; ++i)                    /*  second co-diagonal  */
      d [i]   = h2[i+1] * h1[i] * h1[i+1];
    d [n] = h2[1] * h1[0] * h1[1];

    for (i=1; i<=n-1; ++i)                      /* first co-diagonal  */
      b [i] = h[i] - h2[i] * h1[i] * h3[i-1] - h2[i+1] * h1[i] * h3[i];
    b [n] = h[0] - h2[0] * h1[0] * h3[n-1] - h2[1] * h1[0] * h3[0];

    for (i=1; i<=n-1; ++i)                        /*  main diagonal   */
    {
      k = i - 1;
      a [i] = 2. * (h[k] + h[i])
            + h2[k] * c[k]
            + h2[i] * h3[k] * h3[k]
            + h2[i+1] * c[i];
    }
    a [n] = 2. * (h[n-1] + h[n]) + h2[n-1] * c[n-1]
                         + h2[n] * h3[n-1] * h3[n-1] + h2[1] * c[0];
  }
 /*
     right hand side
 */
  h_var_1 = (fn[1] - fn[0]) * h1[0];
  for (i=1; i<=n-1; ++i, h_var_1=h_var_2)
  {
    h_var_2 = (fn[i+1] - fn[i]) * h1[i];
    rs [i]  = 3. * (h_var_2 - h_var_1);
  }
  rs [n] = 3. * ((fn[1] - fn[0]) * h1[0] - h_var_1);
 /*
     Find coefficients  c[i], i=0, ..., n
     by solving linear system
 */
  if (!rep)
  {
    error = fzyfsy (n, a, b, d, rs, c,            /* LU factorization */
                    &hup[0], &hup[n], &hup[2*n],
                    &hup[3*n], &hup[4*n-2], &hup[5*n-5],
                    &hup[6*n-5], &hup[7*n-5], &hup[8*n-9]);
            /* the 0th entry of the vector denotes the last entry of
               previous vector. This is ok since we do not use it here.
                                                                      */
    if (error != 0)                  return error;
  }
  else
    fzyfsl (n, rs, c, &hup[0], &hup[n]  ,/* for repeated calls no     */
            &hup[2*n], &hup[3*n],        /* factorization !           */
            &hup[4*n-2], &hup[5*n-5],
            &hup[6*n-5], &hup[7*n-5], &hup[8*n-9]);
                                                /* see call of fzyfsy */
  c [0] = c[n];
 /*
     Compute remaining spline coefficients
 */
  a [0] = fn[0] - h2[0] / 3. * h1[0] * (c[1] - c[0])
                + h2[n] / 3. * h1[n-1] * (c[n] - c[n-1]);
  for (i=1; i<=n-1; ++i)
  {
    k = i - 1;
    a [i] = fn[i] - h2[i] / 3. * (c[k] * h1[k] - h3[k] * c[i]
                                             + c[i+1] * h1[i]);
  }
  a [n] = a[0];

  for (i=0; i<=n-1; ++i)
  {
    k = i + 1;
    b [i] = h1[i] * (a[k] - a[i]) - h[i] / 3. * (c[k] + 2. * c[i]);
    d [i] = h1[i] / 3. * (c[k] - c[i]);
  }

  return 0;
}

/* ------------------------------------------------------------------ */
/*.BA*/

int glsptr (int   n,
/*.IX{glsptr}*/
            REAL* x,   REAL* f, REAL* w,
            int   marke_v,
            REAL* px,  REAL* py,
            REAL* a,   REAL* b, REAL* c, REAL* d,
            REAL* phi, REAL* r,
            REAL* phid,
            REAL* help
           )
/***********************************************************************
*                                                                      *
* Compute the coefficients a[i], b[i], c[i] and d[i], i = 0,...,n-1    *
* of a transformed parametric cubic fitting spline for a closed and    *
* smooth curve.                                                        *
.BE*)
* First the program computes from the given nodes (x[i], f[i]) the     *
* transformed polar coordinates (phi[i], r[i]),i=0, ..., n, which      *
* serve as the nodes for computing a non-parametric periodic fiting    *
* spline s(t) according to chapter 11.3 .                              *
*                                                                      *
* With  u := (t-phi[i]),  s(t) has the form:                           *
*                                                                      *
*    s(t) = a[i] + b[i] * u + c[i] * u^2 + d[i] * u^3                  *
*                                                                      *
* where t belongs to the interval [phi[i], phi[i+1]], i=0,...,n-1.     *
*                                                                      *
* Since the nodes phi[i] must be strictly increasing, we will usually  *
* have to compute a translation vector  P = (px,py) and a rotation     *
* of the coordinate system by an angle phid first.                     *
*                                                                      *
* For the strong monotonicity the following are necessary :            *
*   - The point  P must lie in the area describes by (x[i], f[i]) in   *
*     such a way that each ray emanating from P cuts the boundary of   *
*     this region precisely once. P = (px,py) should be supplied by the*
*     user, see  marke_v.                                              *
*   - The points (x[i], f[i]) should be ordered so that their edge     *
*     polygon is traversed counterclockwise from  (x[0], f[0]) to      *
*     (x[n], f[n]). Clearly we must have  (x[n], f[n]) = (x[0], f[0]). *
*                                                                      *
* If these assumptions are met, (phi[i], r[i]) are computed as :       *
*        phi[0] = 0.0                                                  *
*        phi[i] = atan (y'(i)/x'(i)) - phid, i=1, ..., n-1             *
*        phi[n] = 2*PI                                                 *
*        r  [i] = sqrt (X'(i)^2 + Y'(i)^2), i=0, ..., n                *
*        where :  phid  = atan (f[0]/x[0])                             *
*                 X'(i) = x[i] - px,                                   *
*                 Y'(i) = f[i] - py                                    *
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*    Name  Type/size    Meaning                                        *
*   -----------------------------------------------------------------  *
*   int  n              index of last node (n > 5)                     *
*   REAL x [n+1]        nodes : x-values                               *
*   REAL f [n+1]                y-values for  x[i]                     *
*   REAL w [n+1]        weights for  f[i]                              *
*   int  marke_v        Control for origin translation:                *
*                         1 : user supplies P = (px, py)               *
*                         0 : no translation (px = py = 0)             *
*                         2 : the translation vector is computed in    *
*                             glsptr as follows :                      *
*                                px = (xmax+xmin)/2                    *
*                                py = (ymax+ymin)/2                    *
*                                mit: xmax=max(x(i)),                  *
*                                     xmin=min(x(i)),                  *
*                                     ymax=max(f(i)),                  *
*                                     ymin=min(f(i)), i=0,...,n        *
*                             REMARK : This does not ensure that a     *
*                                proper P has been chosen. Check result*
*                                if the return value is  -4            *
*   REAL *px            coordinates of P =(px, py) for  marke_v = 1    *
*   REAL *py                                                           *
*   REAL help [14*n-7]  aux vector                                     *
*                                                                      *
* Output parameters:                                                   *
*                                                                      *
*    Name  Type/size    Meaning                                        *
*   -----------------------------------------------------------------  *
*   REAL a [n+1]        spline coefficients for  s(t)                  *
*   REAL b [n+1]                                                       *
*   REAL c [n+1]                                                       *
*   REAL d [n+1]                                                       *
*   REAL phi [n+1]      (phi,r)[i] are the polar coordinates of the    *
*   REAL r [n+1]        points (x,f)[i] from P starting at angle phid  *
*   REAL *px            P = (px, py)                                   *
*   REAL *py                                                           *
*   REAL *phid          angle of rotation (in radians)                 *
*                                                                      *
* Return value :                                                       *
*                                                                      *
*    0 : all is ok                                                     *
*   -1 : n too small                                                   *
*   -3 : invalid weight w                                              *
*   -4 : not monotonic : phi[i] >= phi[i+1]                            *
*   -5 : x[0] != x[n]  or  f[0] != f[n]  or  w[0] != w[n]              *
*    1 : error in  fzyfsy                                              *
*                                                                      *
* suroutines used :                                                    *
*                                                                      *
*    glsppe : periodic cubic fitting splines                           *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int i, error;
  REAL   xmin, xmax, ymin, ymax, sin_alpha, cos_alpha;

  if (n < 6) return -1;                           /* check input data */
  if ( x[0] != x[n] || f[0] != f[n] || w[0] != w[n] )  return -5;
  for (i=0; i<n+1; i++)
    if (w[i] <= 0.0)  return -3;
  if (marke_v == 0)                                 /* no translation */
  {
    *px = *py = 0.0;
    for (i=0; i<n+1; i++)
    {
      b[i] = x[i];
      c[i] = f[i];
    }
  }
  else                                /* translate origin to  (px,py) */
  {
    if (marke_v == 2)
    {
      xmax = xmin = x[0];
      ymax = ymin = f[0];
      for (i=1; i<n+1; i++)
      {
        xmax = max(xmax,x[i]);
        xmin = min(xmin,x[i]);
        ymax = max(ymax,f[i]);
        ymin = min(ymin,f[i]);
      }
      *px = 0.5*(xmax+xmin);
      *py = 0.5*(ymax+ymin);
    }
    for (i=0; i<n+1; i++)
    {
      b[i] = x[i] - *px;
      c[i] = f[i] - *py;
    }
  }
  for (i=0; i<n+1; i++)                      /* compute r[i]; stop if */
    if ( (r[i] = SQRT(b[i]*b[i]+c[i]*c[i])) == 0.0 )
      return -4;                      /* (px, py) is one of the nodes */

  *phid = ACOS (b[0]/r[0]);            /* Compute rotated coordinates */
  if (c[0] < 0.0)
    *phid = 2*PI - *phid;
  cos_alpha = b[0]/r[0];
  sin_alpha = -c[0]/r[0];
  for (i=0; i<n+1; i++)
  {
    a[i] = b[i]*cos_alpha - c[i]*sin_alpha;
    d[i] = b[i]*sin_alpha + c[i]*cos_alpha;
  }
  phi[0] = 0.0;       /* compute angles phi[i]; stop if not monotonic */
  for (i=1; i<n; i++)
  {
    phi[i] = ACOS(a[i]/r[i]);
    if (d[i] < 0.0)
      phi[i] = 2*PI - phi[i];
    if (phi[i] <= phi[i-1])  return -4;
  }
  phi[n] = 2*PI;
                                       /* Compute spline coefficients */
  error = glsppe (n,phi,r,w,0,a,b,c,d,&help[1],&help[n+2],&help[2*n+3],
                  &help[3*n+4],&help[4*n+4],&help[5*n+4]);
  return error;
}

/* ---------------------------- END glsp.c -------------------------- */
