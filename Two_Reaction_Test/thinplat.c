#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE thinplat.c ------------------------ */

/***********************************************************************
*                                                                      *
* Compute two-dimensional surface splines                              *
* ---------------------------------------                              *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 with 80387                         *
* Source:               equivalent QuickBASIC module                   *
* Author:               Elmar Pohl (QuickBASIC)                        *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Date:                 3.15.1993                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  REAL, ONE, ZERO, POW, sqr, LOG, LOG, */
                        /*       SQRT, SetVec                         */
#include <vmblock.h>    /*  for  vminit, vmalloc, VEKTOR, vmcomplete, */
                        /*       vmfree                               */
#include <linpack.h>    /*  for  sspco, sspsl                         */
#include <thinplat.h>   /*  for  prob2, apprx2                        */


/*--------------------------------------------------------------------*/

static int next2     /* Compute two dimensional monomials ............*/
/*.IX{next2}*/
                 (
                  int i,        /* number of existing monomials ......*/
                  int idx[],    /* idx[j]: Power of X in jth monomial */
                  int idy[],    /* idy[j]: Power of Y in jth monomial */
                  int *ixy      /* Monomial index for next monomial ..*/
                 )              /* multiply by X or Y ? ..............*/

/***********************************************************************
* Construct aux arrays for alpha2() in order to determine all two      *
* variable monomial up to degree M.                                    *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* i    Index of most recently computed two variable monomial           *
* idx  [1..i+1] vector of powers of X in the monomials indexed 1 to i  *
* idy  [1..i+1] vector of powers of Y in the monomials indexed 1 to i  *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* idx  same as above; idx[i+1] contains the power of X of the monomial *
*      with index i+1                                                  *
* idy  as above; idy[i+1] contains the power of Y of the (i+1)st       *
*      monomial                                                        *
* ixy  Index of monomials that must be multiplied by  X or Y , to      *
*      obtain the (i+1)st monomial                                     *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* Flag:                                                                *
* TRUE      :  Monom[i+1] = Monom[ixy] * x                             *
* FALSE     :  Monom[i+1] = Monom[ixy] * y                             *
*                                                                      *
***********************************************************************/
{
  int n,
      j;

  n = idx[i] + idy[i];
  if (idx[i] == 0)
  {
    idx[i + 1] = n + 1;
    idy[i + 1] = 0;
    for (j = 1; j <= i; j++)
      if (idx[j] == n)
        break;

    *ixy = j;
    return TRUE;
  }

  else
  {
    idx[i + 1] = idx[i] - 1;
    idy[i + 1] = idy[i] + 1;
    for (j = 1; j <= i; j++)
      if (idx[j] == idx[i] - 1 && idy[j] == idy[i])
        break;

    *ixy = j;
    return FALSE;
  }
}

/*--------------------------------------------------------------------*/

static int alpha2         /* Compute polynomial part of matrix .......*/
/*.IX{alpha2}*/
                 (
                  int  NX,     /* number of nodes ....................*/
                  int  M,      /* derivative order ...................*/
                  REAL x[],    /* nodes : x-coordinates ..............*/
                  REAL y[],    /*         y-coordinates ..............*/
                  REAL a[]     /* Polynomial part P ..................*/
                 )             /* error code .........................*/

/***********************************************************************
* Aux function for  prob2() to compute two-dimensional interpolating   *
* surface splines. The polynomial part of the system matrix for prob2()*
* is computed here.                                                    *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* NX    number of nodes                                                *
* M     given derivative order                                         *
* x  \  [1..NX] node vectors (X(i),Y(i)), i=1,...,NX                   *
* x  /                                                                 *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* a     [1..(nm*(nm+1))/2] vector with condensed matrix that contains  *
*       polynomial part of system matrix; here  nm := NX+(M*(M+1))/2   *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok                                                          *
* = 1: lack of memory                                                  *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, vminit, vmalloc, VEKTOR, vmcomplete, vmfree, ONE, ZERO, next2, *
* SetVec                                                               *
***********************************************************************/
{
  int  *idx,       /* [1..M*(M+1)/2] vector of powers of X            */
       *idy,       /* [1..M*(M+1)/2] vector of powers of Y            */
       i,          /* Number of current monomial                      */
       j,          /* loop counter for evaluating monomials           */
       ixy,        /* Number of the monomial, from which the current  */
                   /* one is formed by multiplying by X or Y          */
       mal_x,      /* control whether to multiply by X or Y           */
       kl,         /* Index of condensed matrix after which the ixy-th*/
                   /* monomial is stored                              */
       kli;        /* ditto for ith monomial                          */
  REAL *xy;        /* Address of vector x or y                        */
  void *vmblock;   /* List of dynamic allocations                     */

  /* ----------- initialize buffers for aux arrays ------------------ */

  vmblock = vminit();
  idx = (int *)vmalloc(vmblock, VVEKTOR, M * (M + 1) / 2, sizeof(*idx));
  idy = (int *)vmalloc(vmblock, VVEKTOR, M * (M + 1) / 2, sizeof(*idy));
  idx--; idy--;                             /* shift index      !!!!! */
  if (! vmcomplete(vmblock))   /* lack of available memory ?          */
  {
    vmfree(vmblock);
    return 1;
  }


  a += (NX * (NX + 1)) / 2;       /* skip kernel function part G in   */
                                  /* condensed matrix                 */

  SetVec(NX, a + 1, ONE);              /* set up first monomial (= 1) */
  idx[1]    = 0;
  idy[1]    = 0;
  a[NX + 1] = ZERO;

  for (kli = NX + 1, i = 2;  /* Compute the monomials  2,...,M(M+1)/2 */
       i <= M * (M + 1) / 2;
       kli += NX + i, i++)
  {

    mal_x = next2(i - 1, idx, idy, &ixy);   /* find index of monomial */
                                            /* that is multiplied by  */
                                            /* x or y                 */
    kl = (ixy + NX + NX) * (ixy - 1) / 2;
    xy = mal_x ? x : y;
    for (j = 1; j <= NX; j++)
       a[kli + j] = a[kl + j] * xy[j];

    SetVec(i, a + kli + NX + 1, ZERO);     /* set rest of matrix zero */
  }

  vmfree(vmblock);
  return 0;
}

/*--------------------------------------------------------------------*/
/*.BA*/

REAL apprx2    /* Compute functional value of a surface spline .......*/
/*.IX{apprx2}*/
           (
            REAL x0,           /* (x0,y0) = place for evaluation .....*/
            REAL y0,
            int  NX,           /* number of spline nodes .............*/
            int  M,            /* derivative order ...................*/
            REAL x[],          /* nodes: x-coordinates ...............*/
            REAL y[],          /*        y-coordinates ...............*/
            REAL c[]           /* Spline coefficients ................*/
           )                   /* error code .........................*/

/***********************************************************************
* Compute functional value for a two-dimensional interpolating         *
* surface spline                                                       *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0,y0   place of evaluation                                          *
* NX      number of nodes                                              *
* M       derivative order                                             *
* x  \    [1..NX] node vectors  (X(i),Y(i)), i=1, ..., NX              *
* x  /                                                                 *
* c       [1..NX+(M*(M+1))/2] spline coefficient vector                *
*         (output of   prob2())                                        *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* value of spline at  (x0,y0)                                          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, POW, sqr, ZERO, LOG                                            *
.BA*)
***********************************************************************/
/*.BE*/
{
  REAL ap,
       r2;
  int  ix,
       iy,
       i;


  /* --- For M=1,2,3 the polynomial parts are available in closed --- */
  /* --- form, thus taking next to no time to use.                --- */
  /* --- For M>3 we represent the monomials in the form           --- */
  /* --- x^ix * y^iy. Their computation is slower and subject to  --- */
  /* --- rounding errors.                                         --- */

  ap = c[NX + 1];

  switch (M)
  {
    case 1:
      break;                                   /* nothing to be done */
    case 2:
      ap += c[NX + 2] * x0 + c[NX + 3] * y0;
      break;
    case 3:
      ap += (c[NX + 2] + c[NX + 4] * x0 + c[NX + 5] * y0) * x0 +
            (c[NX + 3] + c[NX + 6] * y0) * y0;
      break;
    default:
      ix = 0;
      iy = 0;
      for (i = 2; i <= (M * (M + 1)) / 2; i++)
        if (ix == 0)
        {
          ix =  iy + 1;
          iy =  0;
          ap += c[NX + i] * POW(x0, ix);
        }
        else
        {
          ix--;
          iy++;
          if (ix == 0)
            ap += c[NX + i] * POW(y0, iy);
          else
            ap += c[NX + i] * POW(x0, ix) * POW(y0, iy);
        }
  }

  /* ----- Kernel part  E:                                       ---- */
  /* ----- One could use the function  e2() instead, but this    ---- */
  /* ----- would be much slower than the following direct code.  ---- */
  /* ----- die Ausfuehrung aber erheblich verlangsamen wuerde.   ---- */

  for (i = 1; i <= NX; i++)
  {
    r2 = sqr(x[i] - x0) + sqr(y[i] - y0);
    if (r2 != ZERO)
      ap += c[i] * LOG(r2) * POW(r2, M - 1);
  }

  return ap;
}

/*--------------------------------------------------------------------*/

static REAL e2    /* evaluate the two-dimensional kernel function E ..*/
/*.IX{e2}*/
              (
               REAL x,      /* (x,y) = place of evaluation ...........*/
               REAL y,
               int  M       /* derivative order of spline ............*/
              )             /* value of E at  (x,y) ..................*/

/***********************************************************************
* Evaluate the two-dimensional kernel function E at (x,y)              *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x,y   place of evaluation                                            *
* M:    derivative order                                               *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* value of E at  (x,y)                                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, LOG, POW                                                 *
***********************************************************************/
{
  REAL r2;

  r2 = x * x + y * y;

  if (r2 == ZERO)
    return ZERO;
  else
    return LOG(r2) * POW(r2, M - 1);
}

/*--------------------------------------------------------------------*/

/*.BA*/

void ekreistrafo       /* Transformation to unit circle ..............*/
/*.IX{ekreistrafo}*/
                (
                 REAL x[],    /* original or transformed coordinates .*/
                 REAL y[],
                 int  n       /* number of transformed points ........*/
                )

/***********************************************************************
* Aux function for two-dimensional interpolating surface splines.      *
* It transforms a set of points  (X(i),Y(i)) to the unit circle.       *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x  \  [1..n] coordinate vectors of (X(i),Y(i)), that are to be       *
* y  /  transformed to the unit circle, i=1, ..., n                    *
* n     number of points in set                                        *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* x,y   [1..n] vectors with transformed points                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, sqr, SQRT                                                *
.BA*)
***********************************************************************/
/*.BE*/
{
  REAL xq,
       yq,
       r2,
       r2i,
       r;
  int  i;


  for (xq = yq = ZERO, i = 1; i <= n; i++)
    xq += x[i],
    yq += y[i];
  xq /= n;
  yq /= n;

  for (r2 = ZERO, i = 1; i <= n; i++)
  {
    x[i] -= xq;
    y[i] -= yq;
    r2i = sqr(x[i]) + sqr(y[i]);
    if (r2i > r2)
      r2 = r2i;
  }

  for (r = SQRT(r2), i = 1; i <= n; i++)
    x[i] /= r,
    y[i] /= r;
}

/*--------------------------------------------------------------------*/

static int gamma2   /* Compute kernel part of system matrix ..........*/
/*.IX{gamma2}*/
                 (
                  int  NX,       /* number of spline nodes ...........*/
                  int  M,        /* derivative order .................*/
                  REAL x[],      /* nodes: x-coordinates .............*/
                  REAL y[],      /*        y-coordinates .............*/
                  REAL rho,      /* smoothing parameter (>= 0) .......*/
                  REAL w[],      /* weights ..........................*/
                  REAL a[]       /* Kernel function  G ...............*/
                 )               /* error code .......................*/

/***********************************************************************
* Aux function for prob2() to compute two-dimensional interpolating    *
* or smoothing surface splines. We compute the kernel part G here for  *
* the system matrix of the linear system in  prob2().                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* NX    number of nodes                                                *
* x  \  [1..NX] node vectors  (X(i),Y(i)), i=1,...,NX                  *
* x  /                                                                 *
* M     given derivative order                                         *
* rho   smoothing parameter (>= 0).                                    *
*       For  rho = 0  we compute an interpolating spline,              *
*       for  rho > 0  we compute a smoothing spline. The larger rho is *
*       chosen, the closer the spline gets to the fitting plane.       *
* w     for  rho > 0:                                                  *
*       [1..NX] positive weight vector; nodes with large weights are   *
*       given more prominence in the interpolation.                    *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* a     Vector, which was used to copy the kernel part into the system *
*       matrix                                                         *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all is ok                                                       *
* = 1: for rho > 0: one weight  W(i) is <= 0                           *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, e2, ZERO                                                       *
***********************************************************************/
{
  int i,
      k,
      l;


  for (l = i = 1; i <= NX; i++, l++)
  {
    for (k = 1; k < i; k++, l++)
      a[l] = e2(x[k] - x[i], y[k] - y[i], M);

    /* ----------------- main diagonal ------------------------------ */

    if (rho == ZERO)                           /* Interpolation?      */
      a[l] = ZERO;
    else if (w[i] <= ZERO)                     /* negative weight ?   */
      return 1;
    else                                       /* smoothing spline  ? */
      a[l] = rho / w[i];
  }

  return 0;
}

/*--------------------------------------------------------------------*/
/*.BA*/

int prob2    /* compute two-dimensional surface splines ..............*/
/*.IX{prob2}*/
         (
          int  NX,         /* number of spline nodes .................*/
          REAL x[],        /* nodes: x-coordinates ...................*/
          REAL y[],        /*        y-coordinates ...................*/
          REAL z[],        /* values to be smoothed at (x[i],y[i]) ...*/
          int  M,          /* derivative order .......................*/
          REAL rho,        /* smoothing parameter (>= 0) .............*/
          REAL w[],        /* weights ................................*/
          REAL c[]         /* Spline coefficients ....................*/
         )                 /* error code .............................*/

/***********************************************************************
* Compute two-dimensional surface splines or "Thin Plate Splines" for  *
* arbitrary triples  (X(i),Y(i),Z(i)), i=1, ..., NX. Depending on the  *
* parameter rho one can compute interpolating or smoothing surface     *
* splines. The nodes (X(i),Y(i)) must be distinct. The nodes can be    *
* arranged arbitrarily, as "scattered data". We recommend to transform *
* the nodes  (X(i),Y(i)) to the unit circle using the function         *
* ekreistrafo(). The derivative order should be kept low or the        *
* condition number of the linear system will become too large. Tests   *
* suggest to use a derivative order between 3 and 5. Higher orders led *
* to improvements only very rarely. For increasing number of nodes or  *
* less distance between nodes the condition number is also adversely   *
* affected.                                                            *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* NX    number of nodes                                                *
* x  \  [1..NX] node vectors                                           *
* y  /                                                                 *
* z     [1..NX] vector of function values to be interpolated / smoothed*
* M     given derivative order. The spline is computed continuously    *
*       differentiable up to order  M .                                *
* rho   smoothing parameter (>= 0).                                    *
*       for  rho = 0  we compute an interpolating spline,              *
*       for  rho > 0  we compute a smoothing spline. For large rho the *
*       spline approches the fitting plane for the data.               *
* w     for  rho > 0:                                                  *
*       [1..NX] positive weight vector                                 *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* c     [1..NX+(M*(M+1))/2] spline coefficient vector                  *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all is ok                                                       *
* = 1: system matrix is numerically singular.                          *
* = 2: for rho > 0: one weight  W(i) is non-positive.                  *
* = 3: lack of available memory                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, vminit, vmalloc, VEKTOR, vmcomplete, vmfree, alpha2, gamma2,   *
* sspco, ONE, ZERO, sspsl                                              *
.BA*)
***********************************************************************/
/*.BE*/
{
  REAL *a,        /* [1..(nm*(nm+1))/2] vector of symmetric matrix in */
                  /* in column condensed form                         */
       rcond;     /* estimate for reciprocal of matrix condition #    */
  int  nm,        /* size od system matrix                            */
       fehler,    /* error code from  alpha2() and gamma2()           */
       i,         /* Loop counter                                     */
       *pvt;      /* [1..nm] pivot vector for matrix factorization    */
  void *vmblock;  /* List of dynamic allocations                      */


  nm = NX + (M * (M + 1)) / 2;

  /* ----------- create buffers for aux arrays ---------------------- */

  vmblock = vminit();
  a   = (REAL *)vmalloc(vmblock, VEKTOR,  nm * (nm + 1) / 2, 0);
  pvt = (int *) vmalloc(vmblock, VVEKTOR, nm, sizeof(*pvt));
  a--; pvt--;                               /* shift index      !!!!! */
  if (! vmcomplete(vmblock))                /* lack of memory         */
  {
    vmfree(vmblock);
    return 3;
  }

  /* --------------------- form system matrix ----------------------- */

  fehler = alpha2(NX, M, x, y, a);        /* put polynomial part P    */
                                          /* into upper right corner  */
  if (fehler)                                   /* lack of space ?    */
  {
    vmfree(vmblock);
    return 3;
  }

  fehler = gamma2(NX, M, x, y,    /* kernel part goes to upper left   */
                  rho, w, a);     /* corner, both in condensed form   */
  if (fehler)                                   /* negative weight ?  */
  {
    vmfree(vmblock);
    return 2;
  }

  /* ------------------- factor matrix ------------------------------ */

  if (sspco(a, nm, pvt, &rcond))                 /* lack of space ?   */
  {
    vmfree(vmblock);
    return 3;
  }

  if (ONE + rcond == ONE)                        /* Matrix singular ? */
  {
    vmfree(vmblock);
    return 1;
  }

  /* ------------------ solve linear system ------------------------- */

  for (i = 1; i <= NX; i++)           /* initialize right hand side   */
    c[i] = z[i];
  for (i = NX + 1; i <= nm; i++)
    c[i] = ZERO;

  sspsl(a, nm, pvt, c);               /* solve                        */


  vmfree(vmblock);
  return 0;
}

/* ------------------------- END thinplat.c ------------------------- */
