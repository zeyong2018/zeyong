#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE bspline.c ------------------------ */

/***********************************************************************
*                                                                      *
* Approximation via B splines                                          *
* ---------------------------                                          *
*                                                                      *
* exported functions:                                                  *
*   - bspline():   compute points on an open or closed uniform         *
*                  B spline curve                                      *
*   - bspline2():  compute points on an open or closed uniform         *
*                  B spline surface                                    *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 with 80387                         *
* Source:               FORTRAN (bspline), Turbo Pascal (bspline2)     *
* Authors:              Gisela Engeln-Muellges (FORTRAN),              *
*                       Juergen Dietel (Turbo Pascal)                  *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Date:                 12.22.1993                                     *
*                                                                      *
***********************************************************************/

#include <basis.h>         /*  for  REAL, min, ONE                    */
#include <vmblock.h>       /*  for  vminit, vmalloc, VVEKTOR, MATRIX, */
                           /*       vmcomplete, vmfree                */
#include <bspline.h>       /*  for  bspline, bspline2                 */
/*.BA*/



/*.FE{C 12.4.1}{B Spline Curves}{B Spline Curves}*/

/*.BE*/
/*--------------------------------------------------------------------*/

static void deBoor     /* Compute one point on a B spline curve ......*/
/*.IX{deBoor}*/
                  (
                   REAL *d[],   /* given de Boor points ..............*/
                   int  n,      /* number of de Boor points (>=3) ....*/
                   int  k,      /* Order of the B spline (3<=k <=n) ..*/
                   int  m,      /* dimension of space (2,3) ..........*/
                   int  offen,  /* open curve? .......................*/
                   REAL t0,     /* parameter value of the curve point */
                   int  r,      /* index of t0 node interval .........*/
                   REAL *D[],   /* aux array .........................*/
                   REAL x[]     /* desired curve point ...............*/
                  )

/***********************************************************************
* Compute the point corresponding to the parameter value t0 on the     *
* uniform B spline, using the method of de Boor                        *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* d      [0..n-1,0..m-1] matrix with the de Boor points                *
* n      number of de Boor points (n >= 3)                             *
* k      order of the B spline curve (3 <= k <= n)                     *
* m      Dimension of the space that contains the  de Boor points      *
*        (m >= 2)                                                      *
* offen  flag that shows if the de Boor points define an open curve    *
*        (flag set) or a closed curve (flag cleared)                   *
* t0     Parameter value whose point on the B spline curve is desired  *
* r      Index of node interval that contains t0:                      *
*            t(r) <= t0 <= t(r+1)                                      *
*            (open: k-1 <= r <= n-1,  closed: k-1 <= r <= n+k-2)       *
* D:     [0..k-1,0..m-1] aux array. Not locally allocated in order     *
*        to save repeated allocations and releases.                    *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* x      [0..m-1] vector for the point on the B spline curve           *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, min                                                            *
***********************************************************************/

{

#define t(l)     /* return the value of the lth node of an open or  */ \
  ((! offen) ?   /* closed B spline curve                           */ \
    (l)      :   /* n, k and offen are presumed globally defined.   */ \
  (((l) < k) ?   \
    (k - 1)  :   \
  (((l) < n) ?   \
    (l)      :   \
    n            \
  )))

  int  i,        /* row index in de Boor scheme                       */
       j,        /* column index in de Boor scheme                    */
       l,        /* index for coordinates of a point                  */
       rmk1;     /* r - k + 1                                         */
  REAL alpha;    /* factor for computing a column of the de Boor      */
                 /* scheme from the previous one                      */


  rmk1 = r - k + 1;                   /* copy k adjacent points from  */
                                      /* d to D (column 0 of the de   */
  for (i = min(r, n - 1); i >= rmk1;  /* Boor scheme) in order to     */
       i--)                           /* apply the de Boor algorithm  */
    for (l = 0; l < m; l++)           /* to them                      */
      D[i - rmk1][l] = d[i][l];
                                      /* (closed curve, r >= n:       */
  for (i = n; i <= r; i++)            /* take possible additional     */
    for (l = 0; l < m; l++)           /* de Boor points from the      */
      D[i - rmk1][l] = d[i - n][l];   /* beginning of d)              */

  for (j = 0; j < k - 1; j++)         /* compute column j+1 of the de */
    for (i = k - 1; i > j; i--)       /* Boor scheme from column j    */
    {                                 /* and put result back into D   */
      alpha =       (t0           - (REAL)t(i + rmk1)) /
              (REAL)(t(i + r - j) -       t(i + rmk1));
      for (l = 0; l < m; l++)
        D[i][l] = D[i - 1][l] + alpha *
                  (D[i][l] - D[i - 1][l]);
    }

  for (l = 0; l < m; l++)          /* Now the last element of D is    */
    x[l] = D[k - 1][l];            /* the desired point on the curve. */
}



/*--------------------------------------------------------------------*/
/*.BA*/

int bspline     /* Compute points on a B spline curve ................*/
/*.IX{bspline}*/
           (
            REAL *d[],     /* given de Boor points ...................*/
            int  n,        /* number of de Boor points (>=3) .........*/
            int  k,        /* Order of the B spline (3<=k<=n) ........*/
            int  m,        /* dimension of space (2,3) ...............*/
            int  offen,    /* open curve? ............................*/
            REAL *c[],     /* computed points of the curve ...........*/
            int  *nc       /* maximal or actual number of points .....*/
           )               /* error code .............................*/

/***********************************************************************
* Compute at most nc points of an open or closed uniform               *
* B spline curve of order k                                            *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* d      [0..n-1,0..m-1] matrix with the n given de Boor points        *
* n      number of de Boor points  (n >= 3)                            *
* k      Order of the B spline curve (3 <= k <= n)                     *
* m      Dimension of space for de Boor points  (m >= 2)               *
* offen  flag that shows if the de Boor points define an open curve    *
*        (flag set) or a closed curve (flag cleared)                   *
* nc     maximal number of points on curve                             *
*        (nc-1 >= 2*(n-k+1) (open) resp. nc-1 >= 2*n (closed), that is *
*        nc must be chosen large enough that there are at least two    *
*        points in each node interval)                                 *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* c      [0..ncold-1,0..m-1] matrix with the computed points of the    *
*        curve (Here ncold designates the input value of nc.)          *
* nc     number of actually computed points                            *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all is ok                                                       *
* = 1: invalid input:                                                  *
*      n < 3   or   k < 3   or   k > n   or   m < 2   or   nc < 2*ni   *
*      with  ni := n-k+1 (open)  resp.  ni := n (closed)               *
* = 2: lack of available memory                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, vminit, vmalloc, MATRIX, vmcomplete, vmfree, ONE, deBoor       *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  ni,          /* number of node intervals longer than zero      */
       npi,         /* number of points between two adjacent nodes    */
       r,           /* index for node intervals                       */
       i;           /* variable counting points in one node interval  */
  REAL h,           /* Step size                                      */
       t0,          /* current parameter value for finding a point on */
                    /* curve                                          */
       **D;         /* [0..k-1,0..m-1] aux array                      */
  void *vmblock;    /* List of dynamically allocated vectors/matrices */


  /* ----------- Check input data ----------------------------------- */

  ni = n;
  if (offen)
    ni -= k - 1;

  if (n < 3 || k < 3 || k > n || m < 2 ||     /* wrong n, k, m or not */
      *nc - 1 < 2 * ni)                       /* even two points per  */
    return 1;                                 /* interval possible?   */

  /* ----------- allocate aux array --------------------------------- */

  vmblock = vminit();                         /* initialize storage   */
  D       = (REAL **)vmalloc(vmblock, MATRIX,  k, m);
  if (! vmcomplete(vmblock))                  /* lack of memory?      */
  {
    vmfree(vmblock);                          /* free dynamic memory  */
    return 3;                                 /* report error         */
  }

  npi = (*nc - 1) / ni;
  h   = ONE / (REAL)npi;

  deBoor(d, n, k, m, offen,                   /* first point on curve */
         (REAL)(k - 1), k - 1, D, c[0]);
  *nc = 1;

  for (ni += k - 1, r = k - 1; r < ni; r++)   /* rest of points       */
    for (t0 = (REAL)r + h, i = 0; i < npi; i++, t0 += h, (*nc)++)
      deBoor(d, n, k, m, offen, t0, r, D, c[*nc]);


  vmfree(vmblock);                            /* free dynamic memory  */
  return 0;                                   /* zero: report success */
}
/*.BA*/



/*.FE{C 12.4.2}{B Spline Surfaces}{B Spline Surfaces}*/

/*.BE*/
/* ------------------------------------------------------------------ */

int wKurve      /* compute w curve for parameter v = v0 ..............*/
/*.IX{wKurve}*/
          (
           REAL v0,        /* constant parameter of the w curve ......*/
           int  r,         /* index of the node interval with v0 .....*/
           REAL **d[],     /* given de Boor points ...................*/
           int  m,         /* number of de Boor points in v direction */
           int  n,         /* number of de Boor points in w direction */
           int  k,         /* order of the B spline ..................*/
           int  voffen,    /* open v curves? .........................*/
           int  woffen,    /* open w curves? .........................*/
           REAL *c[],      /* computed points on the w curve .........*/
           int  *nw,       /* max./act. number of points in w direct. */
           REAL *dhilf[],  /* [0..n-1] aux vector \  for              */
           REAL *dv[],     /* [0..m-1] aux vector /  de Boor points ..*/
           REAL *D[]       /* [0..m-1,0..2] aux vector for deBoor() ..*/
          )                /* error code .............................*/

/***********************************************************************
* compute maximal nw points on a w curve belonging to parameter v = v0 *
* of an open or closed uniform B spline surface                        *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* v0      constant v parameter of the spline surface that determines   *
*         the w curve to be computed                                   *
* r       Index of v node interval that contains v0:                   *
*             v(r) <= v0 <= v(r+1)                                     *
*             (open: k-1 <= r <= m-1,  closed: k-1 <= r <= m+k-2)      *
* d       [0..m-1,0..n-1,0..2] matrix with m*n given de Boor points    *
* m       number of de Boor points in v direction (m >= 3)             *
* n       number of de Boor points in w direction (n >= 3)             *
* k       order of the B spline surface (3 <= k <= max(m,n))           *
* voffen  flag that shows if the de Boor points in v direction define  *
*         an open curve (flag set) or a closed curve (flag cleared)    *
* woffen  flag that shows if the de Boor points in w direction define  *
*         an open curve (flag set) or a closed curve (flag cleared)    *
* nw      maximal number of points on w curve                          *
*         (nw-1 >= 2*(n-k+1) (open) resp. nw-1 >= 2*n (closed), that   *
*         is nw must be chosen large enough that there are at least    *
*         two points in each w node interval)                          *
* dhilf   [0..n-1] vector for storing the temporary de Boor points     *
*         on which the w curve is based                                *
* dv      [0..m-1] aux vector that stores the addresses of k           *
*         de Boor points from d in v direction, that is k adjacent     *
*         points from a column of the de Boor matrix d (cyclic if      *
*         necessary)                                                   *
* D       [0..m-1,0..2] aux vector for deBoor()                        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* nw      number of actually computed points on w curve                *
* c       [0..nwold-1,0..2] matrix with the computed points of the     *
*         w curve (Here nwold designates the input value of nw.)       *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* Error code. The following values are possible:                       *
* = 0: all is ok                                                       *
* = 1: invalid input in bspline(); possible is:                        *
*      nw-1 < 2*ni                                                     *
*      with  ni := n-k+1 (w-open)  resp.  ni := n (w-closed)           *
* = 2: lack of memory in bspline()                                     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, min, deBoor, bspline                                           *
***********************************************************************/

{
  int  j, l,         /* loop variables                                */
       fehler;       /* error code from bspline()                     */


  /* compute in dhilf the temporary de Boor points on which the       */
  /* w curve belonging to the constant parameter v = v0 is based      */
  /* that lies in the rth node interval                               */

  for (j = 0; j < n; j++)               /* dhilf[j] exactly depends   */
  {                                     /* on the k de Boor points    */
    for (l = r - k + 1;                 /* in v direction             */
         l <= min(r, m - 1);            /* d[r-k+1,j],...,j[r,j]      */
         l++)                           /* (closed curve: take        */
      dv[l] = d[l][j];                  /* possible additional        */
    for (l = 0; l <= r - m; l++)        /* de Boor points (r >= m)    */
      dv[l] = d[l][j];                  /* from the beginning of the  */
                                        /* jth column of d)           */
    deBoor(dv, m, k, 3, voffen,         /* apply the                  */
           v0, r, D, dhilf[j]);         /* de Boor algorithm to dv    */
  }

  fehler = bspline(dhilf, n, k, 3,      /* w curve for the parameter  */
                   woffen, c, nw);      /* v = v0 in the rth interval */


  return fehler;
}



/*--------------------------------------------------------------------*/
/*.BA*/

int bspline2    /* compute a series of w curves on a B spline surface */
/*.IX{bspline2}*/
            (
             REAL **d[],   /* given de Boor points ...................*/
             int  m,       /* number of de Boor points in v direction */
             int  n,       /* number of de Boor points in w direction */
             int  k,       /* order of the B spline ..................*/
             int  voffen,  /* open v curves? .........................*/
             int  woffen,  /* open w curves? .........................*/
             REAL **c[],   /* computed points on the surface .........*/
             int  *nv,     /* max./act. number of points in v direct. */
             int  *nw      /* max./act. number of points in w direct. */
            )              /* error code .............................*/

/***********************************************************************
* compute maximal nv points in v direction times maximal nw points in  *
* w direction on an open or closed uniform B spline surface.           *
* This is done by computing maximal nv w curves with maximal nw points *
* each. The v direction belongs to the first index of the              *
* de Boor matrix d, the w direction to the second index.               *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* d       [0..m-1,0..n-1,0..2] matrix with m*n given de Boor points    *
* m       number of de Boor points in v direction (m >= 3)             *
* n       number of de Boor points in w direction (n >= 3)             *
* k       order of the B spline surface (3 <= k <= max(m,n))           *
* voffen  flag that shows if the de Boor points in v direction define  *
*         an open curve (flag set) or a closed curve (flag cleared)    *
* woffen  flag that shows if the de Boor points in w direction define  *
*         an open curve (flag set) or a closed curve (flag cleared)    *
* nv      maximal number of surface points in v direction              *
*         (nv-1 >= 2*(m-k+1) (open) resp. nv-1 >= 2*m (closed), that   *
*         is nv must be chosen large enough that there are at least    *
*         two points in each v node interval)                          *
* nw      maximal number of points on a w curve                        *
*         (nw-1 >= 2*(n-k+1) (open) resp. nw-1 >= 2*n (closed), that   *
*         is nw must be chosen large enough that there are at least    *
*         two points in each w node interval)                          *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* nv      number of actually computed v curves                         *
* nw      number of actually computed points on a w curve              *
* c       [0..nvold-1,0..nwold-1,0..2] matrix with the computed points *
*         on the B spline surface (Here nvold and nwold designate the  *
*         input values of nv and nw.)                                  *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* Error code. The following values are possible:                       *
* = 0: all is ok                                                       *
* = 1: invalid input parameters:                                       *
*      m < 3   or   n < 3   or   k < 3   or   k > m   or   nv-1 < 2*ni *
*      with  ni := m-k+1 (v-open)  resp.  ni := m (v-closed)           *
* = 2: error in wKurve()                                               *
* = 3: lack of memory                                                  *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, wKurve, vminit, vmalloc, MATRIX, VVEKTOR, vmcomplete, vmfree,  *
* ONE                                                                  *
.BA*)
***********************************************************************/
/*.BE*/

{
  int  ni,          /* number of v node intervals longer than zero    */
       npi,         /* number of v curve points between two adjacent  */
                    /* nodes of the v node vector                     */
       r,           /* index for v node intervals                     */
       i,           /* variable counting w curves in one              */
                    /* v node interval                                */
       fehler;      /* error code from wKurve                         */
  REAL h,           /* step size in one v node interval               */
       v0,          /* current parameter value in v direction for     */
                    /* finding a w curve on the surface               */

                    /* The following three vectors should be defined  */
                    /* in wKurve(). But in order to save repeated     */
                    /* allocations and releases they were moved here. */
       **dhilf,     /* [0..n-1,0..2] vector for for wKurve() for      */
                    /* storing the temporary de Boor points on which  */
                    /* the w curve is based                           */
       **dv,        /* [0..m-1] aux vector for wKurve() that stores   */
                    /* the addresses of k de Boor points from d in    */
                    /* v direction, that is k adjacent points from a  */
                    /* column of the de Boor matrix d (cyclic if      */
                    /* necessary)                                     */
       **D;         /* [0..m-1,0..2] aux vector for deBoor()          */

  void *vmblock;    /* List of dynamically allocated vectors/matrices */


  /* -------------------- check input parameters -------------------- */

  ni = m;
  if (voffen)
    ni -= k - 1;

  if (m < 3 || n < 3 || k < 3 || k > m ||     /* wrong m, n, k or not */
      *nv - 1 < 2 * ni)                       /* even two points per  */
    return 1;                                 /* interval possible?   */

  npi = (*nv - 1) / ni;
  h   = ONE / (REAL)npi;

  /* --------------- allocate dynamic vectors/matrices -------------- */

  vmblock = vminit();                     /* initialize memory block  */
  dhilf   = (REAL **)vmalloc(vmblock, MATRIX,  n, 3);
  dv      = (REAL **)vmalloc(vmblock, VVEKTOR, m, sizeof(*dv));
  D       = (REAL **)vmalloc(vmblock, MATRIX,  m, 3);
  if (! vmcomplete(vmblock))              /* lack of memory?          */
  {
    vmfree(vmblock);                      /* free dynamic memory      */
    return 3;                             /* report error             */
  }

  /* --------------------- compute first w curve -------------------- */

  fehler = wKurve((REAL)(k - 1), k - 1,   /* w curve for parameter    */
                  d, m, n, k, voffen,     /* v = k-1                  */
                  woffen, c[0], nw,       /* in the (k-1)th interval  */
                  dhilf, dv, D);
  if (fehler)
  {
    vmfree(vmblock);                      /* free dynamic memory      */
    return 1 + fehler;                    /* report error             */
  }
  *nv = 1;

  /* -------------------- compute rest of w curves ------------------ */

  for (r = k - 1, ni += k - 1;     /* for all v node intervals whose  */
       r < ni;                     /* length is larger than zero      */
       r++)
    for (i = 1, v0 = (REAL)r + h;         /* for all points in        */
         i <= npi;                        /* the rth interval         */
         i++, (*nv)++, v0 += h)
    {
      fehler = wKurve(v0, r, d, m, n, k,  /* w curve for parameter    */
                      voffen, woffen,     /* v = v0                   */
                      c[*nv], nw,         /* in the rth interval      */
                      dhilf, dv, D);
      if (fehler)
      {
        vmfree(vmblock);                  /* free dynamic memory      */
        return 1 + fehler;                /* report error             */
      }
    }


  vmfree(vmblock);                        /* free dynamic memory      */
  return 0;                               /* zero: report success     */
}

/* -------------------------- END bspline.c ------------------------- */
