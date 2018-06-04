#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"



/*.BA*/
/*.FE{C 12.1}
     {Interpolating Two-Dimensional Cubic Splines}
     {Interpolating Two-Dimensional Cubic Splines for Constructing
      Smooth Surfaces}*/

/*.BE*/
/* ------------------------- MODULE bikub.c ------------------------- */

#include  <basis.h>    /*  for  REAL, max, NULL, ZERO, ONE, TWO,      */
                       /*       THREE, HALF, mat4x4                   */
#include <vmblock.h>
#include  <u_proto.h>  /*  for  trdiag                                */
#include  <bikub.h>    /*  for  bikub1, bikub2, bikub3, xyintv, bsval */

typedef REAL smallmat [4][4];

static int  bik_st1 (int n, int m, mat4x4** A, REAL* x, REAL* h);
static int  bik_st2 (int n, int m, mat4x4** A, REAL* y, REAL* h);
static int  bik_st3 (int n, int m, mat4x4** A, REAL* h);
static int  bik_st4 (int n, int m, mat4x4** A, REAL* h);
static void bik_st5 (int i, smallmat g, REAL* x);
static void bik_st6 (int i, smallmat g, REAL* x);
static void bik_st5to9 (mat4x4** A, REAL* x,REAL* y,int i,int j);
static int  bik_st12 (int n, int m, mat4x4** A, REAL* x, REAL* h);
static int  bik_st22 (int n, int m, mat4x4** A, REAL* y, REAL* h);
static int  bik_st32 (int n, int m, mat4x4** A, REAL* x, REAL* h);
static int  bik_st13 (int n, int m, mat4x4** A, REAL*** fn);
/*.BA*/

int bikub1 (int      m,
/*.IX{bikub1}*/
            int      n,
            mat4x4   ** mat,
            REAL*    x,
            REAL*    y)
/***********************************************************************
* Computes the coefficients of the bicubic spline of algorithm 12.1.   *
.BE*)
*                                                                      *
* Input parameters :                                                   *
*                                                                      *
*    int  m                     number of x intervals                  *
*    int  n                     number of y intervals                  *
*    REAL mat [m+1][n+1][4][4]  array:                                 *
*                               The following entries are used:        *
*                               for i=0,...,m,j=0,...,n:               *
*                                 mat [i][j][0][0] denotes the function*
*                                 value u [i][j]                       *
*                               for j=0 and j=n:   mat [i][j][1][0]    *
*                                 contains the derivatives   p [i][j]  *
*                               for i=0 and i=m:   mat [i][j][0][1]    *
*                                 contains the derivatives   q [i][j]  *
*                               for i=0 and i=m, and for j=0 and j=n:  *
*                                 mat [i][j][1][1] contains the        *
*                                 derivatives    r [i][j]              *
*    REAL x [m+1]               end points of the x intervals          *
*    REAL y [n+1]               end points of the y intervals          *
*                                                                      *
* Output parameter:                                                    *
*                                                                      *
*   REAL mat [m+1][n+1][4][4]  Spline coefficients                     *
*                                (for i=0,...,m-1, j=0,...,n-1)        *
*                                                                      *
* Return value :                                                       *
*                                                                      *
*    0    : no error                                                   *
*    1-2  : error in trdiag in step_1                                  *
*    3    : lack of memory                                             *
*    4    : x[i] not monotonically increasing                          *
*    5-6  : error in trdiag in step_2                                  *
*    7    : lack of memory for aux vectors                             *
*    8    : y[i] not monotonically increasing                          *
*    9-10 : error in trdiag in step_3                                  *
*   11    : lack of memory for aux vectors                             *
*   12-13 : error in  trdiag in step_4                                 *
*   14    : lack of memory                                             *
*                                                                      *
* sub-routines used :  trdiag                                          *
*                      step1, step_2, step_3, step_4,                  *
*                      step_5, step_6, step_5to9 :                     *
*                      These sub-routines perform the sub-steps of     *
*                      Algorithm 12.1.                                 *
.BA*)
***********************************************************************/
/*.BE*/
{
  int   i, j, error;
  REAL *h1, *h2;
  void *vmblock;

  i = max (n, m) + 1;

  vmblock = vminit();
  h1 = (REAL *)vmalloc(vmblock, VEKTOR, i, 0);
  h2 = (REAL *)vmalloc(vmblock, VEKTOR, i, 0);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 9;
  }

  if ((error = bik_st1 (m,n,mat,x,h1)) != 0) {           goto free_all;}
  if ((error = bik_st2 (m,n,mat,y,h2)) != 0) {error+= 4; goto free_all;}
  if ((error = bik_st3 (m,n,mat,  h1)) != 0) {error+= 8; goto free_all;}
  if ((error = bik_st4 (m,n,mat,  h2)) != 0) {error+=11; goto free_all;}

  for (i=0; i<=m-1; i++)
    for (j=0; j<=n-1; j++)
      bik_st5to9 (mat, x, y, i, j);

free_all:
  vmfree(vmblock);
  return (error);
}

/***********************************************************************
*    Procedure of sub-programs  bik_st1, ... , bik_st4 :               *
*                                                                      *
*    1.  Check monotonicity of  x [i] and/or  y [i]                    *
*    2.  determine columns of the tridiagonal system matrix            *
*    3.  compute the right hand sides, correct first and last equation.*
*        Solve system using trdiag; mark for repeated calls later.     *
*    4.  write solution onto  mat                                      *
*    5.  error codes:                                                  *
*        error = 1, 2 : error  in trdiag,                              *
*        error = 3    : lack of memory for aux vectors                 *
*        error = 4    : Monotonicity violated                          *
***********************************************************************/

#define allocate_abcd(j, k)                                           \
  vmblock = vminit();                                                  \
  a = (REAL *)vmalloc(vmblock, VEKTOR, (j), 0);                        \
  b = (REAL *)vmalloc(vmblock, VEKTOR, (j), 0);                        \
  c = (REAL *)vmalloc(vmblock, VEKTOR, (j), 0);                        \
  d = (REAL *)vmalloc(vmblock, VEKTOR, (j), 0);                        \
  if (! vmcomplete(vmblock))                                           \
  {                                                                    \
    vmfree(vmblock);                                                   \
    return (k);                                                        \
  }                                                                    \
  error = 0;

#define free_them                                                      \
  vmfree(vmblock);                                                     \
  return (error);
static int bik_st1 (int m, int n, mat4x4** mat, REAL* x, REAL* h)
/*.IX{bik\unt st1}*/
{
  int i, j, nm1, rep, error;
  REAL*a, *b, *c, *d;
  void *vmblock;

  allocate_abcd (m-1, 3);

  for (nm1=m-1, i=0; i<=nm1; i++)
  {
    h [i] = x [i+1] - x [i];
    if (h [i] <= ZERO)  { error = 4; goto free_all; }
  }
  for (i=0; i<=m-2; i++)
  {
    b [i] = ONE / h [i];
    c [i] = ONE / h [i+1];
    d [i] = TWO * (b [i] + c [i]);
  }
  for (rep=0, j=0; j<=n; j++, rep=1)
  {
    for (i=0; i<=m-2; i++)
      a [i] = THREE * ((mat [i+1][j][0][0] - mat [i][j][0][0]) /
                                                    (h [i] * h [i])
                    + (mat [i+2][j][0][0] - mat [i+1][j][0][0]) /
                                               (h [i+1] * h [i+1]));
    a [ 0 ] -= mat [0][j][1][0] / h [ 0 ];
    a [m-2] -= mat [m][j][1][0] / h [nm1];

    error = trdiag (m-1, b, d, c, a, rep);
    if (error) break;

    for (i=0; i<=m-2; i++)
      mat [i+1][j][1][0] = a [i];
  }
  free_all: free_them;
}

static int bik_st2 (int m, int n, mat4x4** mat, REAL* y, REAL* h)
/*.IX{bik\unt st2}*/
{
  int i, j, nn1, rep, error;
  REAL*a, *b, *c, *d;
  void *vmblock;

  allocate_abcd (n-1, 3);

  for (nn1=n-1,i=0; i<=nn1; i++)
  {
    h [i] = y [i+1] - y [i];
    if (h [i] <= 0.) { error = 4; goto free_all; }
  }
  for (i=0; i<=n-2; i++)
  {
    b [i] = ONE / h [i];
    c [i] = ONE / h [i+1];
    d [i] = TWO * (b [i] + c [i]);
  }
  for (rep=0,i=0; i<=m; i++,rep=1)
  {
    for (j=0; j<=n-2; j++)
       a [j] = THREE * ((mat [i][j+1][0][0] - mat [i][j][0][0]) /
                                                      (h [j] * h [j])
                     + (mat [i][j+2][0][0] - mat [i][j+1][0][0]) /
                                                 (h [j+1] * h [j+1]));
    a [ 0 ] -= (mat [i][0][0][1] / h [ 0 ]);
    a [n-2] -= (mat [i][n][0][1] / h [nn1]);

    error = trdiag (n-1, b, d, c, a, rep);
    if (error) break;

    for (j=0; j<=n-2; j++)
      mat [i][j+1][0][1] = a [j];
  }
  free_all: free_them;
}

static int bik_st3 (int m, int n, mat4x4** mat, REAL* h)
/*.IX{bik\unt st3}*/
{
  int i, j, rep, error;
  REAL*a, *b, *c, *d;
  void *vmblock;

  allocate_abcd (m-1, 3);

  for (i=0; i<=m-2; i++)
  {
    b [i] = ONE / h [i];
    c [i] = ONE / h [i+1];
    d [i] = TWO * (b [i] + c [i]);
  }
  for (rep=0,j=0; j<=n; j+=n,rep=1)
  {
    for (i=0; i<=m-2; i++)
       a [i] = THREE * ((mat [i+1][j][0][1] - mat [i][j][0][1]) /
                                                      (h [i] * h [i])
                     + (mat [i+2][j][0][1] - mat [i+1][j][0][1]) /
                                                 (h [i+1] * h [i+1]));
    a [ 0 ] -= mat [0][j][1][1] / h [ 0 ];
    a [m-2] -= mat [m][j][1][1] / h [m-1];

    error = trdiag (m-1, b, d, c, a, rep);
    if (error) break;

    for (i=0; i<=m-2; i++)
      mat [i+1][j][1][1] = a [i];
  }
  free_them;
}

static int bik_st4 (int m, int n, mat4x4** mat, REAL* h)
/*.IX{bik\unt st4}*/
{
  int i, j, nn2, rep, error;
  REAL*a, *b, *c, *d;
  void *vmblock;

  allocate_abcd (n-1, 3);

  for (nn2=n-2,i=0; i<=nn2; i++)
  {
    b [i] = ONE / h [i];
    c [i] = ONE / h [i+1];
    d [i] = TWO * (b [i] + c [i]);
  }
  for (rep=0,i=0; i<=m; i++,rep=1)
  {
    for (j=0; j<=nn2; j++)
      a [j] = THREE * ((mat [i][j+1][1][0] - mat [i][j][1][0]) /
                                                     (h [j] * h [j])
                    + (mat [i][j+2][1][0] - mat [i][j+1][1][0]) /
                                                (h [j+1] * h [j+1]));
    a [ 0 ] -= mat [i][0][1][1] / h [ 0 ];
    a [n-2] -= mat [i][n][1][1] / h [n-1];

    error = trdiag (n-1, b, d, c, a, rep);
    if (error) break;

    for (j=0; j<=nn2; j++)
      mat [i][j+1][1][1] = a [j];
  }
  free_them;
}

/***********************************************************************
*   The following performs steps 5 to 9 of algorithm 12.1.             *
***********************************************************************/

static
int  bik_ipt[4][4] = { {0,0,0,0}, {0,0,0,0}, { 2, 1,2, 1}, {3,2, 3,2} };
static
REAL bik_dat[4][4] = { {ONE,    ZERO,ZERO, ZERO}, {ZERO,ONE, ZERO,ZERO},
                       {-THREE,-TWO, THREE,-ONE}, {TWO, ONE,-TWO, ONE}
                     };

static void bik_st5 (int i, smallmat g, REAL* x)
/*.IX{bik\unt st5}*/
{
  int  k, l;
  REAL hpt [4], h;

  h = x [i+1] - x [i];
  for (hpt[0]=ONE,k=1; k<=3; k++)
    hpt [k] = hpt [k-1] / h;
  for (k=0; k<=3; k++)
    for (l=0; l<=3; l++)
      g [k][l] = bik_dat [k][l] * hpt [bik_ipt[k][l]];
}

static void bik_st6 (int i, smallmat g, REAL* x)
/*.IX{bik\unt st6}*/
{
  int  k, l;
  REAL hpt [4], h;

  h = x [i+1] - x [i];
  for (hpt[0]=ONE,k=1; k<=3; k++)
    hpt [k] = hpt [k-1] / h;
  for (k=0; k<=3; k++)
    for (l=0; l<=3; l++)
      g [l][k] = bik_dat [k][l] * hpt [bik_ipt[k][l]];
}

static void bik_st5to9 (mat4x4** mat,REAL* x,REAL* y,int i,int j)
/*.IX{bik\unt st5to9}*/
{
  int      k, l, kl;
  smallmat gx, gyt, w, wgyt;

  bik_st5 (i, gx, x);                                    /*  Step  5  */
  bik_st6 (j, gyt, y);                                   /*  Step  6  */

  for (l=0; l<=1; l++)                                   /*  Step  7  */
    for (k=0; k<=1; k++)
    {
      w [ k ][ l ] = mat [ i ][ j ][k][l];
      w [k+2][ l ] = mat [i+1][ j ][k][l];
      w [ k ][l+2] = mat [ i ][j+1][k][l];
      w [k+2][l+2] = mat [i+1][j+1][k][l];
    }

  for (k=0; k<=3; k++)                                   /*  Step  8  */
    for (l=0; l<=3; l++)
    {
      wgyt [k][l] = 0.;
      for (kl=0; kl<=3; kl++)
        wgyt [k][l] += (w [k][kl] * gyt [kl][l]);
    }

  for (k=0; k<=3; k++)
    for (l=0; l<=3; l++)
    {
      w [k][l] = 0.;
      for (kl=0; kl<=3; kl++)
        w [k][l] += (gx [k][kl] * wgyt [kl][l]);
    }

  for (k=0; k<=3; k++)                                   /*  Step  9  */
    for (l=0; l<=3; l++)
      mat [i][j][k][l] = w [k][l];
}

/* ------------------------------------------------------------------ */
/*.BA*/

int bikub2 (int n, int m, mat4x4** mat, REAL* x, REAL* y)
/*.IX{bikub2}*/

/***********************************************************************
*   Compute the coefficients of a bicubic spline surface without       *
*   boundary data for the partial derivatives  (Algorithm 12.2)        *
.BE*)
* ==================================================================== *
*                                                                      *
*   Input parameters:                                                  *
*   -----------------                                                  *
*                                                                      *
*    Name    Type/size                Meaning                          *
*   ------------------------------------------------------------------ *
*    n       int/---                  number of x intervals            *
*    m       int/---                  number of y intervals            *
*    mat     REAL  /[n+1][m+1][4][4]  mat [i][j][0]] contains the      *
*                                     function values for i=0,...,n,   *
*                                     j=0,...,m.                       *
*                                     mat is used as a pointer array   *
*                                     (****mat).                       *
*    x       REAL  /[n+1]             end points of x intervals        *
*    y       REAL  /[m+1]             end points of y intervals        *
*                                                                      *
*                                                                      *
*   Output parameters:                                                 *
*   -----------------                                                  *
*                                                                      *
*    Name    Type/size                Meaning                          *
*   ------------------------------------------------------------------ *
*    mat     REAL  /[n+1][m+1][4][4]  entries  mat[i][j][k][l]         *
*                                     for                              *
*                                     i=0,...,n-1, j=0,...,m-1;        *
*                                     k=0,1,2,3,   l=0,1,2,3           *
*                                                                      *
*                                                                      *
*   Return value :                                                     *
*   --------------                                                     *
*                                                                      *
*    = 0 : no error                                                    *
*    = 1 : lack of available memory                                    *
*    = 2 : }                                                           *
*    = 3 : } Monotonicity error  in x or y                             *
*    = 4 : }                                                           *
*   >= 5 : error  in bikub1                                            *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   subroutines used :        step_12, step_22, step_32, bikub1        *
*   ------------------                                                 *
*                                                                      *
*   Macros used :             max                                      *
*   -------------                                                      *
*                                                                      *
*   Constants used :          NULL                                     *
*   ----------------                                                   *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int error;
  REAL*h;
  void *vmblock;

  vmblock = vminit();
  h = (REAL *)vmalloc(vmblock, VEKTOR, max(n, m), 0);
  if (! vmcomplete(vmblock))
    return 1;

  if ((error = bik_st12(n,m,mat,x,h)) != 0)
  { vmfree(vmblock); return error+1; }

  if ((error = bik_st22(n,m,mat,y,h)) != 0)
  { vmfree(vmblock); return error+2; }

  if ((error = bik_st32(n,m,mat,x,h)) != 0)
  { vmfree(vmblock); return error+3; }

  vmfree (vmblock);

  if ((error = bikub1 (n, m, mat, x, y)) != 0) return error+4;

  return 0;
}
/***********************************************************************
*   The following subroutines realize the steps 1 to 3 of algorithm    *
*   12.2. First we check monotonicity, then we compute the boundary    *
*   values in mat.                                                     *
***********************************************************************/

static int bik_st12 (int n, int m, mat4x4** mat, REAL* x, REAL* h)
/*.IX{bik\unt st12}*/
{
  int k, l, i;

  for (k=0; k<=1; k++)
    for (l=0; l<=n-2; l+=n-2)
    {
      i = k + l;
      h [i] = x [i+1] - x [i];
      if (h [i] <= ZERO) return (1);
      if (n == 2) break;               /* avoid infinite loop for n=2 */
    }
  for (i=0; i<=m; i++)
  {
    mat [0][i][1][0] =
      (mat [1][i][0][0]-mat [0][i][0][0])
                                   * (ONE/h[0] + HALF/(h[0]+h[1]))
    - (mat [2][i][0][0]-mat [1][i][0][0])
                                   * h[0] / (h[1]*TWO*(h[0]+h[1]));
    mat [n][i][1][0] =
      (mat [ n ][i][0][0]-mat [n-1][i][0][0])
                                   * (HALF/(h[n-2]+h[n-1])+ONE/h[n-1])
    - (mat [n-1][i][0][0]-mat [n-2][i][0][0])
                                  * h[n-1]/(TWO*(h[n-2]+h[n-1])*h[n-2]);
  }
  return (0);
}

static int bik_st22 (int n, int m, mat4x4** mat, REAL* y, REAL* h)
/*.IX{bik\unt st22}*/
{
  int k, l, j;

  for (k=0; k<=1; k++)
    for (l=0; l<=m-2; l+=m-2)
    {
      j = k + l;
      h [j] = y [j+1] - y [j];
      if (h [j] <= ZERO) return (1);
      if (m == 2) break;               /* avoid infinite loop for m=2 */
    }
  for (j=0; j<=n; j++)
  {
    mat [j][0][0][1] =
      (mat [j][1][0][0] - mat [j][0][0][0])
                                   * (ONE/h[0]+HALF/(h[0]+h[1]))
    - (mat [j][2][0][0] - mat [j][1][0][0])
                                   * h[0]/(h[1]*TWO*(h[0]+h[1]));
    mat [j][m][0][1] =
      (mat [j][ m ][0][0] - mat [j][m-1][0][0])
                                   * (HALF/(h[m-2]+h[m-1])+ONE/h[m-1])
    - (mat [j][m-1][0][0] - mat [j][m-2][0][0])
                                  * h[m-1]/(TWO*(h[m-2]+h[m-1])*h[m-2]);
  }
  return 0;
}

static int bik_st32 (int n, int m, mat4x4** mat, REAL* x, REAL* h)
/*.IX{bik\unt st32}*/
{
  int k, l, i;

  for (k=0; k<=1; k++)
    for (l=0; l<=n-2; l+=n-2)
    {
      i = k + l;
      h [i] = x [i+1] - x [i];
      if (h [i] <= ZERO) return (1);
      if (n == 2) break;              /* avoid infinite loop for  n=2 */
    }
  for (i=0; i<=m; i+=m)
  {
    mat [0][i][1][1] =
      (mat [1][i][0][1] - mat [0][i][0][1])
                                   * (ONE/h[0]+HALF/(h[0]+h[1]))
     -(mat [2][i][0][1] - mat [1][i][0][1])
                                   * h[0]/(h[1]*TWO*(h[0]+h[1]));
    mat [n][i][1][1] =
      (mat [ n ][i][0][1] - mat [n-1][i][0][1])
                                   * (HALF/(h[n-2]+h[n-1])+ONE/h[n-1])
    - (mat [n-1][i][0][1] - mat [n-2][i][0][1])
                                  * h[n-1]/(TWO*(h[n-2]+h[n-1])*h[n-2]);
  }
  return 0;
}

/* ------------------------------------------------------------------ */
/*.BA*/

int bikub3 (int n, int m, mat4x4** mat,
/*.IX{bikub3}*/
            REAL* x, REAL* y, REAL*** fn)

/***********************************************************************
* Compute the coefficients of bicubic surface splines for given        *
* functional values and surface normals, see Algorithm 12.3.           *
.BE*)
* ==================================================================== *
*                                                                      *
*   Input parameters:                                                  *
*   -----------------                                                  *
*                                                                      *
*    Name    Type/size                Meaning                          *
*   ------------------------------------------------------------------ *
*    n       int/---                  number of x intervals            *
*    m       int/---                  number of y intervals            *
*    mat     REAL  /[n+1][m+1][4][4]  mat [i][j][0][0] contains the    *
*                                     function values for i=0,...,n,   *
*                                     j=0,...,m, as a pointer array    *
*                                     (****mat).                       *
*    x       REAL  /[n+1]             end points of x intervals        *
*    y       REAL  /[m+1]             end points of y intervals        *
*    fn      REAL  /[n+1][m+1][3]     Normal vectors                   *
*                                                                      *
*                                                                      *
*   Output parameters:                                                 *
*   ------------------                                                 *
*                                                                      *
*    Name    Type/size                Meaning                          *
*   ------------------------------------------------------------------ *
*    mat     REAL  /[n+1][m+1][4][4]  we compute all missing           *
*                                     mat [i][j][k][l] for             *
*                                     i=0,...,n-1, j=0,...,m-1,        *
*                                     k=0,1,2.3, l=0,1,2,3             *
*                                                                      *
*                                                                      *
*   Return value :                                                     *
*   --------------                                                     *
*                                                                      *
*   = 0 : no error                                                     *
*   = 1 : third component fn is zero                                   *
*   = 2 : lack of memory for aux arrays                                *
*   = 3 : error in  step_32                                            *
*   = 4 : Monotonicity not correct                                     *
*   = 5 : }  error  in step_3                                          *
*   = 6 : }                                                            *
*   = 7 : }  error  in step_4                                          *
*   = 8 : }                                                            *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   subroutines used :        step_3, step_4, step_32, step_13         *
*   ------------------                                                 *
*                                                                      *
*   Constant used :           NULL                                     *
*   ---------------                                                    *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int error, i, j;
  REAL*h1, *h2;
  void *vmblock;
                                                   /*  Steps  1 and 2 */
  if ((error = bik_st13 (n, m, mat, fn)) != 0) return error;

  vmblock = vminit();
  h1 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  h2 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 2;
  }

  if ((error = bik_st32 (n, m, mat, x, h1)) != 0)         /*  Step  3 */
  {
    vmfree(vmblock);
    return error + 2;
  }

  for (i=0; i<=n-1; i++)
  {
    h1 [i] = x [i+1] - x [i];
    if (h1 [i] <= ZERO)
    {
      vmfree(vmblock);
      return error + 4;
    }
  }

  for (j=0; j<=m-1; j++)
  {
    h2 [j] = y [j+1] - y [j];
    if (h2 [j] <= ZERO)
    {
      vmfree(vmblock);
      return 4;
    }
  }

  if ((error = bik_st3 (n, m, mat, h1)) != 0)
  {
    vmfree(vmblock);
    return error + 4;
  }

  if ((error = bik_st4 (n, m, mat, h2)) != 0)             /*  Step  4 */
  {
    vmfree(vmblock);
    return error + 6;
  }

  vmfree(vmblock);

  for (i=0; i<=n-1; i++)                            /*  Steps  5 to 9 */
    for (j=0; j<=m-1; j++)
      bik_st5to9 (mat, x, y, i, j);

  return 0;
}

/***********************************************************************
* Realize steps 1 and 2 in Algorithm 12.3                              *
***********************************************************************/
static int bik_st13 (int n, int m, mat4x4** mat, REAL*** fn)
/*.IX{bik\unt st13}*/
{
  int i, j;
  for (i=0; i<=n; i++)
    for (j=0; j<=m; j++)
    {
      if (fn [i][j][2] == ZERO) return (1);
      mat [i][j][1][0] = - fn [i][j][0] / fn [i][j][2];
      mat [i][j][0][1] = - fn [i][j][1] / fn [i][j][2];
    }
  return (0);
}

/* ------------------------------------------------------------------ */
/*.BA*/

int bsval (int      m,
/*.IX{bsval}*/
           int      n,
           mat4x4** mat,
           REAL*    x,
           REAL*    y,
           REAL     xcoord,
           REAL     ycoord,
           REAL*    value
          )
/***********************************************************************
* Compute the value of a bicubic spline at a point                     *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*   int  m                      number of x intervals                  *
*   int  n                      number of y intervals                  *
*   REAL mat [m+1][n+1][4][4]   matrix of coefficients, as supplied by *
*                               bikub1 : a pointer array  (****mat)    *
*   REAL x[m+1]                 end points of x intervals              *
*   REAL y[n+1]                 end points of y intervals              *
*   REAL xcoord                 x-coordinate of point                  *
*   REAL ycoord                 y-coordinate of point                  *
*                                                                      *
*                                                                      *
* Output parameter:                                                    *
*                                                                      *
*   REAL value                  value at point                         *
*                                                                      *
* Return value :                                                       *
*                                                                      *
*   = 0 : no error                                                     *
*   = 1 : point lies outside domain of spline                          *
*                                                                      *
*                                                                      *
* subprograms used :  xyintv                                           *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int  i, j, k, l, error;
  REAL xip [4], yjp [4], xi, yj;

  int xyintv (int m, int n, REAL* x, REAL* y, int* i, int* j,
              REAL* xi, REAL* yj, REAL xcoord, REAL ycoord);

  error = xyintv (m,n,x,y,&i,&j,&xi,&yj,xcoord,ycoord);
  if (error) return (error);

  xip [0] = yjp [0] = 1.;
  *value = ZERO;
  for (k=1; k<=3; k++)
  {
    xip [k] = xip [k-1] * xi;
    yjp [k] = yjp [k-1] * yj;
  }
  for (k=0; k<=3; k++)
    for (l=0; l<=3; l++)
      *value += mat [i][j][k][l] * xip [k] * yjp [l];

  return 0;
}
/*.BA*/

int xyintv (int   m,
/*.IX{xyintv}*/
            int   n,
            REAL* x,
            REAL* y,
            int*  i,
            int*  j,
            REAL* xi,
            REAL* yj,
            REAL  xcoord,
            REAL  ycoord
           )

/***********************************************************************
* Computes the intervals in which the point  (xcoord, ycoord) lies.    *
.BE*)
*                                                                      *
* Input parameters:                                                    *
*                                                                      *
*   int   m                 number of  x intervals                     *
*   int   n                 number of  y intervals                     *
*   REAL  x [m+1]           end points of x intervals                  *
*   REAL  y [n+1]           end points of y intervals                  *
*   REAL  xcoord            x-coordinate of point                      *
*   REAL  ycoord            y-coordinate of point                      *
*                                                                      *
*                                                                      *
* Output parameters:                                                   *
*                                                                      *
*   int  i               Number of x interval, containing  xcoord :    *
*                          x[i] <= xcoord <= x[i+1]                    *
*   int  j               Number of y interval, containing  ycoord :    *
*                          y[j] <= ycoord <= y[j+1]                    *
*   REAL xi              relative x-coordinate: xi = xcoord - x[i]     *
*   REAL yj              relative y-coordinate: yj = ycoord - y[j]     *
*                                                                      *
* Return value :                                                       *
*                                                                      *
*   = 0 : no error                                                     *
*   = 1 : point lies outside domain of spline                          *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{

  int up, low, mid;

  if (xcoord < x [0] || xcoord > x [m]) return 1;
  if (ycoord < y [0] || ycoord > y [n]) return 1;

  low = 0;
  up  = m;
  do
  {
    mid = (up + low) >> 1;
    if      (xcoord < x [ mid ])  up  = mid;
    else if (xcoord > x [mid+1])  low = mid;
  }
  while (xcoord < x [mid] || xcoord > x [mid+1]);

  *i  = mid;
  *xi = xcoord - x [mid];

  low = 0;
  up  = n;
  do
  {
    mid = (up + low) >> 1;
    if      (ycoord < y [ mid ])  up  = mid;
    else if (ycoord > y [mid+1])  low = mid;
  }
  while (ycoord < y [mid] || ycoord > y [mid+1]);

  *j  = mid;
  *yj = ycoord - y [mid];

  return 0;
}

/* --------------------------- END bikub.c -------------------------- */
