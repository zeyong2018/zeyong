#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 4.9}
     {The Matrix Inverse via Exchange Steps}
     {The Matrix Inverse via Exchange Steps}*/

/*.BE*/
/* ------------------------- MODULE fpivot.c ------------------------ */

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>

/*.BA*/

int pivot               /* Find the matrix inverse (Exchange steps) ..*/
/*.IX{pivot}*/
          (
           int     n,             /* size of matrix ..................*/
           REAL *  mat[],         /* input matrix ....................*/
           REAL *  inv[],         /* its inverse .....................*/
           REAL *  s,             /* Check sum .......................*/
           REAL *  cond           /* condition number ................*/
          )
/*====================================================================*
 *                                      -1                            *
 *  pivot computes the inverse matrix  A   of a nonsingular matrix A  *
 *  via exchange steps.                                               *
 *  For stability we use both row and column pivot search.            *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *                                                                    *
 *      Find an explicit inverse of a nonsingular n x n matrix A.     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               size of the matrices mat and inv.                    *
 *      mat      REAL   *mat[n];                                      *
 *               Matrix, stored as one vector                         *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      inv      REAL   *inv[n];                                      *
 *               Inverse  matrix for mat.                             *
 *      s        REAL   *s;                                           *
 *               trace of  ( mat * inv ) - n, this must be near zero. *
 *      cond     REAL   *cond;                                        *
 *               conditions number for the max norm;                  *
 *               For  cond = 1.E+k, k digits of accuracy will be lost *
 *               in  inv.                                             *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      Inverse found                                        *
 *      = 1      n < 1 or other wrong input parameter                 *
 *      = 2      memory exceeded                                      *
 *      = 3      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used :                                                 *
 *   ================                                                 *
 *                                                                    *
 *   vmfree(), vmalloc()                                              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Macros in use :                                                  *
 *   ===============                                                  *
 *      SWAP, ABS                                                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :  NULL, MACH_EPS                                 *
 *   ================                                                 *
 *                                                                    *
 *====================================================================*/
{
  int      *permx,                          /*  row permutationen     */
           *permy;                          /*  column permutationen  */

  register k, j, i, ix, iy;                 /*  loop indices          */
  int      nx = 0, ny = 0;                  /*  Pivot indices         */
  REAL     piv,                             /*  aux variables         */
           tmp, norma, normb,
           faktor, h1, h2;
  void     *vmblock;

  if (n < 1) return (1);                    /* input wrong            */

  if (mat == NULL || inv == NULL) return (1);

  for (k = 0; k < n; k++)
    if (mat[k] == NULL || inv[k] == NULL) return (1);

                                             /*  allocate space       */
  vmblock = vminit();
  permx = (int *)vmalloc(vmblock, VVEKTOR, n, sizeof(*permx));
  permy = (int *)vmalloc(vmblock, VVEKTOR, n, sizeof(*permx));
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 2;
  }

  for (i = 0; i < n; i++)
  {
    permx[i] = permy[i] = -1;             /*  initialize permx, permy */
    for (j = 0; j < n; j++)                  /*  store  mat in inv    */
      inv[i][j] = mat[i][j];
  }

  for (i = 0; i < n; i++)
  {
    for (piv = ZERO, ix = 0; ix < n; ix++)   /*  Search for pivot     */
      if (permx[ix] == -1)
      {
        for (iy = 0; iy < n; iy++)
          if (permy[iy] == -1 && ABS(piv) < ABS(inv[ix][iy]))
          {
            piv = inv[ix][iy];          /* store pivot position and   */
            nx = ix;                    /* pivot element              */
            ny = iy;
          }
      }

    if ( ABS(piv) < MACH_EPS )          /* If pivot too small,        */
    {                                   /* mat is nearly singular     */
      vmfree(vmblock);
      return (3);
    }

    permx[nx] = ny;                     /* exchange pvots             */
    permy[ny] = nx;

    tmp = ONE / piv;                    /* Pivot step .....           */
    for (j = 0; j < n; j++)
      if (j != nx)
      {
        faktor = inv[j][ny] * tmp;
        for (k = 0; k < n; k++)          /* ... outside pivot row and */
                                         /*     column                */

          inv[j][k] -= inv[nx][k] * faktor;

        inv[j][ny] = faktor;             /* ... in der Pivotspalte    */
      }

    for (k = 0; k < n; k++)
      inv[nx][k] *= -tmp;                /* ... in pivot row          */
    inv[nx][ny] = tmp;                   /* ... for the pivot         */

  }   /*  end i */

                              /* reverse row and column interchanges  */

  for (i = 0; i < n; i++)
  {                                       /* find j with permx[j] = i */
    for (j = i; j < n; j++)
      if (permx[j] == i) break;

    if (j != i)                           /* row exchange             */
    {                                     /* only for pointers !      */
      SWAP (REAL*, inv[i], inv[j])
      permx[j]   = permx[i];
      permx[i]   = i;
    }
                                        /* Find j with permy[j] = i   */
    for (j = i; j < n; j++)
      if (permy[j] == i) break;

    if (j != i)
    {
      for (k = 0; k < n; k++)           /* column exchange            */
      {
        SWAP (REAL, inv[k][i], inv[k][j])
      }
      permy[j] = permy[i];
      permy[i] = i;
    }
  }  /* end i */

  vmfree(vmblock);                      /*  free memory               */

  *s = norma = normb = ZERO;
  for (i = 0; i < n; i++)
  {                                     /*  form product mat * inv    */
    h1 = h2 = ZERO;                     /*  and its max norm          */
    for (j = 0; j < n; j++)
    {                                  /* *s must nearly be n = tr(I) */
      for (tmp = ZERO, k = 0; k < n; k++)
        tmp += mat[i][k] * inv[k][j];
      *s += ABS (tmp);

      h1 += ABS (mat[i][j]);
      h2 += ABS (inv[i][j]);
    }
    norma = max (h1, norma);           /* Compute max condition number*/
    normb = max (h2, normb);           /* of mat                      */
    *cond = norma * normb;

  } /* end i */

  *s -= (REAL) n;                      /*   *s nearly n ?             */
  return (0);
}

/* --------------------------- END fpivot.c ------------------------- */
