#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE linpack.c ------------------------ */

/***********************************************************************
*                                                                      *
* Solution of condensed symmetric linear systems of equations          *
* -----------------------------------------------------------          *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 with 80387                         *
* Source:               equivalent QuickBASIC module, FORTRAN 77 module*
* Authors:              Elmar Pohl (QuickBASIC),                       *
*                       Michael Groenheim, Ina Hinze (FORTRAN 77)      *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Literature:           Dongarra/Moler/Bunch/Stewart: LINPACK User's   *
*                       Guide, SIAM, 7th ed. 1989, [DONG88]            *
* Date:                 3.15.1993                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  REAL, FABS, ZERO, SWAP, abs, ONE,    */
                        /*       SQRT, EIGHT, FALSE, TRUE, sign,      */
                        /*       skalprod                             */
#include <vmblock.h>    /*  for  vminit, vmalloc, VEKTOR, vmcomplete, */
                        /*       vmfree                               */
#include <linpack.h>    /*  for  sspfa, sspsl, sspco                  */



/* ------------------------------------------------------------------ */

static void init0_vector     /* initialize vector to be zero .........*/
/*.IX{init0\unt vector}*/
                        (
                         REAL v[],             /* Vector .............*/
                         int  n                /* size of vector .....*/
                        )

/***********************************************************************
* initialize vector v with zeros                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO                                                           *
***********************************************************************/

{
  while (n-- != 0)
    *v++ = ZERO;
}



/* ------------------------------------------------------------------ */

static void incabs_vector/* add absolute value vector of w to v ......*/
/*.IX{incabs\unt vector}*/
                         (
                          REAL v[],      /* vector ...................*/
                          REAL w[],      /* second vector ............*/
                          int  n         /* Vector size ..............*/
                         )

/***********************************************************************
* add the modulus of each component of w to each corresponding         *
* component of v.                                                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  while (n-- != 0)
    *v++ += FABS(*w++);
}



/* ------------------------------------------------------------------ */

static void mul_vector    /* scalar multiply a vector ................*/
/*.IX{mul\unt vector}*/
                      (
                       REAL v[],           /* Vektor .................*/
                       int  n,             /* Vector size ............*/
                       REAL faktor         /* scaling factors ........*/
                      )

/***********************************************************************
* multiply each element v[0],v[1],...,v[n-1] of v by  faktor           *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  while (n-- != 0)
    *v++ *= faktor;
}



/* ------------------------------------------------------------------ */

static void incmul_vector  /* add a multiple of one vector to another */
/*.IX{incmul\unt vector}*/
                         (
                          int  n,         /* Vector size .............*/
                          REAL faktor,    /* factor to be used .......*/
                          REAL quelle[],  /* vector to be added ......*/
                          REAL ziel[]     /* resulting vector ........*/
                         )

/***********************************************************************
* add the elements of the vector  quelle, each multiplied by faktor,   *
* to the vector  ziel                                                  *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  while (n-- != 0)
    *ziel++ += faktor * *quelle++;
}



/* ------------------------------------------------------------------ */

static int indmax    /* find index of largest vector element in magn. */
/*.IX{indmax}*/
                 (
                  REAL v[],             /* Vector ....................*/
                  int  n                /* Vector size ...............*/
                 )                      /* Index + 1 .................*/

/***********************************************************************
* Find the index (+ 1) of the largest entry in a vector v in magnitude *
* Return values are  1,2,...,n  or  0 if the vector is empty  (n < 1). *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, FABS                                                           *
***********************************************************************/

{
  REAL maxwert;             /* largest element in magnitude           */
  int  maxind,              /* Index  + 1                             */
       i;                   /* loop counter                           */

  if (n < 1)
    maxind = 0;
  else
    for (maxwert = FABS(*v++), maxind = 1, i = 2; i <= n; i++, v++)
      if (FABS(*v) > maxwert)
        maxwert = FABS(*v),
        maxind  = i;

  return maxind;
}



/* ------------------------------------------------------------------ */

static void swap_vector    /* swap two vectors .......................*/
/*.IX{swap\unt vector}*/
                       (
                        int  n,                /* Vector size ........*/
                        REAL v[],              /* 1st Vector .........*/
                        REAL w[]               /* 2nd Vector .........*/
                       )

/***********************************************************************
* exchange the elements of v with those of w and vice versa            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, SWAP                                                           *
***********************************************************************/

{
  for ( ; n-- != 0; v++, w++)
    SWAP(REAL, *v, *w);
}



/* ------------------------------------------------------------------ */

static REAL einsnorm     /* compute 1-norm of a REAL vector ..........*/
/*.IX{einsnorm}*/
                    (
                     REAL v[],                  /* Vector ............*/
                     int  n                     /* Vector size .......*/
                    )                           /* 1-Norm ............*/

/***********************************************************************
* Find the sum of the magnitudes of the components of v                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, FABS                                                     *
***********************************************************************/

{
  REAL norm;     /* intermediate sum, finallly becoming the norm of v */

  for (norm = ZERO; n-- != 0; v++)
    norm += FABS(*v);

  return norm;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int sspfa     /* decompose a condensed symmetric matrix ..............*/
/*.IX{sspfa}*/
          (
           REAL ap[],     /* upper triangle of matrix, condensed .....*/
           int  n,        /* size of matrix ..........................*/
           int  pvt[]     /* Pivot vector ............................*/
          )               /* singular pivot blocks ? .................*/

/***********************************************************************
* Factor a real symmetric matrix A given in condensed form using       *
* elimination with symmetric pivoting.                                 *
* This function is usually called from sspco(). It can be used indepen-*
* dently if a condition estimate is not required.                      *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* ap    [1..n*(n+1)/2] vector, the upper triangle of A in condensed    *
*       form (for condensing, see  sspco())                            *
* n     size of matrix A                                               *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ap    a diagonal block matrix D of 1 x 1 or 2 x 2 blocks in condensed*
*       form and multiliers. Tis represents the factorization :        *
*       A  =  U * D * U',  with U a product of permutation matrices    *
*       and upper triangular matrices.                                 *
* pvt   [1..n] vector of pivot indices                                 *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: factorization is ok                                             *
* = k: kth pivot block is singular.                                    *
*      This is not fatal for this function, but on a subsequent call   *
*      of  sspsl(), division by zero is likely.                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ONE, SQRT, EIGHT, ZERO, FABS, indmax, FALSE, TRUE,             *
* swap_vector, SWAP, incmul_vector                                     *
.BA*)
***********************************************************************/
/*.BE*/

{
#define ALPHA  (ONE + SQRT((REAL)17.0)) / EIGHT /* to determine pivot */
                                                /* block size         */
  REAL absakk,  /* modulus of  A[k][k]                                */
       rowmax,  /* largest non zero non diagonal entry in row imax    */
       colmax,  /* largest non zero non diagonal entry in column k    */
       ak,      /* A[k][k] / A[k-1][k]                                */
       akm1,    /* A[k-1][k-1] / A[k-1][k]                            */
       bk,
       bkm1,
       nenner,
       amulk,
       amulkm1;
  int  info,    /* Return value: Index of a singular pivot block      */
       swap,    /* indicates whether swaps were needed                */
       k,       /* column counter of A                                */
       ik,      /* leading index (- 1) of the kth column of A in ap   */
       kstep,   /* size of a pivot block                              */
       ikm1,    /* leading index (- 1) of the (k-1)st column of A     */
       kk,      /* Index of diagonal entry  A[k][k] in ap             */
       km1k,    /* Index of  A[k-1][k] in ap                          */
       imax,    /* Index of largest size element in kth column of A   */
       im = 0,  /* leading index (- 1) of the imaxth column of A in ap*/
       imj,
       j,
       jmax,
       jmim,
       imim,
       ij,
       jj,
       jk,
       jkm1;


  info = 0;

  ik = (n * (n - 1)) / 2;

  for (k = n; k != 0; k -= kstep)
  {
    if (k <= 1)
    {
      pvt[1] = 1;
      if (ap[1] == ZERO)
        info = 1;
      break;
    }

    /* ----- The following statements check which elimination   ----- */
    /* ----- to use. Afterwards kstep contains the size of the  ----- */
    /* ----- pivot block, and swap indicates whether swaps are  ----- */
    /* ----- used.                                              ----- */

    kk     = ik + k;
    absakk = FABS(ap[kk]);

    imax   = indmax(ap + ik + 1, k - 1);       /* find largest non    */
    colmax = FABS(ap[ik + imax]);              /* diagonal entry in   */
                                               /* column  k           */
    if (absakk >= ALPHA * colmax)
      kstep = 1,
      swap  = FALSE;
    else
    {
      rowmax = ZERO;                          /* find largest non     */
      im     = (imax * (imax - 1)) / 2;       /* diagonal element in  */
      imj    = im + 2 * imax;                 /* row imax             */
      for (j = imax + 1; j <= k; imj += j, j++)
        if (FABS(ap[imj]) > rowmax)
          rowmax = FABS(ap[imj]);
      if (imax != 1)
      {
        jmax = indmax(ap + im + 1, imax - 1);
        jmim = jmax + im;
        if (FABS(ap[jmim]) > rowmax)
          rowmax = FABS(ap[jmim]);
      }
      imim = imax + im;
      if (FABS(ap[imim]) >= ALPHA * rowmax)
        kstep = 1,
        swap  = TRUE;
      else if (absakk >= ALPHA * colmax * (colmax / rowmax))
        kstep = 1,
        swap  = FALSE;
      else
        kstep = 2,
        swap  = (imax != k - 1);
    }

    if (absakk == ZERO && colmax == ZERO)  /* nur Nullen in Spalte k? */

      info   = k,        /* report a singular pivot block             */
      pvt[k] = k;        /* start next loop                           */

    else if (kstep == 1)                          /* 1x1 pivot block? */
    {
      if (swap)
      {
        swap_vector(imax, ap + im + 1, ap + ik + 1);
        imj = ik + imax;
        for (jj = imax; jj <= k; jj++)
        {
          j  = k + imax - jj;
          jk = ik + j;
          SWAP(REAL, ap[jk], ap[imj]);
          imj -= j - 1;
        }
      }

      ij = ik - (k - 1);                               /* Elimination */
      for (jj = 1; jj < k; ij -= j - 1, jj++)
      {
        j      = k - jj;
        jk     = ik + j;
        amulk  = -ap[jk] / ap[kk];
        incmul_vector(j, amulk, ap + ik + 1, ap + ij + 1);
        ap[jk] = amulk;
      }

      pvt[k] = swap ? imax : k;                  /* set pivot index   */
    }

    else                                          /* 2x2 pivot block? */
    {
      km1k = ik + k - 1;
      ikm1 = ik - (k - 1);

      if (swap)
      {
        swap_vector(imax, ap + im + 1, ap + ikm1 + 1);
        imj = ikm1 + imax;
        for (jj = imax; jj < k; jj++, imj -= j - 1)
        {
          j    = k - 1 + imax - jj;
          jkm1 = ikm1 + j;
          SWAP(REAL, ap[jkm1], ap[imj]);
        }
        SWAP(REAL, ap[km1k], ap[ik + imax]);
      }

      if (k != 2)                                      /* Elimination */
      {
        ak      = ap[kk] / ap[km1k];
        akm1    = ap[ik] / ap[km1k];
        nenner  = ONE - ak * akm1;
        ij      = ik - (k - 1) - (k - 2);
        for (jj = 1; jj < k - 1; jj++, ij -= j - 1)
        {
          j       = k - 1 - jj;
          jk      = ik + j;
          bk      = ap[jk] / ap[km1k];
          jkm1    = ikm1 + j;
          bkm1    = ap[jkm1] / ap[km1k];
          amulk   = (akm1 * bk - bkm1) / nenner;
          amulkm1 = (ak * bkm1 - bk)   / nenner;
          incmul_vector(j, amulk, ap + ik + 1, ap + ij + 1);
          incmul_vector(j, amulkm1, ap + ikm1 + 1, ap + ij + 1);
          ap[jk]    = amulk;
          ap[jkm1]  = amulkm1;
        }
      }

      pvt[k]     = swap ? (-imax) : (1 - k);   /* set pivot indices   */
      pvt[k - 1] = pvt[k];
    }

    ik -= k - 1;
    if (kstep == 2)
      ik -= k - 2;
  }


  return info;
}



/* ------------------------------------------------------------------ */
/*.BA*/

void sspsl    /* Solve linear system for a symmetric condensed matrix */
/*.IX{sspsl}*/
          (
           REAL ap[],      /* Vector with condensed factorization ....*/
           int  n,         /* size of matrix .........................*/
           int  pvt[],     /* Pivot indices ..........................*/
           REAL b[]        /* right hand side/solution vector ........*/
          )

/***********************************************************************
* Solve the linear system  A * x = b  for a real symmetric matrix A in *
* condensed form. This function needs to know the factorization of A   *
* from  sspco() or sspfa().                                            *
.BE*)
*                                                                      *
* Eingabeparameter:                                                    *
* =================                                                    *
* ap    [1..n*(n+1)/2] vector with the factorization                   *
*       (Output from  sspco() or sspfa())                              *
* n     size of matrix                                                 *
* pvt   [1..n] pivot vector from  sspco() or sspfa()                   *
* b     [1..n] right hand side vector                                  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* b     solution vector                                                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, SWAP, incmul_vector, abs, skalprod                             *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL ak,          /* A[k][k] / A[k-1][k]                            */
       akm1,        /* A[k-1][k-1] / A[k-1][k]                        */
       akm1k,       /* A[k-1][k]                                      */
       bk,
       bkm1,
       nenner;
  int  k,           /* column counter of A                            */
       ik,          /* leading index (- 1) for kth column of A in ap  */
       ikm1,        /* leading index (- 1) of (k-1)st column of A     */
       kk,          /* Index of diagonal element A[k][k] in ap        */
       pvtind,      /* current pivot index                            */
       ikp1;


  ik = (n * (n - 1)) / 2;
  for (k = n; k != 0; )
  {
    kk = ik + k;
    if (pvt[k] >= 0)                              /* 1x1 pivot block? */
    {
      if (k != 1)
      {
        pvtind = pvt[k];
        if (pvtind != k)
          SWAP(REAL, b[k], b[pvtind]);
        incmul_vector(k - 1, b[k], ap + ik + 1, b + 1);
      }
      b[k] /= ap[kk];
      k--;
      ik -= k;
    }

    else                                          /* 2x2 pivot block? */
    {
      ikm1 = ik - (k - 1);
      if (k != 2)
      {
        pvtind = abs(pvt[k]);
        if (pvtind != k - 1)
          SWAP(REAL, b[k - 1], b[pvtind]);
        incmul_vector(k - 2, b[k],     ap + ik + 1,   b + 1);
        incmul_vector(k - 2, b[k - 1], ap + ikm1 + 1, b + 1);
      }
      akm1k    = ap[ik + k - 1];
      ak       = ap[kk]   / akm1k;
      akm1     = ap[ik]   / akm1k;
      bk       = b[k]     / akm1k;
      bkm1     = b[k - 1] / akm1k;
      nenner   =  ak * akm1 - ONE;
      b[k]     =  (akm1 * bk - bkm1) / nenner;
      b[k - 1] =  (ak * bkm1 - bk)   / nenner;
      k        -= 2;
      ik       -= (k + 1) + k;
    }
  }

  for (ik = 0, k = 1; k <= n; )
    if (pvt[k] >= 0)                              /* 1x1 pivot block? */
    {
      if (k != 1)
      {
        b[k]   += skalprod(ap + ik + 1, b + 1, k - 1);
        pvtind =  pvt[k];
        if (pvtind != k)
          SWAP(REAL, b[k], b[pvtind]);
      }
      ik += k;
      k++;
    }

    else                                          /* 2x2 pivot block? */
    {
      if (k != 1)
      {
        b[k]     += skalprod(ap + ik + 1, b + 1, k - 1);
        ikp1     =  ik + k;
        b[k + 1] += skalprod(ap + ikp1 + 1, b + 1, k - 1);
        pvtind   =  abs(pvt[k]);
        if (pvtind != k)
          SWAP(REAL, b[k], b[pvtind]);
      }
      ik += k + k + 1;
      k  += 2;
    }
}



/* ------------------------------------------------------------------ */

static void ux_b_loesen  /* solve lin. system w. lower triang. matrix */
/*.IX{ux\unt b\unt loesen}*/
                       (
                        REAL ap[],    /* condensed matrix ............*/
                        int  n,       /* size of matrix ..............*/
                        int  pvt[],   /* Pivot indices ...............*/
                        REAL b[]      /* right hand side/solution ....*/
                       )

/***********************************************************************
* Solve the linear system  U * X = B  for U unit lower triangular.     *
* Such systems occur in sspco() as  U' * Y = W  and  U' * Z = V .      *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* ap    [1..n*(n+1)/2] vector of condense and factored matrix          *
*       A  =  U * D * U'  (Output from sspfa())                        *
* n     size of matrix  A                                              *
* pvt   [1..n] pivot vector from  spfa()                               *
* b     [1..n] right hand side                                         *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* b     solution                                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, skalprod, SWAP                                                 *
***********************************************************************/

{
  int k,             /* column counter for matrix                     */
      ik,            /* leading index (- 1) of kth column of A in ap  */
      kstep,         /* step size when moving through columns of A:   */
                     /* 1 or 2, i.e. size of pivot block              */
      pvtind;        /* Pivot index |pvt[k]|                          */


  for (ik = 0, k = 1; k <= n; k += kstep)
  {
    kstep = (pvt[k] >= 0) ? 1 : 2;
    if (k != 1)
    {
      b[k] += skalprod(ap + ik + 1, b + 1, k - 1);
      if (kstep == 2)
        b[k + 1] += skalprod(ap + ik + k + 1, b + 1, k - 1);
      if ((pvtind = abs(pvt[k])) != k)
        SWAP(REAL, b[k], b[pvtind]);
    }
    ik += k;
    if (kstep == 2)
      ik += k + 1;
  }
}



/* ------------------------------------------------------------------ */
/*.BA*/

int sspco    /* factor condensed symmetric matrix, estimate condition */
/*.IX{sspco}*/
         (
          REAL ap[],     /* upper triangle of matrix, condensed ......*/
          int  n,        /* size of matrix ...........................*/
          int  pvt[],    /* Pivot indices ............................*/
          REAL *rcond    /* estimate for reciprocal of condition # ...*/
         )               /* error code ...............................*/

/***********************************************************************
* Factor a real symmetric matrix A in condensed form and estmate its   *
* condition number.                                                    *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* ap    [1..n*(n+1)/2] vector with the upper triangle of A in condensed*
*       form. In condensing, the columns of the upper triangle are     *
*       stored in string form in the one-dimensional vector ap.        *
*       The following code would achieve this condensation of A in ap :*
*                   for (j = 1, k = 0; j <= n; j++)                    *
*                     for (i = 1; i <= j; i++)                         *
*                       ap[++k] = A[i][j];                             *
* n     size of matrix                                                 *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ap    contains the block diagonal D with blocks of size 2 or 1 and   *
*       corresponding multipliers. This corresponds to the             *
*       factorization   A = U * D * U',  where U is the product of     *
*       permutation and upper triangular matrices.                     *
* pvt   [1..n] pivot vector                                            *
* rcond estimate for the reciprocal of the matrix condition number:    *
*       Perturbations of the system of order Epsilon will result in    *
*       perturbations of the solution of order Epsilon / rcond.        *
*       If rcond satisfies   1 + rcond  =  1, the system matrix is     *
*       judged to be numerically singular.                             *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all is ok                                                       *
* = 1: lack of available memory                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, vminit, vmalloc, VEKTOR, vmcomplete, vmfree, sumabs, FABS,     *
* ZERO, sspfa, ONE, init0_vector, SWAP, sign, incmul_vector,           *
* mul_vector, abs, incabs_vector, einsnorm                             *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL *z,          /* [1..n] vector for various right hand sides     */
                    /* used in condition estimation                   */
       anorm,       /* column sum norm of  A                          */
       ek,          /* current right hand side for linear system      */
                    /* U * D * W = E                                  */
       ynorm,       /* Norm(Y) / Norm(Z)                              */
       s,           /* Scaling factor for vectors                     */
       ak,          /* A[k][k] / A[k-1][k]                            */
       akm1,        /* A[k-1][k-1] / A[k-1][k]                        */
       akm1k,       /* A[k-1][k]                                      */
       bk,
       bkm1,
       nenner;
  int  k,           /* column counter of A in norm/condition work     */
       ik,          /* leading index (- 1) of kth column of A in ap   */
       kstep,       /* Step sizes for covering A; 1 or 2, depending   */
                    /* on pivot block sizes                           */
       pvtind,      /* current pivot index  |pvt[k]|                  */
       ikm1,        /* leading index (- 1) of (k-1)th column of A     */
       kk,          /* Index for diagonal element A[k][k] in ap       */
       k1,          /* leading index for kth column of A in ap        */
       kps;         /* Index of z entry for swapping with z[pvtind]   */
                    /* if needed                                      */
  void *vmblock;    /* List of dynamic allocations                    */


  /* ----------- allocate aux storage ------------------------------- */

  vmblock = vminit();                 /* initialize storage buffers   */
  z = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  z--;                                      /* shift index  !!!!!     */
  if (! vmcomplete(vmblock))                /* lack of storage ?      */
  {
    vmfree(vmblock);
    return 1;
  }


  for (k1 = 1, k = 1; k <= n; k1 += k, k++)  /* find column sum norm  */
  {                                          /* of A                  */
    z[k] = einsnorm(ap + k1, k);             /* z[k]: kth column sum  */
    if (k >= 2)
      incabs_vector(z + 1, ap + k1, n);
  }
  for (anorm = ZERO, k = 1; k <= n; k++)              /* find max     */
    if (z[k] > anorm)                                 /* column sum   */
      anorm = z[k];


  sspfa(ap, n, pvt);                              /* decompose matrix */


  /* ---------------- estimate condition of A ----------------------- */

  /* rcond  =  1 / (Norm(A) * norm estimate of A^{-1} ).              */
  /* The norm of the inverse is estimated by   Norm(Z) / Norm(Y),     */
  /* where  A * Z = Y  and  A * Y = E.  The components of E are       */
  /* selected so that the entries of W in  U * D * W = E  grow        */
  /* maximally.                                                       */


  /* --------------- first solve   U * D * W = E  ------------------- */

  init0_vector(z + 1, n);              /* initialize z = zero vector  */
  ek = ONE;
  ik = (n * (n - 1)) / 2;

  for (k = n; k != 0; )
  {
    kk     = ik + k;
    ikm1   = ik - (k - 1);
    kstep  = (pvt[k] >= 0) ? 1 : 2;
    pvtind = abs(pvt[k]);
    kps    = k + 1 - kstep;
    if (pvtind != kps)
      SWAP(REAL, z[kps], z[pvtind]);
    if (z[k] != ZERO)
      ek = sign(ek, z[k]);
    z[k] += ek;
    incmul_vector(k - kstep, z[k], ap + ik + 1, z + 1);

    if (kstep == 2)                               /* 2x2 pivot block? */
    {
      if (z[k - 1] != ZERO)
        ek = sign(ek, z[k - 1]);
      z[k - 1] += ek;
      incmul_vector(k - kstep, z[k - 1], ap + ikm1 + 1, z + 1);
      akm1k    = ap[ik + k - 1];
      ak       = ap[kk]   / akm1k;
      akm1     = ap[ik]   / akm1k;
      bk       = z[k]     / akm1k;
      bkm1     = z[k - 1] / akm1k;
      nenner   = ak * akm1 - ONE;
      z[k]     = (akm1 * bk - bkm1) / nenner;
      z[k - 1] = (ak * bkm1 - bk)   / nenner;
    }

    else                                          /* 1x1 pivot block? */
    {
      if (FABS(z[k]) > FABS(ap[kk]))
      {
        s = FABS(ap[kk]) / FABS(z[k]);
        mul_vector(z + 1, n, s);
        ek *= s;
      }
      if (ap[kk] == ZERO)
        z[k] = ONE;
      else
        z[k] /= ap[kk];
    }

    k  -= kstep;
    ik -= k;
    if (kstep == 2)
      ik -= k + 1;
  }

  s = 1 / einsnorm(z + 1, n);                          /* normalize W */
  mul_vector(z + 1, n, s);


  /* --------------------- Solve  U' * Y = W  ----------------------- */

  ux_b_loesen(ap, n, pvt, z);

  s = ONE / einsnorm(z + 1, n);                        /* normalize Y */
  mul_vector(z + 1, n, s);


  /* --------------------- Solve   U * D * V = Y  ------------------- */

  ynorm = ONE;
  ik    = (n * (n - 1)) / 2;

  for (k = n; k != 0; )
  {
    kk    = ik + k;
    ikm1  = ik - (k - 1);
    kstep = (pvt[k] >= 0) ? 1 : 2;
    if (k != kstep)
    {
      pvtind  = abs(pvt[k]);
      kps     = k + 1 - kstep;
      if (pvtind != kps)
        SWAP(REAL, z[kps], z[pvtind]);
      incmul_vector(k - kstep, z[k], ap + ik + 1, z + 1);
      if (kstep == 2)
        incmul_vector(k - kstep, z[k - 1], ap + ikm1 + 1, z + 1);
    }

    if (kstep == 2)                               /* 2x2 pivot block? */
    {
      akm1k    = ap[ik + k - 1];
      ak       = ap[kk]   / akm1k;
      akm1     = ap[ik]   / akm1k;
      bk       = z[k]     / akm1k;
      bkm1     = z[k - 1] / akm1k;
      nenner   = ak * akm1 - ONE;
      z[k]     = (akm1 * bk - bkm1) / nenner;
      z[k - 1] = (ak * bkm1 - bk)   / nenner;
    }

    else                                          /* 1x1 pivot block? */
    {
      if (FABS(z[k]) > FABS(ap[kk]))
      {
        s = FABS(ap[kk]) / FABS(z[k]);
        mul_vector(z + 1, n, s);
        ynorm *= s;
      }
      if (ap[kk] == ZERO)
        z[k] = ONE;
      else
        z[k] /= ap[kk];
    }

    k  -= kstep;
    ik -= k;
    if (kstep == 2)
      ik -= k + 1;
  }

  s = ONE / einsnorm(z + 1, n);                        /* normalize V */
  mul_vector(z + 1, n, s);
  ynorm *= s;


  /* --------------------- Solve  U' * Z = V  ----------------------- */

  ux_b_loesen(ap, n, pvt, z);

  s = ONE / einsnorm(z + 1, n);                        /* normalize Z */
  mul_vector(z + 1, n, s);
  ynorm *= s;


  if (anorm == ZERO)
    *rcond = ZERO;
  else
    *rcond = ynorm / anorm;


#ifdef DEBUG
  fprintf(stderr, "rcond = %"LZP"g\n", *rcond);
#endif
  vmfree(vmblock);
  return 0;
}

/* -------------------------- END linpack.c ------------------------- */
