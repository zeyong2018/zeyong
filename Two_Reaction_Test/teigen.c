#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*----------------------------  TEST for eigen  ----------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>

/*====================================================================*/
/*  Calculation of eigenvalues and eigenvectors                       */
/*====================================================================*/
/*    Call   :  teigen < mat                                          */
/*                                                                    */
/*  As an example, the input file has the form:                       */
/*  5                                                                 */
/*  1.0 2.0 3.0  1.0  2.0                                             */
/*  0.1 0.3 0.0  2.0  1.0                                             */
/*  3.0 1.0 9.0  4.0  2.0                                             */
/*  0.0 0.0 0.0 -1.0  0.0                                             */
/*  0.1 0.1 .01  7.0  1.0                                             */
/*                                                                    */
/*  The first number (5) is the dimension of the matrix followed by   */
/*  the 1st to 5th row of the 5 x 5 matrix (in our case).             */
/*====================================================================*/

int main (int argc, char *argv[])
{
  REAL   **mat,              /* input  matrix                         */
         **a,                /* copy of the input matrix              */
         **ev = NULL,        /* Eigenvectors if vec <> 0              */
         *wr,                /* Eigenvalues (Real part)               */
         *wi;                /* Eigenvalues (Imaginary parts)         */

  int n,
      *cnt,                /* Iteration counter                       */
      rc,                  /* Return Code                             */
      vec = 1,             /* flag for eigenvectors (=0 -> none)      */
      fehler,              /* error code from `umleiten'              */
      ortho = 0,           /* flag for orthogonal                     */
                           /* Hessenberg reduction                    */
      ev_norm = 1;         /* flag for normalization of Eigenvectors  */

  register i, j, k;
  REAL   v, w, norm;
  void   *vmblock;

  if ((fehler = umleiten(argc, argv))   /* assign input and output    */
      != 0)                             /* file to standard input and */
    return fehler;  /* 1 or 2 */        /* standard output            */

  WriteHead ("Eigenvalues and Eigenvectors");

  if (scanf ("%d", &n) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError ("Dimension must be > 0", 0,  __FILE__, __LINE__);
    return 1;
  }

  /* Allocate Memory .................................................*/

  vmblock = vminit();
  mat  = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  a    = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  wr   = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  wi   = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  cnt  = (int  *) vmalloc(vmblock, VVEKTOR, n, sizeof(*cnt));

  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  if (vec)                              /* Only for eigenvectors .....*/
  {
    ev = (REAL **)vmalloc(vmblock, MATRIX, n, n);
    if (! vmcomplete(vmblock))
    {
      LogError ("No Memory", 0, __FILE__, __LINE__);
      return 1;
    }
  }

  rc = ReadMat (n, n, mat);
  if (rc)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("\nDimension of the input matrix = %d\n", n);
  printf ("Input matrix:\n");

  WriteMat (n, n, mat);
  CopyMat (n, n, mat, a);

  printf("\n");

  rc = eigen (vec, ortho, ev_norm, n, mat, ev, wr, wi, cnt);

  if (rc != 0)                                  /*  ERROR !!!         */
  {
    LogError ("eigen", rc,  __FILE__, __LINE__);
    for (i = 0; i < n; i++) printf (" %d ", cnt[i]);
    return 1;
  }
                                 /*  If vec != 0, print eigenvectors  */
  if (vec)
  {
    if (ev_norm)
      printf("Normalized Eigenvectors:\n");
    else
      printf("not normalized Eigenvectors:\n");
    WriteMat (n, n, ev);
  }

  printf ("\nEigenvalues:\t\t\t\t\t\tIterations:\n\n");
  for (i = 0; i < n; i++)
  {
    printf (FORMAT_2016LE, wr[i]);
    printf ("+ ");
    printf (FORMAT_2016LE, wi[i]);
    printf ("* i\t%4d\n", cnt[i]);
  }
                      /* Check result: sum of L1 norms of             */
                      /* Matrix*Eigenvector - Eigenvalue*Eigenvector  */
  if (vec)            /* (this must be nearly 0).                     */
  {
    for (norm = ZERO, k = 0; k < n; k++)
    {
      if (wi[k] == ZERO)
      {
        for (i = 0; i < n; i++)
        {
          for (w = ZERO, j = 0; j < n; j++)
            w += a[i][j] * ev[j][k];
          w -= wr[k] * ev[i][k];
          norm += ABS (w);
        }
      }
      else
      {
        for (i = 0; i < n; i++)
        {
          for (w = ZERO, j = 0; j < n; j++)
            w += a[i][j] * ev[j][k];
          w -= wr[k] * ev[i][k] - wi[k] * ev[i][k+1];
          for (v = ZERO, j = 0; j < n; j++)
            v += a[i][j] * ev[j][k+1];
          v -= wr[k] * ev[i][k+1] + wi[k] * ev[i][k];
          norm += 2.0 * SQRT (v*v + w*w);
        }
        k++;
      }
    }
    printf ("\nCheck sum = ");
    printf (FORMAT_LE, norm);
    printf ("(must be approximately 0)\n\n");
  }

  WriteEnd ();

  return (0);
}
