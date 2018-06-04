#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test Gauss method                                                 */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  REAL **a, **lu, *b, *x, sum;
  int j, k, n, cas, rc, *perm, signd;
  void *vmblock;

  /* --- assign input file to standard input ------------------------ */
  if (argc >= 2)                           /* at least one entry ?    */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error in opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Gauss Method");

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

  vmblock = vminit();
  a    = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  b    = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  x    = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));
  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  rc = ReadMat (n, n, a);
  if (rc)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("Dimension of the input matrix = %d\n", n);
  printf ("Input matrix:\n");

  WriteMat (n, n, a);
  CopyMat (n, n, a, lu);

  printf ("Transposed Inverse:\n\n");

  sum = ZERO;
  cas = 0;
  for (j = 0; j < n; j++)
  {
    SetVec (n, b, ZERO);
    b[j] = ONE;

    if ((rc = gauss (cas, n, lu, lu, perm, b, x, &signd)) != 0)
    {
      LogError ("gauss", rc,  __FILE__, __LINE__);
      return 1;
    }
    else
    {
      WriteVec (n, x);
      for (k = 0; k < n; k++) sum += a[j][k] * x[k];
    }
    cas = 2;
  }   /* end for j */

  WriteHead ("Gauss Method for multiple right hand sides");

  /* set lu equal to identity matrix .................................*/
  SetMat (n, n, lu, ZERO);
  for (j = 0; j < n; j++) lu[j][j] = ONE;

  rc = mgauss (n, n, a, lu);
  if (rc)
  {
    LogError ("mgauss", rc,  __FILE__, __LINE__);
    return 1;
  }

  printf ("Inverse:\n");
  WriteMat (n, n, lu);

  printf ("Determinant = ");
  printf (FORMAT_LE, det (n, a));
  printf ("\n\nTrace of (Matrix * Inverse) = ");
  printf (FORMAT_2016LE, sum);

  WriteHead ("Solution with iterative refinement steps:");

  rc = gauss (1, n, a, lu, perm, b, x, &signd);    /* Decomp. ........*/
  if (rc)
  {
    LogError ("gauss: mod 1", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("Transposed Inverse:\n\n");
  for (j = 0; j < n; j++)
  {
    SetVec (n, b, ZERO);
    b[j] = ONE;

    rc = gauss (3, n, a, lu, perm, b, x, &signd);    /* Solution .....*/
    if (rc > 0)
    {
      LogError ("gauss: mod 3", 0,  __FILE__, __LINE__);
      return 1;
    }

    if (rc) printf ("rc = %d\n", rc);

    WriteVec (n, x);
  }
  WriteEnd ();

  return 0;
}
