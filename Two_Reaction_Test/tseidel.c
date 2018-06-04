#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test program for Gauss Seidel method                              */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  REAL **a, **c, *b, *x, *residu;
  REAL omega = (REAL) 1.5;
  int n, j, rc, krit = 0, iter;
  void *vmblock;

  /* --- assign the input file to the standard input file ----------- */
  if (argc >= 2)                             /* at least one entry ?  */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Gauss Seidel method");

  if (scanf ("%d ", &n) <= 0)
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
  a      = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  c      = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  b      = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
  x      = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
  residu = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);

  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  if (ReadMat (n, n, a) != 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("Dimension of the input matrix = %d\n", n);
  printf ("Input matrix:\n");

  WriteMat (n, n, a);
  CopyMat (n, n, a, c);

  printf ("\nTransposed Inverse:\n\n");

  for (j = 0; j < n; j++)
  {
    SetVec (n, x, ZERO);
    SetVec (n, b, ZERO);
    b[j] = ONE;
    CopyMat (n, n, c, a);

    if ( (rc = seidel (krit, n, a, b, omega, x, residu, &iter)) != 0)
    {
      LogError ("seidel", rc,  __FILE__, __LINE__);
      return 1;
    }
    WriteVec (n, x);
  }

  printf("\nResidual vector:\n");
  WriteVec (n, residu);

  WriteEnd ();

  return (0);
}
