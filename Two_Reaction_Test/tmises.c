#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*====================================================================*/
/*  Test program for fmises                                           */
/*====================================================================*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  REAL ew, **c, *x;
  int n, rc;
  void *vmblock;

  /* --- assign an input file to the standard input file ------------ */
  if (argc >= 2)                           /* at least one entry ?    */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "Error when opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Vector iteration method");

  if (scanf ("%d", &n) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError ("Dimension must be > 0", 0, __FILE__, __LINE__);
    return 1;
  }

  vmblock = vminit();
  c = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  x = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  rc = ReadMat (n, n, c);
  if (rc)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("\nDimension of the input matrix = %d \n",n);
  printf ("Input matrix:\n");
  WriteMat (n, n, c);

  rc = mises (n, c, x, &ew);

  if (rc)                                       /*  ERROR !!!         */
  {
    LogError ("mises", rc,  __FILE__, __LINE__);
    return 1;
  }

  printf ("\nEigenvector:\n");

  WriteVec (n, x);
  printf ("\nMaximum modulus eigenvalue:\n");
  printf (FORMAT_2016LF, ew);

  WriteEnd ();

  return (0);
}
