#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*-------------------------------------------------------------------*/
/* Test Householder's Method                                         */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>

int main (int argc, char *argv[])
{
  REAL **a, **inv, *b, sum;
  int n, i, j, k, rc;
  void *vmblock;

  /* --- assign the input file to the standard input file ----------- */
  if (argc >= 2)                           /* at least one entry   ?  */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error in opening %s!\n", argv[1]);
      return 2;
    }

  sum = ZERO;
  WriteHead ("Householder Method");

  if (scanf ("%d", &n) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  vmblock = vminit();
  a   = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  inv = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  b   = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
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

  printf ("\nTransposed Inverse:\n\n");

  for (j = 0; j < n; j++)
  {
    CopyMat (n, n, a, inv);

    SetVec (n, b, ZERO);
    b[j] = ONE;

    rc = house (n, n, inv, b);
    if (rc)
    {
      LogError ("house", rc,  __FILE__, __LINE__);
      return 1;
    }

    WriteVec (n, b);
    for (k = 0; k < n; k++) sum += a[j][k] * b[k];
  }

  printf ("\nTrace of (Matrix*Inverse) = ");
  printf (FORMAT_2016LE, sum);
  printf ("\n\n");

  SetMat (n, n, inv, ZERO);
  for (i = 0; i < n; i++) inv[i][i] = ONE;

  rc = mhouse (n, n, n, a, inv);
  if (rc)
  {
    LogError ("mhouse", rc,  __FILE__, __LINE__);
    return 1;
  }

  printf ("\nHouseholder method with nonunique solutions\n\n");
  printf ("Inverse:\n");

  WriteMat (n, n, inv);

  WriteEnd ();

  return (0);
}
