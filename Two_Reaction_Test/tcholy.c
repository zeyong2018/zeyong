#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test Cholesky Method (fcholy)                                     */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


/* Only diagonal and sub-diagonal elements will allocated and read in */


int main (int argc, char *argv[])
/*--------------------------------------------------------------------*/
/* Test Cholesky method for symmetric and positive definite matrices  */
/*--------------------------------------------------------------------*/
{
  REAL d, **a, *b, *x, **org, trace;
  double tmp;
  int n, i, j, cas, rc;
  void *vmblock;

  /* --- assign input file to standard input file ------------------- */
  if (argc >= 2)                           /* at least one entry?     */
    if (freopen(argv[1], "r", stdin) == NULL)     /* open input file  */
    {
      fprintf(stderr, "error opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Cholesky Method");

  if (scanf ("%d", &n) <= 0)
  {
    LogError ("Input Stream", 0, __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError ("n must be > 0", 0, __FILE__, __LINE__);
    return 1;
  }

  vmblock = vminit();
  a   = (REAL **)vmalloc(vmblock, UMATRIX, n, 0);
  org = (REAL **)vmalloc(vmblock, UMATRIX, n, 0);
  x   = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  b   = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);

  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  for (i = 0; i < n; i++)
  {

    for (j = 0; j <= i; j++)
    {
      if (scanf (FORMAT_IN, &tmp) <= 0)
      {
        LogError ("Input Stream", 0, __FILE__, __LINE__);
        return 1;
      }
      org[i][j] = (REAL)tmp;
    }
  }

  printf("Input Matrix:\n\n");
  for (i = 0; i < n; i++)
  {

    for (j = 0; j <= i; j++)
    {
      a[i][j] = org[i][j];
      printf (FORMAT_126LF, a[i][j]);
    }
    printf("\n");
  }

  printf ("\nInverse:\n\n");

  cas = 0;
  trace = ZERO;

  for (j = 0; j < n; j++)
  {
    SetVec (n, b, ZERO);
    b[j] = ONE;
    rc = choly (cas, n, a, b, x);
    if (rc == 0)
    {
      for (i = 0; i < n; i++)
      {
        printf (FORMAT_126LF, x[i]);
        if (j >= i)
          trace += org[j][i] * x[i];
        else
          trace += org[i][j] * x[i];
      }
    }
    else
    {
      LogError ("choly", rc,  __FILE__, __LINE__);
      return 1;
    }
    cas = 2;
    printf ("\n");
  }

  for (d = ONE, i = 0; i < n; i++) d *= SQR (a[i][i]);

  printf("\nDeterminant = ");
  printf (FORMAT_LE, d);
  printf ("\nTRACE of (matrix * inverse) = ");
  printf (FORMAT_2016LE, trace);

  WriteEnd ();

  return (0);
}
