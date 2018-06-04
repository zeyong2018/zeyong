#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*--------------------------------------------------------------------*/
/* Test frame for tridiagonal matrices (tdiag)                        */
/*--------------------------------------------------------------------*/


#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (void)
{
  int n, i,j;
  int rep, rc;
  REAL *l, *d, *u, *b, determ = ONE;
  void *vmblock;

  WriteHead ("Triangular Matrix (tdiag)");

  n = 10;
  vmblock = vminit();
  l = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  d = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  u = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  b = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  for (i = 0; i < n; i++)
  {
    l[i] = u[i] = ONE;
    d[i] = (REAL)10.0;
  }

  printf ("Input Matrix (in condensed form):\n\n");
  for (i = 0; i < n; i++)
  {
    printf (FORMAT_126LF, l[i]);
    printf (FORMAT_126LF, d[i]);
    printf (FORMAT_126LF, u[i]);
    printf ("\n");
  }

  printf ("\nInverse:\n\n");
  rep = 0;
  for (j = 0; j < n; j++)
  {
    SetVec (n, b, ZERO);
    b[j] = ONE;
    rc = trdiag (n, l, d, u, b, rep);
    if (rc)
    {
      LogError ("trdiag", rc,  __FILE__, __LINE__);
      return 1;
    }

    WriteVec (n, b);
    rep = 1;
  }

  for (i = 0; i < n; i++) determ *= d[i];

  printf ("\nDeterminant = ");
  printf (FORMAT_LE, determ);

  WriteEnd ();

  return (0);
}
