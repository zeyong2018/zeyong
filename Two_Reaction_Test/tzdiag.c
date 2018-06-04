#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*--------------------------------------------------------------------*/
/* Test frame for tzdiag (cyclic tridiagonal matrix)                  */
/*--------------------------------------------------------------------*/

#include <basis.h>
#include <u_proto.h>

#define dim 4


int main (void)
{
  int n = dim, i,j;
  int rep, rc;
  REAL l[dim], d[dim], u[dim], ricol[dim], lowrow[dim],
       b[dim], determ = ONE;


  WriteHead ("Cylic triangular Matrix (using tzdiag)");

  l[0] = l[1] = l[2] = l[3] = -ONE;
  u[0] = u[1] = u[2] = u[3] = -ONE;
  d[0] = d[1] = d[2] = d[3] = (REAL)4.0;
  lowrow[0] = -ONE; ricol[0] = -ONE;

  for (i = 1; i < n - 1; i++)
  {
    l[i] = u[i] = ONE;
    d[i] = ONE;
  }

  d[n-1] = ONE; u[0] = ONE; l[n-1] = ONE; d[0] = (REAL)1.00001;

  printf ("Input Matrix (in condensed form):\n\n");
  for (i = 0; i < n; i++)
  {
    printf (FORMAT_126LF, l[i]);
    printf (FORMAT_126LF, d[i]);
    printf (FORMAT_126LF, u[i]);
    printf ("\n");
  }

  printf ("Inverse:\n\n");
  rep = 0;
  for (j = 0; j < n; j++)
  {
    SetVec (n, b, ZERO);
    b[j] = ONE;
    rc = tzdiag(n, l, d, u, lowrow, ricol, b, rep);
    if (rc)
    {
      LogError ("tzdiag", rc,  __FILE__, __LINE__);
      return 1;
    }
    WriteVec (n, b);
    rep = 1;
  }

  for (determ = ONE, i = 0; i < n; i++) determ *= d[i];

  printf ("\nDeterminant = ");
  printf (FORMAT_LE, determ);

  WriteEnd ();

  return(0);
}
