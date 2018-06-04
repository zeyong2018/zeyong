#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*-------------------------------------------------------------------*/
/* Test unpack                                                       */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>

int main (void)
{
  int ud , ld, n, i, j, rc, dim;
  REAL **packmat, *row;
  double tmp;
  void *vmblock;

  WriteHead ("Test unpacking a condensed matrix");

  if (scanf("%d", &n) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError ("n must be > 0", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (scanf("%d", &ld) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (ld < 0)
  {
    LogError ("ld must be > 0", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (scanf("%d", &ud) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (ud < 0)
  {
    LogError ("ud must be > 0", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (ld + ud >= n)
  {
    LogError ("ld + ud must be < n", 0, __FILE__, __LINE__);
    return (1);
  }

  dim = ld + ud + 1;

  vmblock = vminit();
  packmat = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  row     = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < dim; j++)
    {
      scanf(FORMAT_IN, &tmp);
      packmat[i][j] = (REAL) tmp;
    }
  }

  printf ("%d %d %d\n", n, ld, ud);
  for (i = 0; i < n; i++)
  {
    rc = unpack (n, ld, ud, i, packmat[i], row);
    if (rc)
    {
      LogError ("unpack", rc, __FILE__, __LINE__);
      return 1;
    }
    WriteVec (n, row);
  }

  WriteEnd ();

  return (0);
}
