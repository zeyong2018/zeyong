#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test Gauss method for band matrix                                 */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  int ud, ld, n, i, rc, *perm, dsign, dim, mod;
  REAL **packmat, **save, *b, determ;
  void *vmblock = NULL;

  /* --- assign the input file to the standard input file  ---------- */
  if (argc >= 2)                               /* at least on entry?  */
    if (freopen(argv[1], "r", stdin) == NULL)     /* open input file  */
    {
      fprintf(stderr, "error in opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Banded matrix");

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

  printf ("Dimension: %d\n",      n );
  printf ("Subdiagonals: %d\n",   ld);
  printf ("Superdiagonals: %d\n", ud);

  if (ld + ud >= n)
  {
    LogError ("ld + ud must be < n", 0, __FILE__, __LINE__);
    return (1);
  }

  dim = ld + ud + 1 + min (ld, ud);

  vmblock = vminit();
  packmat = (REAL **)vmalloc(vmblock, MATRIX,  n, dim);
  save    = (REAL **)vmalloc(vmblock, MATRIX,  n, dim);
  b       = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  perm    = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));

  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  if (ReadMat (n, ld + ud + 1, packmat))
  {
    LogError ("Input Stream", 0, __FILE__, __LINE__);
    return 1;
  }

  CopyMat (n, ld + ud + 1, packmat, save);

  printf ("Input matrix:\n");

  WriteMat (n, ld + ud + 1, packmat);

  printf ("\nBand without pivot:\n\n");
  printf ("Inverse (transposed)\n\n");

  mod = 0;
  for (i = 0; i < n; i++)
  {
    SetVec (n, b, ZERO);
    b[i] = ONE;

    rc = bando (mod, n, ld, ud, packmat, b);
    if (rc)
    {
      LogError ("bando", rc, __FILE__, __LINE__);
      return 1;
    }
    WriteVec (n, b);
    mod = 2;
  }

  determ = ONE;
  for (i = 0; i < n; i++) determ *= packmat[i][ld];

  printf ("\nDeterminant = ");
  printf (FORMAT_LE, determ);

  WriteHead ("band with pivot");

  printf ("Inverse (transposed)\n\n");

  mod = 0;
  for (i = 0; i < n; i++)
  {
    SetVec (n, b, ZERO);
    b[i] = ONE;

    rc = band (mod, n, ld, ud, save, b, perm, &dsign);
    if (rc)
    {
      LogError ("band", rc, __FILE__, __LINE__);
      return 1;
    }
    WriteVec (n, b);
    mod = 2;
  }

  determ = (REAL) dsign;
  for (i = 0; i < n; i++) determ *= packmat[i][ld];

  printf ("\nDeterminant = ");
  printf (FORMAT_LE, determ);

  WriteEnd ();

  return (0);
}
