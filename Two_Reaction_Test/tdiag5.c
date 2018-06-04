#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test 5 diagonal matrix solver diag5                               */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  int  i, j, ud, ld, n, rc, cas;
  REAL *ld2, *ld1, *d, *ud1, *ud2, *b, determ;
  double dum1, dum2, dum3, dum4, dum5;
  void *vmblock;

  /* --- assign an input file to the standard input ----------------- */
  if (argc >= 2)                           /* at least one entry   ?  */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("5 diagonal Matrices");

  if (scanf ("%d", &n)  < 1 ||
      scanf ("%d", &ld) < 1 ||
      scanf ("%d", &ud) < 1)
  {
    LogError ("Input Stream", 0, __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError ("Dimension must be > 0", 0, __FILE__, __LINE__);
    return 1;
  }

  if (ld != 2 || ud != 2)
  {
    LogError ("Wrong input dimensions", 0,  __FILE__, __LINE__);
    return 1;
  }

  vmblock = vminit();
  ld2 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  ld1 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  d   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  ud1 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  ud2 = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  b   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  printf ("\nInput Matrix (condensed):\n\n");
  for (j = 0; j < n; j++)
  {
    if ( (scanf (FORMAT_IN, &dum1) < 1) ||
         (scanf (FORMAT_IN, &dum2) < 1) ||
         (scanf (FORMAT_IN, &dum3) < 1) ||
         (scanf (FORMAT_IN, &dum4) < 1) ||
         (scanf (FORMAT_IN, &dum5) < 1)    )
    {
      LogError ("Input stream", 0,  __FILE__, __LINE__);
      return 1;
    }

    ld2[j] = (REAL) dum1;
    ld1[j] = (REAL) dum2;
    d[j]   = (REAL) dum3;
    ud1[j] = (REAL) dum4;
    ud2[j] = (REAL) dum5;

    printf (FORMAT_LF, ld2[j]);
    printf (FORMAT_LF, ld1[j]);
    printf (FORMAT_LF, d[j]  );
    printf (FORMAT_LF, ud1[j]);
    printf (FORMAT_LF, ud2[j]);
    printf ("\n");
  }

  printf ("\nInverse (transposed):\n\n");

  cas = 0;
  for (i = 0; i < n; i++)
  {
    SetVec (n, b, ZERO);
    b[i] = ONE;

    rc = diag5 (cas, n, ld2, ld1, d, ud1, ud2, b);
    if (rc)
    {
      LogError ("diag5", rc,  __FILE__, __LINE__);
      return 1;
    }

    WriteVec (n, b);
    cas = 2;
  }

  determ = ONE;
  for (i = 0; i < n; i++) determ *= d[i];

  printf ("\nDeterminant = ");
  printf (FORMAT_LE, determ);

  WriteEnd ();

  return (0);
}
