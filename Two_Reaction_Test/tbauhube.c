#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test program for Bauhuber's method                                */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  REAL *ar, *ai, *rootr, *rooti, *val;
  int i, j, n, rc, skala;
  void *vmblock;

  /* --- assign a possible input file to the standard input file ---- */
  if (argc >= 2)                           /* at least one variable?  */
    if (freopen(argv[1], "r", stdin) == NULL)      /* open input file */
    {
      fprintf(stderr, "Error opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Bauhuber's method");

  if (scanf ("%d", &n) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError ("Polynomial degree must be > 0", 0,  __FILE__, __LINE__);
    return 1;
  }

  vmblock = vminit();
  ar    = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  ai    = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  rootr = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  rooti = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  val   = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);

  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  rc = ReadVec (n + 1, ar);
  if (rc)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  SetVec (n + 1, ai, ZERO);

  printf ("Polynomial coefficients:\n");
  for (i = 0; i <= n; i++)
  {
    printf("\n a[%2d] = ", i);
    printf (FORMAT_2010LF, ar[i]);
    printf (FORMAT_2010LF, ai[i]);
  }

  skala = 0;
  for (j = 0; j < 2; j++)
  {
    rc = bauhub (0, skala, n, ar, ai, rootr, rooti, val);
    if (rc == 0)
    {
      if (skala == 0)
        printf ("\n\nRoots (without scaling)\n");
      else
        printf ("\n\nRoots (with scaling)\n");

      printf ("    No. Real part               Imaginary part         "
              "function value\n");

      for (i = 0; i < n; i++)
      {
        printf (" \n %4d  ", i);
        printf (FORMAT_2016LE, rootr[i]);
        printf (FORMAT_2016LE, rooti[i]);
        printf (FORMAT_LE, val[i]);
      }
      printf ("\n");
    }
    else
      LogError ("bauhube", rc, __FILE__, __LINE__);

    skala = 1;
  }

  WriteEnd ();

  return (0);
}
