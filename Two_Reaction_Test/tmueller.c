#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test program for mueller's method                                 */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  REAL *a, *zreal, *zimag, fre, fim;
  int i, j, n, rc, scale;
  void *vmblock;

  /* --- assign an eventuel input file to the standard input file  -- */
  if (argc >= 2)                           /* at least one variable?  */
    if (freopen(argv[1], "r", stdin) == NULL)      /* open input file */
    {
      fprintf(stderr, "Error opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Mueller's Method");

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
  a     = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  zreal = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);
  zimag = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);

  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  rc = ReadVec (n + 1, a);
  if (rc)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("Coefficients (ascending): \n");

  for (i = 0; i <= n; i++)
  {
    printf("\n a[%2d] = ", i);
    printf (FORMAT_2010LF, a[i]);
  }

  scale = 0;
  for (j = 0; j < 2; j++)
  {
    rc = mueller (n, a, scale, zreal, zimag);
    if (rc == 0)
    {
      if (scale == 0)
        printf ("\n\nWithout scaling\n");
      else
        printf ("\n\nWith scaling\n");

      printf("Roots  Real part                 Imaginary part\t\tf "
             "value\n");
      for (i = 0; i < n; i++)
      {
        fmval (n, 0, a, a[n], zreal[i], zimag[i], &fre, &fim);
        printf (" \n %4d  ", i);
        printf (FORMAT_2016LE, zreal[i]);
        printf (FORMAT_2016LE, zimag[i]);
        printf (FORMAT_LE, fre);
      }
      printf("\n");
    }
    else
    {
      LogError ("mueller", rc,  __FILE__, __LINE__);
      return 1;
    }
    scale = 1;
  }

  WriteEnd ();

  return (0);
}
