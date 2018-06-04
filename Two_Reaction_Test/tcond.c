#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*-------------------------------------------------------------------*/
/* Test estimation programs for condition number                     */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  REAL **a, **inv, suma = ZERO, sumi = ZERO, norma, normi;
  int i, j, n, rc;
  void *vmblock;

  /* --- assign the input file to standard input file --------------- */
  if (argc >= 2)                             /* at least one entry ?  */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error opening %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Condition numbers:");

  if (scanf ("%d", &n) <= 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  if (n < 1)
  {
    LogError ("Dimension must be > 0", 0,  __FILE__, __LINE__);
    return 1;
  }

  vmblock = vminit();
  a   = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  inv = (REAL **)vmalloc(vmblock, MATRIX, n, n);

  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  rc = ReadMat (n, n, a);
  if (rc != 0)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("Dimension of the input matrix = %d\n", n);
  printf ("Input matrix:\n");

  WriteMat (n, n, a);

  /* use mgauss to determine the inverse .............................*/
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (i == j) inv[i][j] = ONE;
        else inv[i][j] = ZERO;

  rc = mgauss (n, n, a, inv);
  if (rc)
  {
    LogError ("mgauss", rc,  __FILE__, __LINE__);
    return 1;
  }

  /* Calculate |A| * |INVERSE of A|, | | means the maximum norm ......*/
  norma = ZERO;
  normi = ZERO;
  for (i = 0; i < n; i++)
  {
    suma = ZERO;
    sumi = ZERO;
    for (j = 0; j < n; j++)
    {
      suma += ABS (a[i][j]);
      sumi += ABS (inv[i][j]);
    }
    if (suma > norma) norma = suma;
    if (sumi > normi) normi = sumi;
  }

  printf ("Calcutated |A| * |INVERSE of A| = ");
  printf (FORMAT_LE, suma * sumi);
  printf ("\nEstimations:\n");
  printf ("  Cline estimate            = ");
  printf (FORMAT_LE, ccond (n, a));
  printf ("\n  Forsythe/Moler estimate   = ");
  printf (FORMAT_LE, fcond (n, a));
  printf ("\n\nHadamard condition number     = ");
  printf (FORMAT_LE, hcond (n, a));
  printf ("\n\n");

  WriteEnd ();

  return (0);
}
