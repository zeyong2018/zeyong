#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test pivot method pivot                                           */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  REAL **a, **inv, sum, cond;
  int n, rc;
  void *vmblock;

  /* --- assign the input file to the standard input file ----------- */
  if (argc >= 2)                           /* at least one entry ?    */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error in opening input file %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Pivot Method");

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

  printf ("Dimension of the input matrix = %d\n", n);
  printf ("Input matrix:\n");

  rc = ReadMat (n, n, a);
  if (rc)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  WriteMat (n, n, a);

  sum = ZERO;
  rc = pivot (n, a, inv, &sum, &cond);
  if (rc == 0)
  {
    printf ("\nInverse:\n");
    WriteMat (n, n, inv);
    printf ("\nCheck sum = ");
    printf (FORMAT_2016LE, sum);
    printf ("\ncondition number = ");
    printf (FORMAT_LF, cond);
  }
  else
  {
    LogError ("pivot", rc,  __FILE__, __LINE__);
    return 1;
  }

  WriteEnd ();

  return (0);
}
