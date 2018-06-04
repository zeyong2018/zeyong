#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test program for newpoly                                          */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


int main (int argc, char *argv[])
{
  int  rc, n, iter;
  REAL x, f, *a;
  void *vmblock;

  /* --- assign a possible input file to the standard input       -- */
  if (argc >= 2)                           /*  at least one argument?*/
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file     */
    {
      fprintf(stderr, "Error opening  %s!\n", argv[1]);
      return 2;
    }

  WriteHead ("Newton method for polynomials");

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
  a = (REAL *)vmalloc(vmblock, VEKTOR, n + 1, 0);

  rc = ReadVec (n + 1, a);
  if (rc)
  {
    LogError ("Input stream", 0,  __FILE__, __LINE__);
    return 1;
  }

  printf ("Polynomial coefficients (ascending):\n");
  WriteVec (n, a);

  x = (REAL) 0.9;

  printf ("\nStarting value:    "); printf (FORMAT_LF, x);

  rc = newpoly (n, a, &x, &f, &iter);

  printf ("\nReturn code:       % d\n", rc);
  printf ("Root:              "); printf (FORMAT_2016LF, x);
  printf ("\nFunction value:    "); printf (FORMAT_LE, f);
  printf ("\nIterations:        % d\n", iter);

  WriteEnd ();

  return (0);
}
