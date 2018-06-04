#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test program for newmod                                           */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <u_proto.h>
#include <tfunc1.h>


int main (void)
{
  int rc, j, iter, i;
  REAL x, f;

  WriteHead ("Newmod: modified Newton method");

  for (i = 1; i <= 6; i++)
  {
    x = TWO;

    switch (i)
    {
      case 1 :   rc = newmod (f1, f1s, f1ss, &x, &f, &iter, &j); break;
      case 2 :   rc = newmod (f2, f2s, f2ss, &x, &f, &iter, &j); break;
      case 3 :   rc = newmod (f3, f3s, f3ss, &x, &f, &iter, &j); break;
      case 4 :   rc = newmod (f4, f4s, f4ss, &x, &f, &iter, &j); break;
      case 5 :   rc = newmod (f5, f5s, f5ss, &x, &f, &iter, &j); break;
      case 6 :   rc = newmod (f6, f6s, f6ss, &x, &f, &iter, &j); break;
      default:   printf ("Invalid Function\n"); return (1);
    }

    printf ("\n Function number: %d", i);
    printf ("\n Return code:             % d",  rc);
    printf ("\n Root x:                  "); printf (FORMAT_2016LF, x);
    printf ("\n with multiplicity:       % d",  j);
    printf ("\n function value at x:     "); printf (FORMAT_LE,  f);
    printf ("\n Iterations:              % d",  iter);
    printf("\n");
  }
  WriteEnd ();

  return (0);
}
