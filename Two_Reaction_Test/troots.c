#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test program for roots                                            */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <u_proto.h>
#include <tfunc1.h>


int main (void)
{
  REAL (*func)(REAL), x1, x2, ff2;
  int method, rc, iter, quadex = 1, fctno;

  WriteHead ("Root finder");

  method = 4;

  for (fctno = 1; fctno <= 6; fctno++)
  {
    printf (" Function number = %d\n", fctno);
    switch (fctno)
    {
      case 1: x1 =   (REAL)0; x2 = (REAL)2; func = f1; break;
      case 2: x1 =   (REAL)0; x2 = (REAL)2; func = f2; break;
      case 3: x1 =  -(REAL)1; x2 = (REAL)2; func = f3; break;
      case 4: x1 =   (REAL)0; x2 = (REAL)2; func = f4; break;
      case 5: x1 =  -(REAL)1; x2 = (REAL)2; func = f5; break;
      case 6: x1 = (REAL)0.5; x2 = (REAL)2; func = f6; break;
      default: printf ("Invalid function\n"); return (1);
    }

    rc = roots (method, func, quadex, &x1, &x2, &ff2, &iter);


    printf (" Return code:             % d\n", rc);
    printf (" Root:                    "); printf (FORMAT_2016LF, x2);
    printf ("\n Function value           "); printf (FORMAT_LE, ff2);
    printf ("\n Iterations:              % d\n\n", iter);
  }

  WriteEnd ();
  return (0);
}
