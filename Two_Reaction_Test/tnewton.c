#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*-------------------------------------------------------------------*/
/* Test program for newton                                           */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <u_proto.h>
#include <tfunc1.h>

int main (void)
{
  int i, rc, iter;
  REAL x, f;
  REAL (*fct) (REAL), (*fctd) (REAL);
  char *text;

  WriteHead ("Newton method for real valued, nonlinear functions");

  for (i = 1; i <= 6; i++)
  {
    switch (i)
    {
      case 1: x = ONE;
              fct  = f1;
              fctd = f1s;
              text = "f = 5 * x - exp (x)";
              break;
      case 2: x = TWO;
              fct  = f2;
              fctd = f2s;
              text = "f = (((((x-6)*x+15)*x-20)*x+15)*x-6)*x+1";
              break;
      case 3: x = ONE;
              fct  = f3;
              fctd = f3s;
              text = "f = sin (x)";
              break;
      case 4: x = ONE;
              fct  = f4;
              fctd = f4s;
              text = "f = 1 + sin (x)";
              break;
      case 5: x = TWO;
              fct  = f5;
              fctd = f5s;
              text = "f = exp(x) -(1.0 + x + x*x*0.5)";
              break;
      case 6: x = TWO;
              fct  = f6;
              fctd = f6s;
              text = "f = (x-1.0)*(x-1.0)*( sin(PI*x)  - "
                     "log(fabs(2.0*x/(x+1.0)))";
              break;
      default: return (1);
    }

    printf ("%s\n", text);
    printf ("Starting value x0 = ");
    printf (FORMAT_LF, x);
    printf ("\n");

    rc = newton (fct, fctd, &x, &f, &iter);

    printf ("Return code       = % d\n", rc);
    printf ("Root              = "); printf (FORMAT_2016LF, x);
    printf ("\nFunction value    = "); printf (FORMAT_LE, f);
    printf ("\nIterations        = % d\n\n", iter);
  }

  WriteEnd ();
  return (0);
}
