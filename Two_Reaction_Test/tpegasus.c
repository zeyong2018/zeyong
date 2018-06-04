#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*-------------------------------------------------------------------*/
/* Test programm for fpegasus                                        */
/* The module tfunc1 must be included                                */
/*-------------------------------------------------------------------*/

#include <basis.h>
#include <u_proto.h>
#include <tfunc1.h>


int main (void)
{
  REAL x1, x2, ff2;
  int rc, iter, fno;

  WriteHead ("Pegasus method");

  for (fno = 1; fno <= 6; fno++)
  {
    printf (" Function number = %d\n", fno);
    switch (fno)
    {
      case 1:   x1 = (REAL)0; x2 = (REAL)2;
                rc = pegasus (f1, &x1, &x2, &ff2, &iter);
                break;
      case 2:   x1 = (REAL)0; x2 = (REAL) 2;
                rc = pegasus (f2, &x1, &x2, &ff2, &iter);
                break;
      case 3:   x1 = -(REAL)1; x2 = (REAL)2;
                rc = pegasus (f3, &x1, &x2, &ff2, &iter);
                break;
      case 4:   x1 = (REAL)0; x2 = (REAL)2;
                rc = pegasus (f4, &x1, &x2, &ff2, &iter);
                break;
      case 5:   x1 = -(REAL)1; x2 = (REAL)2;
                rc = pegasus (f4, &x1, &x2, &ff2, &iter);
                break;
      case 6:   x1 = (REAL)0.5; x2 = (REAL)2;
                rc = pegasus (f4, &x1, &x2, &ff2, &iter);
                break;
      default:  return (1);
    }

    printf (" Return code:             % d\n", rc);
    printf (" Root:                    "); printf (FORMAT_LF, x2);
    printf ("\n Function value           "); printf (FORMAT_LE, ff2);
    printf ("\n Iterations:              % d\n\n", iter);
  }

  WriteEnd ();
  return (0);
}
