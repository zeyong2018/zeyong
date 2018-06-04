#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test program: Numerical differentation according to  Romberg         *
***********************************************************************/
#include <basis.h>
#include <difrom.h>
#include <stdio.h>

REAL Hyperbel (REAL x)
{
  return (1 / x);
}

void aus (int  n,   REAL prec,  REAL h,
          REAL x,   int  error, REAL schaetz,
          REAL res, int  nend,  REAL hend)
{
  printf ("\ninput parameters:\n\n");
  printf ("x     = %20.13"LZP"e\n"
          "eps   = %20.13"LZP"e\n"
          "n     = %2d\n"
          "h     = %20.13"LZP"e\n"
          "error = %d\n",
           x, prec, n, h, error);
  if (error != 1)
  {
    printf ("\nresults:\n\n");
    printf ("er_app= %20.13"LZP"e\n"
            "res   = %20.13"LZP"e\n"
            "nend  = %2d\n"
            "hend  = %20.13"LZP"e\n",
             schaetz, res, nend, hend);
  }
}

#define input(Question,EinFormat,AusFormat,Variable)                  \
   ( printf (Question"\n"),           /* Write out question        */ \
     scanf  (EinFormat,    &Variable),/* read value                */ \
     printf (AusFormat"\n", Variable))/* Reflection for file input */ \

int main (int argc, char *argv[])
{
  int  error, nend, n;
  REAL x, prec, h;
  REAL res, schaetz, hend;

  /* --- assign input fiel to standard file ------------------------- */
  if (argc >= 2)                              /* at least one file ?  */
    if (freopen(argv[1], "r", stdin) == NULL) /* open input file      */
    {
      fprintf(stderr, "error opening %s!\n", argv[1]);
      return 2;
    }

  printf ("Numerical differentation according to Romberg\n\n");
  printf ("Test function   f(x) = 1/x\n\n");

  input ("put in x-value at which you want to evaluate derivative :",
         "%"LZS"f", "%"LZP"f", x);

  input ("put in desired accuracy :", "%"LZS"f", "%"LZP"f", prec);

  input("maximal number of columns in Romberg scheme :", "%d", "%d", n);

  input ("starting step size :", "%"LZS"f", "%"LZP"f", h);

  error = difrom (Hyperbel,x,prec,n,h,&res,&schaetz,&nend,&hend);

  aus (n,prec,h,x,error,schaetz,res,nend,hend);

  return 0;
}

#if 0
T E S T  E X A M P L E
======================

Input:

         0.12
         0.0000000001
         4
         0.005


Output (on 286 Compaq portable III and 386er SICOMP-PC 32-05):

         Numerical differentation according to Romberg

         Test function   f(x) = 1/x

         put in x-value at which you want to evaluate derivative :
         0.120000
         put in desired accuruacy :
         0.000000
         maximal number of columns in Romberg scheme :
         4
         starting step size :
         0.005000

         input parameters:

         x     =  1.2000000000000e-01
         eps   =  1.0000000000000e-10
         n     =  4
         h     =  5.0000000000000e-03
         error =  0

         results:

         er_app=  8.8917317953019e-11
         res   = -6.9444444444444e+01
         nend  =  3
         hend  =  6.2500000000000e-04

#endif
