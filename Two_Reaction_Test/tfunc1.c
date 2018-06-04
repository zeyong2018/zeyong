#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*--------------------------------------------------------------------*/
/* Definition of test functions for newton, pegasus and roots         */
/* algorithms                                                         */
/*                                                                    */
/*--------------------------------------------------------------------*/

#include <basis.h>
#include <u_proto.h>
#include <tfunc1.h>


REAL f1 (REAL x)             /* f(x) = 5*x -exp(x)  */
{
  return ((REAL)5.0 * x - EXP (x));
}

REAL f1s (REAL x)
{
  return ((REAL)5.0 - EXP (x));
}

REAL f1ss (REAL x)
{
  return (-EXP(x));
}
/*--------------------------------------------------------------------*/

REAL f2 (REAL x)                     /* f(x) = (x-1) hoch 6  */
{ return ( (((((x-6)*x+15)*x-20)*x+15)*x-6)*x+1 ); }

REAL f2s (REAL x)
{ return( ((((6*x-30)*x+60)*x-60)*x+30)*x-6 ); }

REAL f2ss (REAL x)
{ return ( (((30*x-120)*x+180)*x-120)*x+30 ); }

/*--------------------------------------------------------------------*/

REAL f3 (REAL x)              /* f(x) = sin(x)  */
{
  return (SIN(x));
}

REAL f3s (REAL x)
{
  return (COS(x));
}

REAL f3ss (REAL x)
{
  return (-SIN(x));
}

/*--------------------------------------------------------------------*/

REAL f4 (REAL x)              /* f(x) = 1.0 + sin(x)  */
{
  return (ONE + SIN (x));
}

REAL f4s (REAL x)
{
  return (COS(x));
}

REAL f4ss (REAL x)
{
  return (-SIN(x));
}

/*--------------------------------------------------------------------*/

REAL f5 (REAL x)      /* f(x) = exp(x) -(1.0 + x + x*x/2)  */
{
  return (EXP (x) - (ONE + x + x * x * HALF));
}

REAL f5s (REAL x)
{
  return (EXP (x) - ONE - x);
}

REAL f5ss (REAL x)
{
  return (EXP (x) - ONE);
}

/*--------------------------------------------------------------------*/

#ifndef PI
#define PI ((REAL)4.0 * ATAN (1.0))
#endif

REAL f6 (REAL x)    /* f(x) = sqr(x-1)*(SIN(PI*x)-ln(2x/(x+1))  */
{
  REAL f, temp;

  temp = TWO * ABS ((x / (x + ONE)));
  f = (x - ONE) * (x - ONE) * (SIN (PI * x)  - LOG (temp));
  return (f);
}

REAL f6s  (REAL x)
{
  REAL f, temp;

  temp = TWO * ABS ((x / (x + ONE)));

  f = TWO * (x - ONE) * (SIN (PI * x) - LOG (temp));

  f += (x - ONE) * (x - ONE) * (PI * COS (PI * x)
       - ONE / x / (x + ONE));

  return (f);
}

REAL f6ss (REAL x)
{
  REAL f, temp;

  temp = TWO * ABS ((x / (x + ONE)));

  f = TWO * (SIN (PI * x) - LOG (temp));

  f += (REAL)4.0 * (x - ONE) * (PI * COS (PI * x)
        - ONE / x / (x + ONE));

  f += (x - ONE) * (x - ONE)
        * ((TWO * x + ONE) / x / x / (x + ONE) / (x + ONE)
        - PI * PI * SIN (PI * x));

  return (f);
}
