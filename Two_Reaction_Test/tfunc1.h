/* --------------------- DECLARATIONS tfunc1.h ---------------------- */


/*--------------------------------------------------------------------*/
/* Declaration of test functions for newton, pegasus and roots        */
/*                                                                    */
/*                                                                    */
/*--------------------------------------------------------------------*/

#ifndef TFUNC1_H_INCLUDED

#define TFUNC1_H_INCLUDED

REAL   f1   (REAL x);               /* f(x) = x -exp(x) ..............*/
REAL   f1s  (REAL   x);
REAL   f1ss (REAL   x);

REAL   f2   (REAL   x);             /* f(x) = (x-1) hoch 6 ...........*/
REAL   f2s  (REAL   x);
REAL   f2ss (REAL   x);

REAL   f3   (REAL   x);             /* f(x) = sin(x) .................*/
REAL   f3s  (REAL   x);
REAL   f3ss (REAL   x);

REAL   f4   (REAL   x);             /* f(x) = 1.0 + sin(x) ...........*/
REAL   f4s  (REAL   x);
REAL   f4ss (REAL   x);

REAL   f5   (REAL   x);          /* f(x) = exp(x) -(1.0 + x + x*x/2) .*/
REAL   f5s  (REAL   x);
REAL   f5ss (REAL   x);

REAL   f6   (REAL   x);   /* f(x) = sqr(x-1)*(sin(PI*x)-ln(2x/(x+1)) .*/
REAL   f6s  (REAL   x);
REAL   f6ss (REAL   x);

#endif

/* -------------------------- END tfunc1.h -------------------------- */
