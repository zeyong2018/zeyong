#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Here we define examples for systems of real functions for linear     *
* approximation together with an alpha-numerical description of them.  *
*                                                                      *
* When running the main program and using the function                 *
* linansf_waehlen(), one can choose between various examples of the    *
* type bsptyp2. The individual examples are declared as global for the *
* module.                                                              *
*                                                                      *
* The pointer phi inside  bsptyp2 points to a C function after call    *
* of linansf_waehlen(), which computes the functional value for the    *
* test function numbered by  i inside the set of test functions phi () *
*                                                                  j   *
* at the node x and return the functional value there.                 *
* The pointer  phi_text gives the address of the string which describes*
* the system of functions alpha-numerically.                           *
* The number of functions in the test system is given by n.            *
* If n < 0, at least in theory, one might call infinitely many such    *
* test functions.                                                      *
***********************************************************************/

#include <basis.h>      /*  for  POW, EXP, COS, SIN, fprintf, stderr, */
                        /*       exit, REAL, ZERO                     */
#include <happrans.h>   /*  for  bsptyp2, linansf_waehlen             */



/***********************************************************************
* Example  0:                                                          *
* ===========                                                          *
* test functions phi (.): Monomials    phi (x) = x^j,  j = 0,1,2,...   *
*                   j                     j                            *
***********************************************************************/
static REAL phi0(int i, REAL x)
{
  return POW(x, (REAL)i);
}

static char *phi0_text(void) /* alpha-numerical description of system */
{
  return
    "PHI(x) = c0 + c1 * x + c2 * x^2 + ... + cn * x^n\n";
}



/***********************************************************************
* Example  1:                                                          *
* ===========                                                          *
* test functions   phi (.): Exponential powers     phi (x) = exp(j*x), *
*                     j                               j                *
*                                                  j = 0,1,2,...       *
***********************************************************************/
static REAL phi1(int i, REAL x)
{
  return EXP(i * x);
}

static char *phi1_text(void)
{
  return
    "PHI(x) = c0 + c1 * exp(x) + c2 * exp(2*x) + ... + cn * exp(n*x)\n";
}



/***********************************************************************
* Example  2:                                                          *
* ===========                                                j         *
* test functions   phi (.): Sine powers      phi (x) = sin(x) ,        *
*                     j                         j                      *
*                                            j = 0,1,2,...             *
***********************************************************************/
static REAL phi2(int i, REAL x)
{
  return POW(SIN(x), (REAL)i);
}

static char *phi2_text(void)      /* system of test functions as text */
{
  return
    "PHI(x) = c0 + c1 * sin(x) + c2 * sin(x)^2 + ... + cn * sin(x)^n\n";
}



/***********************************************************************
* Example  3:                                                          *
* ===========                                                          *
*                                                              j       *
* test dunctions   phi (.): Cosine powers      phi (x) = cos(x) ,      *
*                     j                           j                    *
*                                              j = 0,1,2,...           *
***********************************************************************/
static REAL phi3(int i, REAL x)
{
  return POW(COS(x), (REAL)i);
}

static char *phi3_text(void)      /* system of test functions as text */
{
  return
    "PHI(x) = c0 + c1 * cos(x) + c2 * cos(x)^2 + ... + cn * cos(x)^n\n";
}



/***********************************************************************
* Example  4:                                                          *
* ===========                                                          *
* two separate functions:             phi (x) = x*x, phi (x) = exp(x)  *
*                                        0              1              *
***********************************************************************/
static REAL phi4(int i, REAL x)
{
  switch (i)
  {
    case 0:  return x * x;
    case 1:  return EXP(x);
    default: fprintf(stderr, "error in phi4(), happrans.c:  "
                             "one test function not defined %d!\n", i);
             exit(99);
             return ZERO;
  }
}

static char *phi4_text(void)      /* system of test functions as text */
{
  return
    "PHI(x) = c0 * x + c1 * exp(x)\n";
}



static bsptyp2 beispiel[] =   /* Vector, to address all above examples*/
  {{ 0, phi0, phi0_text },
   { 0, phi1, phi1_text },
   { 0, phi2, phi2_text },
   { 0, phi3, phi3_text },
   { 2, phi4, phi4_text }
  };



/***********************************************************************
* Using this function, one of the examples in beipiel can be selected  *
* for testing.                                                         *
* If the variable nummer lies inside the index set for beispiel, then  *
* the address of vector component is returned; otherwise a zero return *
* signifies the error.                                                 *
***********************************************************************/

bsptyp2 *linansf_waehlen(unsigned int nummer)

{
  if (nummer >= sizeof(beispiel) / sizeof(*beispiel))
    return NULL;                         /* invalid example called    */

  return &beispiel[nummer];
}
