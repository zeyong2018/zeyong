#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* Test program for Brown's method                                      *
*                                                                      *
* This program defines several examples of nonlinear systems of        *
* equations. Specifically we describe a C function for computing as    *
* well as a C function which verbally describes the system.            *
*                                                                      *
* By using the function nlgls_waehlen(), which supplies a pointer for  *
* an object of the proper type bsptyp1, one can select one of these    *
* for a run in the calling program. The example functions are realized *
* here as global only for this module.                                 *
*                                                                      *
* The pointer fkt inside an object of the type bsptyp1 designates      *
* (upon calling  nlgls_waehlen() ) a C function, which computes the    *
* value f of the kth function of the system at x.                      *
* The function pointer fkt_text designates the address of a string,    *
* that describes the system in readable form.                          *
* The number of equations of the example is given by n.                *
* When calling  fkt, k may have values 0,..,n-1, otherwise fkt returns *
* an error code.                                                       *
***********************************************************************/

#include <basis.h>     /*  for  REAL, ZERO, SIN, COS, PI, EIGHT, sqr, */
                       /*       FOUR, LOG, ONE, TWO, HALF, MACH_EPS,  */
                       /*       TEN, POW, SIX, EXP, FIVE, THREE, NULL */
#include <brown.h>     /*  for  nlgls                                 */
#include <nlglstst.h>  /*  for  bsptyp1, nlgls_waehlen                */



/* ------------------------------------------------------------------ */

static int fkt0(int k, REAL x[], REAL *f);
static int fkt1(int k, REAL x[], REAL *f);
static int fkt2(int k, REAL x[], REAL *f);
static int fkt3(int k, REAL x[], REAL *f);
static int fkt4(int k, REAL x[], REAL *f);
static int fkt5(int k, REAL x[], REAL *f);
static int fkt6(int k, REAL x[], REAL *f);
static int fkt7(int k, REAL x[], REAL *f);
static int fkt8(int k, REAL x[], REAL *f);
static int fkt9(int k, REAL x[], REAL *f);
static int fkt10(int k, REAL x[], REAL *f);
static char *fkt0_text(void);
static char *fkt1_text(void);
static char *fkt2_text(void);
static char *fkt3_text(void);
static char *fkt4_text(void);
static char *fkt5_text(void);
static char *fkt6_text(void);
static char *fkt7_text(void);
static char *fkt8_text(void);
static char *fkt9_text(void);
static char *fkt10_text(void);

static bsptyp1 beispiel[] =              /* Vector, that registers all*/
  {{ 2, fkt0,  fkt0_text  },             /* examples below            */
   { 3, fkt1,  fkt1_text  },
   { 2, fkt2,  fkt2_text  },
   { 4, fkt3,  fkt3_text  },
   { 2, fkt4,  fkt4_text  },
   { 2, fkt5,  fkt5_text  },
   { 2, fkt6,  fkt6_text  },
   { 2, fkt7,  fkt7_text  },
   { 2, fkt8,  fkt8_text  },
   { 2, fkt9,  fkt9_text  },
   { 2, fkt10, fkt10_text },
  };



/* ------------------------------------------------------------------ */

static int fkt0(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 0                               *
***********************************************************************/

{
  switch (k)
  {
    case 0:  *f = sqr(x[0]) - x[1] - ONE;
             break;
    case 1:  *f = sqr(x[0] - TWO) + sqr(x[1] - HALF) - ONE;
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt0_text(void)

{
  return
    "x0^2 - x1 - 1                  =  0\n"
    "(x0 - 2)^2 + (x1 - 0.5)^2 - 1  =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt1(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 1                               *
***********************************************************************/

{
  if (x[2] == ZERO)
    return -1;

  switch (k)
  {
    case 0:  *f = x[0] - x[1] + SIN(x[0]) + COS(x[1]) - PI / EIGHT;
             break;
    case 1:  *f = sqr(x[0]) - FOUR * sqr(x[1]) - LOG(sqr(x[2]));
             break;
    case 2:  *f = x[2] - ONE + TWO * x[0] - FOUR * x[1];
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt1_text(void)

{
  return
    "x0 - x1 + sin(x0) + cos(x1) - pi / 8  =  0\n"
    "x0^2 - 4 * x1^2 - ln(x2^2)            =  0\n"
    "x2 - 1 + 2 * x0 - 4 * x1              =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt2(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 2                               *
***********************************************************************/

{
  if (x[0] <= ZERO)
    return -1;

  switch (k)
  {
    case 0:  *f = sqr(x[0]) + x[0] + sqr(x[1]) - TWO;
             break;
    case 1:  *f = sqr(x[0]) + x[1] - sqr(x[1]) - ONE + LOG(x[0]);
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt2_text(void)

{
  return
    "x0^2 + x0 + x1^2 - 2           =  0\n"
    "x0^2 + x1 - x1^2 - 1 + ln(x0)  =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt3(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 3                               *
***********************************************************************/

{
  if (x[0] < MACH_EPS)
    return -1;
  if (x[1] + x[2] < MACH_EPS)
    return -1;

  switch (k)
  {
    case 0:  *f = TEN * x[0] + x[1] + x[2] + x[3] - (REAL)20.0 +
                  sqr(SIN(x[0])) + sqr(COS(x[1]));
             break;
    case 1:  *f = x[0] + (REAL)20.0 * x[1] + x[2] + x[3] -
                  (REAL)48.0 + ONE / POW(x[0], SIX);
             break;
    case 2:  *f = sqr(x[0] + x[1]) + (REAL)30.0 * x[2] + x[3] -
                  (REAL)97.0 + LOG(x[0]) + LOG(x[1] + x[2]);
             break;
    case 3:  *f = x[0] + x[1] + x[2] + (REAL)40.0 * x[3] -
                  (REAL)166.0 + sqr(x[0]);
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt3_text(void)

{
  return
    "10 * x0 + x1 + x2 + x3 - 20 + sin(x0)^2 + cos(x1)^2     =  0\n"
    "x0 + 20 * x1 + x2 + x3 - 48 + 1 / x0^6                  =  0\n"
    "(x0 + x1)^2 + 30 * x2 + x3 - 97 + ln(x0) + ln(x1 + x2)  =  0\n"
    "x0 + x1 + x2 + 40 * x3 - 166 + x0^2                     =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt4(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 4                               *
***********************************************************************/

{
  switch (k)
  {
    case 0:  *f = x[0] - EXP(SIN(x[1]));
             break;
    case 1:  *f = PI * EXP(SIN(PI * x[0])) + x[1];
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt4_text(void)

{
  return
    "x0 - exp(sin(x1))            =  0\n"
    "pi * exp(sin(pi * x0)) + x1  =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt5(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 5                               *
***********************************************************************/

{
  switch (k)
  {
    case 0:  *f = x[0] * x[0] * x[0] + ONE;
             break;
    case 1:  *f = x[1] * x[1] * x[1] + ONE;
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt5_text(void)

{
  return
    "x0^3 + 1  =  0\n"
    "x1^3 + 1  =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt6(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 6                               *
***********************************************************************/

{
  switch (k)
  {
    case 0:  *f = SIN(x[0]) - x[1];
             break;
    case 1:  *f = x[0] - COS(x[1]);
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt6_text(void)

{
  return
    "sin(x0) - x1  =  0\n"
    "x0 - cos(x1)  =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt7(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 7                               *
***********************************************************************/

{
  int  i;
  REAL fi;

  if (k < 0 || k >= beispiel[7].n)     /* invalid function number?    */
    return -1;

  for (fi = ZERO, i = 0; i < beispiel[7].n; i++)
    fi += x[i];

  *f = fi - x[k] - k - ONE + sqr(x[k]);

  return 0;
}


static char *fkt7_text(void)

{
  return
    "(x0 + ... +x(n-1)) - x0     - 0     - 1 + x0^2      =  0\n"
    "(x0 + ... +x(n-1)) - x1     - 1     - 1 + x1^2      =  0\n"
    "(x0 + ... +x(n-1)) - x2     - 2     - 1 + x2^2      =  0\n"
    "...                                                 =  ...\n"
    "(x0 + ... +x(n-1)) - x(n-1) - (n-1) - 1 + x(n-1)^2  =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt8(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 8                               *
***********************************************************************/

{
  switch (k)
  {
    case 0:  *f = sqr(x[0]) + sqr(x[1]) - FIVE;
             break;
    case 1:  *f = x[0] + x[0] - x[1];
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt8_text(void)

{
  return
    "x0^2 + x1^2 - 5  =  0\n"
    "2 * x0 - x1      =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt9(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 9                               *
***********************************************************************/

{
  switch (k)
  {
    case 0:  *f = HALF * sqr(x[0] - ONE) + (REAL)0.25 * sqr(x[1]) - ONE;
             break;
    case 1:  *f = x[1] - SIN(x[0]);
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt9_text(void)

{
  return
    "(x0 - 1)^2 / 2 + x1^2 / 4 - 1  =  0\n"
    "x1 - sin(x0)                   =  0\n";
}



/* ------------------------------------------------------------------ */

static int fkt10(int k, REAL x[], REAL *f)

/***********************************************************************
*                         Test example 10                              *
***********************************************************************/

{
  switch (k)
  {
    case 0:  *f = EXP(sqr(x[0]) + sqr(x[1])) - THREE;
             break;
    case 1:  *f = x[0] + x[1] - SIN(THREE * (x[0] + x[1]));
             break;
    default: return -1;
  }

  return 0;
}


static char *fkt10_text(void)

{
  return
    "exp(x0^2 + x1^2) - 3          =  0\n"
    "x0 + x1 - sin(3 * (x0 + x1))  =  0\n";
}



/* ------------------------------------------------------------------ */

bsptyp1 *nlgls_waehlen(int nummer)

/***********************************************************************
* Using this function, one can select one of the examples in the vector*
* beispiel for testing.                                                *
* If the value of nummer lies inside the valid range of indices for    *
* beispiel, the address of the corresponding vector element is produced*
* as the return value, otherwise the zero return indicates an error.   *
***********************************************************************/

{
  if (nummer < 0 || nummer >= sizeof(beispiel) / sizeof(*beispiel))
    return NULL;                    /* invalid number for the example */

  return &beispiel[nummer];
}
