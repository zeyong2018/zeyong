#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* This program lists model problems for nonlinear least squares with   *
* real functions.                                                      *
* For each example we give a mathematical and an alpha-mumeric code.   *
*                                                                      *
* When running the governing program and using the function            *
* nlinansf_waehlen(), one can point to various examples of the type    *
* bsptyp3. The individual examples are declared as global for the      *
* module                                                               *
*                                                                      *
* The pointer PHI inside  bsptyp3 points to a C function after call    *
* of nlinansf_waehlen(), which computes the functional value for the   *
* least square function defined by the coefficients c[i] at the node x *
* and returns the functional value there.                              *
* PHI serves as model function for a nonlinear discrete approximation  *
* whose optimal coefficients c[i] are to be determined. PHI has to be  *
* suitably chosen by the user for the problem at hand.                 *
* ABL computes the partial derivatives of PHI with respect to its      *
* coefficients at x. The pointer PHI_text gives the address of the     *
* string of code that describes the model function alpha-numerically.  *
* The number of coeficients for the model function is given by n.      *
*                                                                      *
***********************************************************************/

#include <basis.h>                  /*  for  EXP, basis               */
#include <nlappphi.h>               /*  for  bsptyp3, linansf_waehlen */



/*--------------------------------------------------------------------*/

static REAL PHI0(REAL *c, REAL x)

/***********************************************************************
* compute the functional value of a real function given by the         *
* parameters c[i] at x and return the function value                   *
***********************************************************************/

{
                                      /*                    bx        */
  return c[0] * EXP(c[1] * x);        /* Example   y = a * e          */
                                      /*                              */
                                      /*           a = c[0], b = c[1] */
}



/*--------------------------------------------------------------------*/

static void ABL0(REAL x, REAL *c, REAL *f)

/***********************************************************************
* compute the partial derivatives f[i] of the model function PHI0 at x *
* with respect to its coefficients c[i]                                *
***********************************************************************/

{
                           /*        bx      dy    bx    dy        bx */
  f[0] = EXP(c[1] * x);    /* y = a e    =>  -- = e   ,  -- = a x e   */
                           /*                da          db           */
  f[1] = c[0] * x * f[0];  /*                                         */
                           /* a = c[0], b = c[1]                      */
}



/*--------------------------------------------------------------------*/

static char *PHI0_text(void)

/***********************************************************************
* alpha-numeric code of  PHI0                                          *
***********************************************************************/

{
  return
    "PHI(x) = c0 + exp(c1 * x)\n";
}



/*--------------------------------------------------------------------*/

static bsptyp3 beispiel[] =     /* Vector, which contains all model   */
  {{ 2, PHI0, ABL0, PHI0_text } /* function examples                  */
  };



/*--------------------------------------------------------------------*/

bsptyp3 *nliansf_waehlen(unsigned int nummer)

/***********************************************************************
* This function is used to select one example registered in the vector *
* beispiel for test purposes.                                          *
* If nummer lies in the valid domain of  beispiel, the address of the  *
* respective element of beispiel is returned; otherwise the zero       *
* address signifies an error.                                          *
***********************************************************************/

{
  if (nummer >= sizeof(beispiel) / sizeof(*beispiel))
    return NULL;                     /* invalid model function number */

  return &beispiel[nummer];
}
