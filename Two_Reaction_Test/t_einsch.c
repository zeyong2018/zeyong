#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
*                                                                      *
* Examples for right hand sides of first order ordinary differential   *
* equations; in alpha-numeric form and the theoretical solution if     *
* known.                                                               *
* When running the main test program the function dgl_waehlen() will   *
* chose among the examples.                                            *
*                                                                      *
***********************************************************************/

#include <basis.h>                      /*  for  LOG, REAL, ONE, FABS */
#include <t_einsch.h>                   /*  for  bsptyp, dgl_waehlen  */



/* right hand side f of the differential equation y' = f(x, y) */

static REAL y0_strich(REAL x, REAL y)

{
  return y / x + (ONE / LOG(x));
}

/*exact solution of above DE depending on initial condition y(x0) = y0*/

static REAL y0_exakt(REAL x0, REAL y0, REAL x)

{
  return x * (y0 / x0 + (LOG(FABS(LOG(x))) - LOG(FABS(LOG(x0)))));
}



/* right hand side of above DE in alpha-mumeric form */

static char *dgl0_text(void)

{
  return "y / x + (1.0 / log(x))";
}

static bsptyp beispiel[] =           /* Vector with registers for all */
  {{ y0_strich, dgl0_text, y0_exakt }/* DEs defined above             */
  };                                 /* zero pointer in case exact    */
                                     /* solution is not known         */


/***********************************************************************
* This function allows to pick one of the differential equations in the*
* vector beispiel for test purposes.                                   *
* If the parameter  nummer lies in the valid range for  beispiel, the  *
* address of the corresponding vector entry is returned, otherwise the *
* zero pointer indicates the error.                                    *
***********************************************************************/

bsptyp *dgl_waehlen(int nummer)

{
  if (nummer < 0 || nummer >= sizeof(beispiel) / sizeof(*beispiel))
    return NULL;                        /* invalid number for example */

  return &beispiel[nummer];
}
