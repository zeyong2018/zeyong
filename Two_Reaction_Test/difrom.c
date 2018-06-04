#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE difrom.c ------------------------- */

#include <basis.h>
#include <vmblock.h>
#include <stdio.h>
#include  <math.h>
/*.BA*/

int difrom (REAL  func (REAL),
/*.IX{difrom}*/
            REAL  x0,
            REAL  eps,
            int   n,
            REAL  h,
            REAL* res,
            REAL* er_app,
            int*  nend,
            REAL* hend)
/***********************************************************************
*  Computes an approximation for the derivative of func at x0 using the*
*  ROMBERG method.                                                     *
.BE*)
*                                                                      *
*  Input parameters :                                                  *
*                                                                      *
*    REAL   func (REAL)   Name of function to be differentiated        *
*    REAL   x0            place where derivative is to be found        *
*    REAL   eps           desired accuracy                             *
*    int    n             max. number of columns in the Romberg scheme *
*                         (n > 1)                                      *
*    REAL   h             initial step size                            *
*                                                                      *
*  Output parameters:                                                  *
*                                                                      *
*    REAL   *res          approximate derivative                       *
*    REAL   *er_app       error estimate for  res                      *
*    int    *nend         number of columns actually used in scheme    *
*    REAL   *hend         final step size                              *
*                                                                      *
*  Return value :                                                      *
*                                                                      *
*    0:   no error: er_app < eps                                       *
*    1:   n < 1  or  eps <= 0  or  h < MACH_EPS                        *
*    2:   desired accuracy not reached after n steps                   *
*    3:   step size drooped below  MACH_EPS                            *
*    4:   lack of sufficient memory                                    *
*                                                                      *
*  constant used  :     MACH_EPS                                       *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  register i, j;
  int      m, error;
  REAL     h2, d1, d2, *d;
  void     *vmblock;

  if (n <= 1 || eps <= 0.0 || h < MACH_EPS)            /* check input */
    return (1);
                                               /* allocate aux vector */
  vmblock = vminit();
  d = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
    return (4);

  h2 = 2.0 * h;
  d [0] = (func (x0 + h) - func (x0 - h)) / h2;

/***********************************************************************
* This loop runs until the maximum of Romberg rows is filled or until  *
* the step size use drops below the machine constant or if the desired *
* accuracy is reached.                                                 *
***********************************************************************/

  for (error = 2, j = 1; j <= n - 1; j++)
  {
    d [j] = 0.;
    d1    = d [0];
    h2    = h;
    h    *= 0.5;
    *nend = j;

    if (h < MACH_EPS)
    {
      error = 3;                              /* step size less than  */
      break;                                  /* machine constant     */
    }

    d [0] = (func (x0 + h) - func (x0 - h)) / h2;

    for (m = 4, i = 1; i <= j; i++, m *= 4)
    {
      d2 = d [i];
      d [i] = (m * d [i-1] - d1) / (m-1);
      d1 = d2;
    }

    *er_app = FABS (d [j] - d [j-1]);

    if (*er_app < eps)
    {
      error = 0;                          /* desired accuracy reached */
      break;
    }
  }

  *res  = d [*nend];                   /* save final values           */
  *hend = h;

  vmfree(vmblock);                     /* free storage                */

  return (error);
}

/* -------------------------- END difrom.c -------------------------- */
