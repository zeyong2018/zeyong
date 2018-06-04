#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 15.10}{Romberg Integration}{Romberg Integration}*/

/*.BE*/
/* ------------------------- MODULE quarom.c ------------------------ */

#ifndef QuaRom_C
#define QuaRom_C

#include <basis.h>
#include <vmblock.h>
#include   <gax.h>
/*.BA*/

int QuaRom (
/*.IX{QuaRom}*/
            REAL   a,
            REAL   b,
            REAL   eps,
            int*  _m,
            REAL* _h,
            REAL   f(REAL),
            REAL* _Qwert,
            REAL* _Sch
           )
/***********************************************************************
* Evaluate the integral of  f(x) over the interval [a,b] using         *
* Romberg Integration.                                                 *
.BE*)
*                                                                      *
* Parameters :                                                         *
*   REAL   a, b       Interval of integration                          *
*   REAL   eps        accuracy bound for error estimate                *
*   int*   m          maximal number of columns and rows in Romberg    *
*                     Romberg scheme: m > 1                            *
*                     on output: actual number of Romberg columns      *
*   REAL*  h          initial step size  (b-a) / k, k integer          *
*                     on output: step size at the end                  *
*   REAL   f (REAL  ) Function                                         *
*   REAL   *Qwert     value for integral                               *
*   REAL   *Schaetz   error estimate for  Qwert                        *
*                                                                      *
* Return value :                                                       *
*   0:                o.k.                                             *
*   1:                accuracy not reached after m steps               *
*   2:                m improper                                       *
*   3:                eps improper                                     *
*   4:                lack of memory                                   *
*                                                                      *
* Author              Uli Eggermann, 2.2.1991                          *
.BA*)
***********************************************************************/
/*.BE*/
{
  #define  m     (*_m)                  /* serves the readability:    */
  #define  h     (*_h)                  /* m, h, Qwert and Sch are    */
  #define  Qwert (*_Qwert)              /* used as references         */
  #define  Sch   (*_Sch)                /* ("*" not needed)           */

  REAL   *L, L0 = ZERO, Lk = ZERO;
  int    pot4, N0, j, k;
  void   *vmblock;
/***********************************************************************
*        check input                                                   *
***********************************************************************/
  if (a == b) { Qwert = Sch = h = ZERO, m = 0; return (0); }
  if (m < 2)       return (2);
  if (eps <= ZERO) return (3);
/***********************************************************************
*        allocate aux arrays                                           *
***********************************************************************/
  vmblock = vminit();
  L = (REAL *)vmalloc(vmblock, VEKTOR, m + 1, 0);
  if (! vmcomplete(vmblock))
    return (4);
/***********************************************************************
*        adjust step size                                              *
***********************************************************************/
  h = min (FABS (h), FABS (b-a));          /* max length of interval  */
  if (b < a)  h = -h;                      /* direction  ?            */
  if (h == ZERO) h = b - a;                /* maximal step size       */
  N0 = (int) ((b - a) / h + 0.5);          /* round!                  */
/***********************************************************************
*        first Romberg row (one entry only)                            *
***********************************************************************/
  L[0] = .5 * (f (a) + f (b));
  for (k = 1; k < N0; k++) L[0] += f (a + k * h);
  L[0] *= h;
/***********************************************************************
*        lower Romberg rows                                            *
***********************************************************************/
  for (j = 1; j < m; j++)
  {
    h   /= 2;
    L0   = L[0];
    L[0] = L[j] = ZERO;
    for (k = 0; k < N0; k++)                    /* further quadrature */
      L[0] += f (a + (2 * k + 1) * h);
    L[0] = .5 * L0 + h * L[0];
    pot4 = 1;
    for (k = 1; k <= j; k++)
    {
      pot4 *= 4;
      Lk    = L[k];
      L[k]  = (pot4 * L[k-1] - L0) / (pot4 - 1);
      L0    = Lk;
    }
    Sch = L[j] - L[j-1];                                /* estimate   */
    if (FABS (Sch) < eps) break;              /* accuracy sufficient ?*/
    N0  *= 2;
  }
/***********************************************************************
*        Return                                                        *
***********************************************************************/
   Qwert = L[j];                   /* Quadrature value                */

   k = (j == m);                   /* for k: "accuracy not met"       */
   m = j;                          /* # of Romberg rows               */
   vmfree(vmblock);                /* free storage                    */

   return (k);                     /* return                          */

   #undef  m
   #undef  h
   #undef  Qwert
   #undef  Sch
}                                                           /* QuaRom */
#endif

/* -------------------------- END quarom.c -------------------------- */
