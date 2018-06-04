#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE einschr.c ------------------------ */

/***********************************************************************
*                                                                      *
* Solve a differential equation via one-step methods                   *
* --------------------------------------------------                   *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Jobst Hoffmann                                 *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 3.10.1992                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>         /*  for  min, FALSE, TRUE, FABS, MACH_EPS, */
                           /*       dglfnk, REAL, ZERO, ONE, TWO,     */
                           /*       THREE, SIX, EIGHT                 */
#include <einschr.h>       /*  for  dglesv                            */

/* ------------------------------------------------------------------ */

static REAL euler_cauchy
/*.IX{euler\unt cauchy}*/
                        (
                         REAL   x0,
                         REAL   y0,
                         REAL   h,
                         dglfnk dgl,
                         int    neu_f0
                        )

/***********************************************************************
* Perform one Euler-Cauchy step                                        *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0:     x-value                                                      *
* y0      y-value                                                      *
* h       step size                                                    *
* dgl     pointer to the function that evaluates the right hand side   *
*         of the differential equation                                 *
* neu_f0: Flag :                                                       *
*         != 0 : The slope at  (x0,y0) must be computed                *
*          = 0 : The slope at  (x0,y0) is known from a previous call   *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* Approximate solution of the differential equation at x0 + h from the *
* Euler-Cauchy polygonal method.                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, dglfnk                                                         *
***********************************************************************/

{
  static REAL f0;

  if (neu_f0)
    f0 = (*dgl)(x0, y0);

  return y0 + h * f0;
}

/* ------------------------------------------------------------------ */

static REAL heun
/*.IX{heun}*/
                (
                 REAL   x0,
                 REAL   y0,
                 REAL   h,
                 dglfnk dgl,
                 int    neu_f0
                )

/***********************************************************************
* Perform one Heun step of integration                                 *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0:     x-value                                                      *
* y0      y-value                                                      *
* h       step size                                                    *
* dgl     pointer to the function that evaluates the right hand side   *
*         of the differential equation                                 *
* neu_f0: Flag :                                                       *
*         != 0 : The slope at  (x0,y0) must be computed                *
*          = 0 : The slope at  (x0,y0) is known from a previous call   *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* Approximate solution of the differential equation at x0 + h from the *
* Heun method.                                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, dglfnk                                                         *
***********************************************************************/

{
  static REAL f0;
  REAL        hilf1,
              f1;

  if (neu_f0)
    f0 = (*dgl)(x0, y0);
  hilf1 = y0 + h * f0;
  f1 = (*dgl)(x0 + h, hilf1);

  return y0 + HALF * h * (f0 + f1);
}

/* ------------------------------------------------------------------ */

static REAL runge_kutta
/*.IX{runge\unt kutta}*/
                       (
                        REAL   x0,
                        REAL   y0,
                        REAL   h,
                        dglfnk dgl,
                        int    neu_f0
                       )

/***********************************************************************
* Perform one Runge-Kutta step                                         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0:     x-value                                                      *
* y0      y-value                                                      *
* h       step size                                                    *
* dgl     pointer to the function that evaluates the right hand side   *
*         of the differential equation                                 *
* neu_f0: Flag :                                                       *
*         != 0 : The slope at  (x0,y0) must be computed                *
*          = 0 : The slope at  (x0,y0) is known from a previous call   *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* Approximate solution of the differential equation at x0 + h from the *
* Runge-Kutta method.                                                  *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, dglfnk, TWO, SIX                                               *
***********************************************************************/

{
  static REAL f0;
  REAL        k1, k2, k3, k4;

  if (neu_f0)
    f0 = (*dgl)(x0, y0);
  k1 = h * f0;
  k2 = h * (*dgl)(x0 + HALF * h, y0 + HALF * k1);
  k3 = h * (*dgl)(x0 + HALF * h, y0 + HALF * k2);
  k4 = h * (*dgl)(x0 + h, y0 + k3);
  return y0 + (k1 + TWO * (k2 + k3) + k4) / SIX;
}

/* ------------------------------------------------------------------ */

typedef struct { REAL x, y; } puffertyp[5];        /* define type of  */
/*.IX{puffertyp}*/
                                                   /* circular buffer */



/* ------------------------------------------------------------------ */

static REAL interpolieren
/*.IX{interpolieren}*/
                         (
                          puffertyp stuetzpunkt,
                          int       anzahl,
                          REAL      x)

/***********************************************************************
* Use Newton interpolation to find the value of the interpolating      *
* polynomial through anzahl many nodes at x.                           *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* stuetzpunkt: Vector of nodes                                         *
* anzahl:      number of nodes (maximally 5)                           *
* x:           value used for evaluation                               *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* Functional value of interpolating polynomial at x                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* puffertyp, REAL                                                      *
***********************************************************************/

{
  REAL   delta[5],  /* Vector of divided differences                  */
         y;         /* Function value of interpolating polynomial at x*/
  int    i, j;      /* running variables                              */


  for (i = 0; i < anzahl; i++)                /* store nodes in delta */
    delta[i] = stuetzpunkt[i].y;

  for (i = 1; i < anzahl; i++)         /* compute divided differences */
    for (j = anzahl - 1; j >= i; j--)       /* and thus the polynomial*/
      delta[j] =  delta[j - 1] - delta[j],  /* coefficients           */
      delta[j] /= stuetzpunkt[j - i].x -
                  stuetzpunkt[j].x;

  for (y = delta[anzahl - 1], i = anzahl - 2; i >= 0; i--)
    y *= x - stuetzpunkt[i].x,
    y += delta[i];

  return y;
}

/* ------------------------------------------------------------------ */

#define F0_NEU                  1  /* compute slope at (x,y)          */
#define F0_ALT                  0  /* slope at (x,y) is known         */

#define MAXIMALE_AUFRUFZAHL  5000   /* maximal number of evaluations  */
                                    /* right hand side                */

typedef REAL (*verfahrenstyp)(REAL   x0,     /* Type of function used */
/*.IX{verfahrenstyp}*/
                              REAL   y0,     /* for integration       */
                              REAL   h,
                              dglfnk dgl,
                              int    neu_f0
                             );


/* ------------------------------------------------------------------ */
/*.BA*/

int dglesv         /* One-step methods for 1st order DEs .............*/
/*.IX{dglesv}*/
          (
           REAL   *x,        /* initial/final x-value ................*/
           REAL   *y,        /* ditto for y ..........................*/
           dglfnk dgl,       /* righ thand side of DE ................*/
           REAL   xend,      /* desired final x-value ................*/
           REAL   *h,        /* initial/final step size ..............*/
           REAL   epsabs,    /* absolute error bound .................*/
           REAL   epsrel,    /* relative error bound .................*/
           int    intpol,    /* use interpolation ? ..................*/
           int    methode,   /* select method  (0, 1, 2) .............*/
           int    rand,      /* xend not reached ? ...................*/
           int    neu        /* pass on old data ? ...................*/
          )                  /* error code ...........................*/

/***********************************************************************
* Compute the solution of the first order ordinary differential        *
* equation                                                             *
*                                                                      *
*     y' = f(x,y)      with initial value      y(x0) = y0              *
*                                                                      *
* at xend. One may choose among three methods :                        *
*                                                                      *
*   1. Euler-Cauchy polygonal method,                                  *
*   2. method of Heun and                                              *
*   3. explicit Runge-Kutta method.                                    *
.BE*)
*                                                                      *
* Here one can                                                         *
*                                                                      *
*   a) use values from the last call of this function, or              *
*   b) start anew.                                                     *
*                                                                      *
* Re a) Here one can choose to find the value of the solution at xend  *
*       by interpolation or can continue the integration at the final  *
*       x-value that was reached with the last call.                   *
*                                                                      *
* The interpolation data is arranged in a cicular buffer. This buffer  *
* must be filled with data before an interpolation of the proper       *
* order can be performed. This might entail that in the first step     *
* we compute data way beyond the desired end point.                    *
*                                                                      *
* This program uses automatic step size control by either halving,     *
* or doubling the current step size depending on how the error bounds  *
* are met.                                                             *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x:        initial x-value                                            *
* y:        initial y-value at x                                       *
* dgl:      pointer to function that evaluates the right hand side     *
* xend:     final x-value for which a value is desired for the         *
*           solution of the DE; x < xend                               *
* h:        starting step size                                         *
* epsabs:\  error bounds, nonnegative. We perform the mixed test :     *
* epsrel:/      |local error|  <=  |Y| * epsrel + epsabs;              *
*           if epsrel = 0 we test only the absolute error;             *
*           if epsabs = 0 we test the relative error only.             *
* intpol:   Flag (used only when  neu = 0) :                           *
*            = 0 : continue with integration in next call              *
*           != 0 : try interpolating to xend from the current data     *
* methode:  select desired method :                                    *
*           = 0:  Euler-Cauchy                                         *
*           = 1:  Heun                                                 *
*           = 2:  explicit Runge-Kutta                                 *
* rand:     indicates whether the end point xend may be passed :       *
*            = 0 : xend may be passed                                  *
*           != 0 : do not pass xend                                    *
* neu:      indicates whether old data should be saved for new call :  *
*            = 0 : save old data                                       *
*           != 0 : no old data                                         *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x: final x-value reached (usually x = xend)                          *
* y: final y-value at x                                                *
* h: final step size (may be saved for next call)                      *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: integration went to xend.                                         *
* 1: solution at xend was found using interpolation.                   *
* 2: too many function evaluations                                     *
* 3: step size became less than 8 * machine constant.                  *
*    Increase h  and error bounds in subsequent calls.                 *
* 4: epsabs < 0  or  epsrel < 0 or both equal to zero.                 *
* 5: xend = x.                                                         *
* 6: imvalid number for method                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* puffertyp, interpolieren, euler_cauchy, heun, runge_kutta, REAL,     *
* dglfnk, min, FABS, MACH_EPS, FALSE, TRUE, ZERO, ONE, TWO, THREE,     *
* EIGHT                                                                *
.BA*)
***********************************************************************/
/*.BE*/

{
  int
    pufferlaenge,    /* number of nodes for interpolation             */
    schaetzung_ok,   /* Flag for error estimate                       */
    aufrufzahl;      /* number of calls of   dgl()                    */
  REAL
    x0,              /* (x0,y0) is the initial value of the problem   */
    y0,
    yhilf,           /* y-value for one step of size h from (x0,y0)   */
    y1,              /* y-value for one step of size 2*h from (x0,y0) */
    y2,              /* y-value for one step, size h from (x0+h,yhilf)*/
    fehlerschaetzung,/* error estimate derived from  y1 and y2        */
    schaetzfaktor,   /* constant for adjusting error estimates        */
    hh;              /* aux variable for step size                    */
  verfahrenstyp
    mein_verfahren;  /* pointer to method being used                  */
  static verfahrenstyp      /* Vector of pointers to the functions    */
    verfahren[] =           /* that realize the three methods         */
      {euler_cauchy, heun,
       runge_kutta};
  static REAL         /* Vector with estimate factors for the methods */
    schtzfkt[] =
      { ONE, THREE, (REAL)15.0 };
  static puffertyp       /* circular buffer                           */
    puffer;
  static int
    puflaeng[] =         /* Vector with lengths for the circular      */
      {2, 3, 5},         /* buffer needed for each method             */
    pufferindex =        /* Index of last entry (x,y) in circular     */
      0,                 /* buffer                                    */
    puffer_leer =        /* Flag indicating buffer is empty           */
      TRUE,
    puffer_voll =        /* Flag indicating buffer is full            */
      FALSE,
    aufruf_neu[] =       /* Vector with number of function evaluations*/
      {2, 5, 11},        /* that one step with each method uses       */
    aufruf_add[] =       /* Vector with number of additional function */
      {1, 3, 7};         /* evaluations for each method, if step is   */
                         /* not successful                            */

  /* --------------- check input parameters  ------------------------ */

  if (epsabs < ZERO || epsrel < ZERO ||
      (epsabs == ZERO && epsrel == ZERO))
    return 4;
  if (xend == *x)
    return 5;
  if (methode < 0 || methode >= 3)
    return 6;


  if (neu)                     /* do not use old contents of buffer ? */
  {
    pufferindex = 0;                    /* clear buffer and fill anew */
    puffer_voll = FALSE;                /* starting with initial      */
    puffer[0].x = x0 = *x;              /* values (x,y)               */
    puffer[0].y = y0 = *y;
    puffer_leer = FALSE;
  }

  else                                   /* continue and use buffer ? */
  {
    if (puffer_leer)                  /* buffer empty (first call) ?  */
      puffer[0].x = x0 = *x,                  /* start filling buffer */
      puffer[0].y = y0 = *y,
      puffer_leer = FALSE;
    x0 = puffer[pufferindex].x,      /* take newest entry from buffer */
    y0 = puffer[pufferindex].y;

    if (intpol && puffer_voll &&    /* We interpolate if desired, the */
        ((x0 > xend &&              /* buffer is full and the desired */
          *h > ZERO) ||             /* final x-value xend lies within */
         (x0 < xend &&              /* the nodes of the buffer        */
          *h < ZERO)))
    {
      *x = xend;
      *y = interpolieren(puffer, puflaeng[methode], xend);
      return 1;
    }
  }

  /* --------------------- initialize variablen --------------------- */

  hh = FABS(*h);                        /* correct step size          */
  if (rand)
    hh = min(hh, FABS(xend - *x));
  *h   = (xend > *x) ? hh : -hh;

  mein_verfahren = verfahren[methode];              /* choose method  */
  pufferlaenge   = puflaeng[methode];
  schaetzfaktor  = schtzfkt[methode];
  aufrufzahl     = 0;


  for ( ; ; )                                      /* try to solve DE */
  {
    aufrufzahl += aufruf_neu[methode];
    if (aufrufzahl > MAXIMALE_AUFRUFZAHL)          /* too many calls  */
    {                                              /* of  dgl()?      */
      *x = x0;
      *y = y0;
      return 2;
    }

    if (FABS(*h) < MACH_EPS * EIGHT * FABS(*x))     /* step size too  */
    {                                               /* small ?        */
      *x = x0;
      *y = y0;
      return 3;
    }

    /* ------------------ perform one step -------------------------- */

    yhilf = (*mein_verfahren)(x0,      y0,    *h,       dgl, F0_NEU);
    y1    = (*mein_verfahren)(x0,      y0,    *h * TWO, dgl, F0_ALT);
    y2    = (*mein_verfahren)(x0 + *h, yhilf, *h,       dgl, F0_NEU);


    for (schaetzung_ok = FALSE; ! schaetzung_ok; )      /* estimate   */
    {                                                   /* the error  */
      fehlerschaetzung = FABS(y2 - y1) / schaetzfaktor;

      /* If the error estimate turns out to be too large, we repeat   */
      /* the step with half the step size.                            */

      if (fehlerschaetzung > epsabs + FABS(y2) * epsrel)
      {
        *h *= HALF;
        aufrufzahl += aufruf_add[methode];
        if (aufrufzahl > MAXIMALE_AUFRUFZAHL)      /* too many calls  */
        {                                          /* of dgl() ?      */
          *x = x0;
          *y = y0;
          return 2;
        }
        if (FABS(*h) < MACH_EPS * EIGHT * FABS(*x))   /* step size    */
        {                                             /* too small ?  */
          *x = x0;
          *y = y0;
          return 3;
        }
        y1    = yhilf;
        yhilf = (*mein_verfahren)(x0,      y0,    *h, dgl, F0_ALT);
        y2    = (*mein_verfahren)(x0 + *h, yhilf, *h, dgl, F0_NEU);
      }
      else
        schaetzung_ok = TRUE;
    }

    x0 += *h * TWO;                        /* compute better estimate */
    y0 = ((schaetzfaktor + ONE) * y2 - y1) / schaetzfaktor;

    if (pufferindex + 2 == pufferlaenge)   /* buffer full ?           */
      puffer_voll = TRUE;                  /* record this fact        */
    pufferindex++;                         /* store new acceptable    */
    pufferindex %= pufferlaenge;           /* approximation (x0,y0)   */
    puffer[pufferindex].x = x0;            /* in the buffer           */
    puffer[pufferindex].y = y0;

    /* - If both approximations  y1 and y2 do not differ much, the    */
    /* - step size is doubled for the next step                       */

    if (fehlerschaetzung < (REAL)0.0666 * (epsabs + FABS(y2) * epsrel))
      *h *= TWO;

    if (rand && xend != x0) /* passing right end point not allowed ?  */
    {
      *h = min(FABS(*h), HALF * FABS(xend - x0));
      if (xend <= x0)
        *h = -*h;
    }

    if (x0 == xend)               /* precisely at desired end point ? */
    {                                /* Done !                        */
      *x = xend;
      *y = y0;
      return 0;
    }

    if (puffer_voll &&                /* passed end point and buffer  */
        ((x0 > xend && *h > ZERO) ||  /* full ?                       */
         (x0 < xend && *h < ZERO)))
    {                                 /* compute solution at  xend    */
      *x = xend;                      /* using interpolation          */
      *y = interpolieren(puffer, pufferlaenge, xend);
      return 1;
    }

  }                                       /* End of big infinite loop */
}

/* -------------------------- END einschr.c ------------------------- */
