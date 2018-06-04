#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
*                                                                      *
* Test examples for solving two point boundary value problems of first *
* order via the shooting method.                                       *
*                                                                      *
* Each example specifies the right hand side of a first order ordinary *
* system of DEs, followed by an alpha-numeric description of the       *
* example, the boundary condition and an alpha-numeric description of  *
* it.                                                                  *
*                                                                      *
* Using the function rwp_waehlen() the user may choose amongst examples*
* registered in the vector 'beispiel'.                                 *
*                                                                      *
***********************************************************************/

#include <basis.h>    /*  for  SIN, EXP, COS, COSH, NULL, REAL, HALF, */
                      /*       ONE, TWO, FOUR, NINE, THREE            */
#include <t_rwp.h>    /*  for  bsptyp, rwp_waehlen                    */



/* ---------------------- Boundary value problem 0 ------------------ */

/***********************************************************************
* Compute the value f of the right hand side of an explicit DE system  *
* at (x,y).                                                            *
***********************************************************************/
static void dgl0(REAL x, REAL *y, REAL *f)
{
  f[0] = y[0] * y[1] + COS(x) - HALF * SIN(TWO * x);
  f[1] = y[0] * y[0] + y[1] * y[1] - (ONE + SIN(x));
}

/***********************************************************************
* alpha-numeric description of above function                          *
***********************************************************************/
static char *dgltxt0(void)
{
  return
    "y1' = y1 * y2 + cos(x) - 0.5 * sin(2.0*x)\n"
    "y2' = y1 * y1 + y2 * y2 - (1 + sin(x))\n";
}

/***********************************************************************
* Compute value r for the boundary condition of the two point problem  *
* at the boundary points ya and yb.                                    *
***********************************************************************/
static void rand0(REAL *ya, REAL *yb, REAL *r)
{
  r[0] = ya[0] + ya[1] - ONE;
  r[1] = yb[1] - yb[0] - ONE;
}

/***********************************************************************
* alpha-numeric description of the boundary condition                  *
***********************************************************************/
static char *randtxt0(void)
{
  return
    "y1(a) + y2(a) - 1 = 0\n"
    "y2(b) - y1(b) - 1 = 0\n";
}



/* ---------------------- Boundary value problem 1 ------------------ */

/***********************************************************************
* Compute the value f of the right hand side of an explicit DE system  *
* at (x,y).                                                            *
***********************************************************************/
static void dgl1(REAL x, REAL *y, REAL *f)
{
  f[0] = y[1];
  f[1] = (ONE + y[1] * y[1]) / y[0];
}

/***********************************************************************
* alpha-numeric description of above function                          *
***********************************************************************/
static char *dgltxt1(void)
{
  return
    "y1' = y2\n"
    "y2' = (1 + y2 * y2) / y1\n";
}

/***********************************************************************
* Compute value r for the boundary condition of the two point problem  *
* at the boundary points ya and yb.                                    *
***********************************************************************/
static void rand1(REAL *ya, REAL *yb, REAL *r)
{
  r[0] = ya[0] - COSH(-ONE);
  r[1] = yb[0] - COSH(ONE);
}

/***********************************************************************
* alpha-numeric description of the boundary condition                  *
***********************************************************************/
static char *randtxt1(void)
{
  return
    "y1(a) - cosh(-1) = 0\n"
    "y1(b) - cosh(1)  = 0\n";
}



/* ---------------------- Boundary value problem 2 ------------------ */

/***********************************************************************
* Compute the value f of the right hand side of an explicit DE system  *
* at (x,y).                                                            *
***********************************************************************/
static void dgl2(REAL x, REAL *y, REAL *f)
{
  f[0] = y[1];
  f[1] = -FOUR * y[0] + EXP(x);
}

/***********************************************************************
* alpha-numeric description of above function                          *
***********************************************************************/
static char *dgltxt2(void)
{
  return
    "y1' = y2\n"
    "y2' = -4 * y1 + exp(x)\n";
}

/***********************************************************************
* Compute value r for the boundary condition of the two point problem  *
* at the boundary points ya and yb.                                    *
***********************************************************************/
static void rand2(REAL *ya, REAL *yb, REAL *r)
{
  r[0] = ya[0];
  r[1] = yb[1];
}

/***********************************************************************
* alpha-numeric description of the boundary condition                  *
***********************************************************************/
static char *randtxt2(void)
{
  return
    "y1(a) = 0\n"
    "y2(b) = 0\n";
}



/* ---------------------- Boundary value problem 3 ------------------ */

/***********************************************************************
* Compute the value f of the right hand side of an explicit DE system  *
* at (x,y).                                                            *
***********************************************************************/
static void dgl3(REAL x, REAL *y, REAL *f)
{
  f[0] = y[1];
  f[1] = y[2];
  f[2] = y[3];
  f[3] = y[4];
  f[4] = ((REAL)45.0 * y[2] * y[3] * y[4] -
          (REAL)40.0 * y[3] * y[3] * y[3]) / (NINE * y[2] * y[2]);
}

/***********************************************************************
* alpha-numeric description of above function                          *
***********************************************************************/
static char *dgltxt3(void)
{
  return
    "y1' = y2\n"
    "y2' = y3\n"
    "y3' = y4\n"
    "y4' = y5\n"
    "y5' = (45 * y3 * y4 * y5 - 40 * y4 * y4 * y4) / (9 * y3 * y3)\n";
}

/***********************************************************************
* Compute value r for the boundary condition of the two point problem  *
* at the boundary points ya and yb.                                    *
***********************************************************************/
static void rand3(REAL *ya, REAL *yb, REAL *r)
{
  r[0] = ya[0] - ONE;
  r[1] = yb[0] * (ONE - TWO * ya[1]);
  r[2] = yb[0] * yb[1] - HALF;
  r[3] = ONE / THREE * ya[2] + HALF * yb[1];
  r[4] = ONE / THREE * ya[3] - (REAL)0.25 * yb[0] * yb[1];
}

/***********************************************************************
* alpha-numeric description of the boundary condition                  *
***********************************************************************/
static char *randtxt3(void)
{
  return
    "y1(a) - 1 = 0\n"
    "y1(b) * (1 - 2 * y2(a))       = 0\n"
    "y1(b) * y2(b) - 0.5           = 0\n"
    "y3(a) / 3 + y2(b) / 2         = 0\n"
    "y4(a) / 3 - y1(b) * y2(b) / 4 = 0\n";
}



/* ---------------------- Boundary value problem 4 ------------------ */

/***********************************************************************
* Compute the value f of the right hand side of an explicit DE system  *
* at (x,y).                                                            *
***********************************************************************/
static void dgl4(REAL x, REAL *y, REAL *f)
{
  f[0] = y[1];
#if defined(BC3)
  _fpreset();
#endif
/*#ifndef BC3*/
  f[1] = -y[0] * COSH(x);
/*#else*/
  /* COSH() wird an dieser Stelle deshalb nicht verwendet, weil  */
  /* coshl() aus der Bibliothek von Borland C++ 3.0 einen Fehler */
  /* zu enthalten scheint: Mit den Eingabewerten  x = 0.1  und   */
  /* y = (0.177696409473669979,1.77696409473664004)  in dgl4()   */
  /* liefert coshl() folgenden Laufzeitfehler:                   */
  /* "Floating point error: stack fault."                        */
  /* Wenn man aber mit _fpreset() die Gleitkommabibliothek von   */
  /* BC3 reinitialisiert, klappt alles auch mit coshl(). Daher   */
  /* wurde die Ersatzfunktion hier vorlaeufig auskommentiert     */
  /* (bis zum naechsten Fehler...).                              */
/*  f[1] = -y[0] * (EXP(x) + EXP(-x)) / TWO;*/
/*#endif*/
}

/***********************************************************************
* alpha-numeric description of above function                          *
***********************************************************************/
static char *dgltxt4(void)
{
  return
    "y1' = y2\n"
    "y2' = -y1 * cosh(x)\n";
}

/***********************************************************************
* Compute value r for the boundary condition of the two point problem  *
* at the boundary points ya and yb.                                    *
***********************************************************************/
static void rand4(REAL *ya, REAL *yb, REAL *r)
{
  r[0] = ya[0];
  r[1] = yb[0] - ONE;
}

/***********************************************************************
* alpha-numeric description of the boundary condition                  *
***********************************************************************/
static char *randtxt4(void)
{
  return
    "y1(a)     = 0\n"
    "y1(b) - 1 = 0\n";
}



/* ---------------------- Boundary value problem 5 ------------------ */
/* ---------------------- (das aus dem FNUM-Test) ------------------- */

/***********************************************************************
* Compute the value f of the right hand side of an explicit DE system  *
* at (x,y).                                                            *
***********************************************************************/
static void dgl5(REAL x, REAL *y, REAL *f)
{
  f[0] = y[1];
  f[1] = -y[0] * y[0] * y[0];
}

/***********************************************************************
* alpha-numeric description of above function                          *
***********************************************************************/
static char *dgltxt5(void)
{
  return
    "y1' = y2\n"
    "y2' = -y1^3\n";
}

/***********************************************************************
* Compute value r for the boundary condition of the two point problem  *
* at the boundary points ya and yb.                                    *
***********************************************************************/
static void rand5(REAL *ya, REAL *yb, REAL *r)
{
  r[0] = ya[0];
  r[1] = yb[0];
}

/***********************************************************************
* alpha-numeric description of the boundary condition                  *
***********************************************************************/
static char *randtxt5(void)
{
  return
    "y1(a) = 0\n"
    "y1(b) = 0\n";
}



/* ------------------------------------------------------------------ */

static bsptyp beispiel[] =                /* Vector that registers the*/
  {{ 2, dgl0, dgltxt0, rand0, randtxt0 }, /* above examples           */
   { 2, dgl1, dgltxt1, rand1, randtxt1 },
   { 2, dgl2, dgltxt2, rand2, randtxt2 },
   { 5, dgl3, dgltxt3, rand3, randtxt3 },
   { 2, dgl4, dgltxt4, rand4, randtxt4 },
   { 2, dgl5, dgltxt5, rand5, randtxt5 }
  };



/***********************************************************************
* The below function allows to choose one of the test examples         *
* registered in the vector `beispiel`, provided the chosen number is   *
* valid.                                                               *
* if 'nummer' lies inside the valid indices of 'beispiel', the address *
* of the coresponding vector entry is returned, otherwise the error    *
* is reported.                                                         *
***********************************************************************/

bsptyp *rwp_waehlen(int nummer)

{
  if (nummer < 0 || nummer >= sizeof(beispiel) / sizeof(*beispiel))
    return NULL;                    /* invalid number for the example */

  return &beispiel[nummer];
}
