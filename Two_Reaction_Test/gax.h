/*.BA*/
/*.KA{C 15}{Numerical Integration}{Numerical Integration}*/
/*.BE*/
/* ----------------------- DECLARATIONS gax.h ----------------------- */

/***********************************************************************
* GAX.h:  Include file for GAX                         Egg, 02.14.1991 *
***********************************************************************/
#ifndef GAX_H_INCLUDED
#define GAX_H_INCLUDED

typedef enum
{
  GxTRAPEZ,
  GxGAUSSJORDAN,
  GxCLENSHAWCURTIS,
  GxROMBERG,
  GxNEWTONCOTES
} temptyp;

/*
#define GxTRAPEZ           0
#define GxGAUSSJORDAN      1
#define GxCLENSHAWCURTIS   2
#define GxROMBERG          3
#define GxNEWTONCOTES      4
*/

typedef  struct            /* interval information */
{
  REAL FSch;       /* error estimate               */
  REAL QRes;       /* quadrature result            */
  REAL Links;      /* left hand interval boundary  */
  REAL Rechts;     /* right hand interval boundary */
} IntInfo;

extern REAL* LegCoef;          /* for array of coeffizients           */
extern int   LegGrad;          /* for length of array of coefficients */

int  ClenCurtStGew  (int n, REAL* Null, REAL* Gew);
int  ClenCurt       (REAL f(REAL), int m, REAL* t, int n,
                     REAL* Ak, REAL* Tk, REAL* Resultat);
int  LegendreCoeff  (REAL* s, int Grad, REAL *Leg_r);
REAL LegPolWert     (REAL x, REAL *LegCoef, int LegGrad);
int  LegendreNullst (REAL* Root, int Grad, REAL *LegCoef, REAL *Leg_r);
int  AdaQuaStGew    (int Grad, REAL* xNull, int Flag, REAL* zGew);
int  OrtogP         (int n, REAL* Integrale,
                     REAL* StStelle, REAL* Gewicht);
int  QuaRom         (REAL a, REAL b, REAL eps, int* n,
                     REAL* h, REAL f(REAL), REAL* Qwert, REAL* Sch);
int  QuaNeC         (REAL von, REAL bis, int AnzInt, int nrV,
                     REAL f(REAL), REAL* Qwert, int* Ordnung,
                     REAL* Schaetz);

int  GaxT           (REAL*Zerlegung, int N, REAL RelError, REAL f(REAL),
                     int* nMax, int Methode, int Anz, REAL* QuadResult,
                     REAL* AbsFehl, IntInfo*** Intervalle,
                     void *vmblock);

void GaxF           (void *vmblock);

int  Gax            (REAL*Zerlegung, int N, REAL RelError, REAL f(REAL),
                     int* nMax, int Methode, int Anz, REAL* QuadResult,
                     REAL* AbsFehl);
#endif

/* --------------------------- END gax.h ---------------------------- */
