#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 15.12}{Adaptive Quadrature Methods}
              {Adaptive Quadrature Methods}*/

/*.BE*/
/* --------------------------- MODUL gax.c -------------------------- */

#include <basis.h>
#include <vmblock.h>
#include <gax.h>          /* Prototypes and Constants for "Methode"   */

#ifdef __DJGPP__          /* DJGPP 2 with GNU CC 2.7.2 and LDOUBLE?   */
#if __DJGPP__ == 2        /* To prevent the compiler from crashing    */
#ifdef LDOUBLE            /* during compilation, in this special case */
#undef FABS               /* FABS has to be replaced by ABS.          */
#define FABS  ABS
#endif
#endif
#endif

#ifdef __EMX__            /* NT 0.9b v2 with GNU CC 2.7.2 and LDOUBLE?*/
#ifdef __WIN32__          /* To prevent the compiler from crashing    */
#ifdef LDOUBLE            /* during compilation, in this special case */
#undef FABS               /* FABS has to be replaced by ABS.          */
#define FABS ABS
#endif
#endif
#endif

/*.BA*/

int Gax (
/*.IX{Gax}*/
         REAL* Zerlegung,
         int   N,
         REAL  RelError,
         REAL  Fkt (REAL),
         int*  nMax,
         int   Methode,
         int   MethPar,
         REAL* QuadResult,
         REAL* AbsFehl
        )
/***********************************************************************
* Computes the integral of the function Fkt(x) over the interval       *
* (a,b) using the method specified in "Methode".                       *
*                                                                      *
* The interval (a,b) must be partitioned by the user suitably          *
*   a = partition [0], ..., partition [N] = b.                         *
*                                                                      *
* We compute the integral and local error estimate for each subinterval*
*                                                                      *
* Thereafter the subintervals with excessive errors are automatically  *
* halved and the integrals recomputed.                                 *
*                                                                      *
* The program Gax stops, when either the relative overall error        *
* has become less than RelError or when the maximally allowed number   *
* of subintervals (nMax) has been reached.                             *
.BE*)
*                                                                      *
* Parameters:                                                          *
*   REAL    Zerlegung []  partition of  (a,b)                          *
*   int     N             number of subintervals of  (a,b)             *
*   REAL    RelError      maximally allowed relative error             *
*   REAL    Fkt (REAL)    integrating function                         *
*   int     *nMax         maximal allowed number of subintervals       *
*   int     Methode       Number for formula used (see "GAX.H", below) *
*   int     MethPar       Parameter for different methods (see below)  *
*   REAL    *QuadResult   result of quadrature                         *
*   REAL    *AbsFehl      error for absolute error                     *
*   IntInfo ***Fields     Information for test (see GaxT)              *
*                                                                      *
* Method          | MethPar     | Meaning                              *
* ----------------|-------------|------------------------------------- *
* GxTRAPEZ        | -           | Trapezoidal rule                     *
* GxGAUSSJORDAN   | m>=2, m<=20 | m point Gauss formula                *
* GxCLENSHAWCURTIS| m>0, even   | Clenshaw-Curtis with m+1 weights     *
* GxROMBERG       | -           | Romberg method                       *
* GxNEWTONCOTES   | m>=2, m<=7  | Newton-Cotes with m nodes            *
*                                                                      *
* global variables that might be used depending on Methode:            *
*   StStellen [], Gewichte []                                          *
*                                                                      *
* Author: Uli Eggermann, 10.3.1991                                     *
.BA*)
***********************************************************************/
/*.BE*/
{
  IntInfo **T;                                /* allocate             */
  int  g;
  void *vmblock = NULL;
  vmblock = vminit();
  g = GaxT (Zerlegung, N, RelError,
            Fkt, nMax, Methode, MethPar,      /* central routine      */
            QuadResult, AbsFehl, &T, vmblock);
  GaxF (vmblock);                             /* NECESSARY after GaxT */
  return (g);
}

/***********************************************************************
* global variables :                                                   *
***********************************************************************/

REAL *StStellen;       /* Weights and nodes for                       */
REAL *Gewichte;        /*              GAUSSJORDAN and CLENSHAWCURTIS */

/***********************************************************************
* local  routines (Prototypes):                                        *
***********************************************************************/

void GxSort     (IntInfo** F, int N);
int  GxQuad     (IntInfo* Info, REAL f(REAL), int Methode, int Anz);
int  GxRiEp     (REAL  qWert, REAL* L, REAL* Leps, REAL* EpsExe,
                 REAL* EpsL, int M, REAL* Eps, REAL* AbsFehl,
                 REAL* qExe);
int  GxClenCurt (REAL f(REAL), int n, IntInfo* Info );
void GxGauss    (REAL f(REAL), int Anz, IntInfo* Info );
REAL GxxGauss   (REAL a, REAL b, int Anz, REAL f(REAL));

int GaxT (
/*.IX{GaxT}*/
          REAL*  Zerlegung,
          int    N,
          REAL   RelError,
          REAL   Fkt(REAL),
          int*  _FieldCount,
          int    Methode,
          int    MethPar,
          REAL* _QuadResult,
          REAL* _AbsFehl,
          IntInfo*** _Fields,
          void   *vmblock
         )
/***********************************************************************
* int GaxT (... , pointer)    (other parameters as in Gax, see above)  *
*                                                                      *
* GaxT performs the same task as Gax, except for one extra parameter:  *
* A pointer to an intervall information array pointer.                 *
* This supplies subintervals and can be used later for test purposes.  *
*                                                                      *
* The routine GaxF must be calld subsequently in order to free storage.*
***********************************************************************/
{
  #define Fields     (*_Fields)        /* Used as references          */
  #define FieldCount (*_FieldCount)
  #define QuadResult (*_QuadResult)    /* Redefined for better        */
  #define AbsFehl    (*_AbsFehl)       /* readability                 */

  IntInfo *CurrInt, *new1, *new2;
  int     i, Maximum = FieldCount, RetVal = 0, St_Gew;

  /*********************************************************************
  * Check input paramaters                                             *
  *********************************************************************/
  FieldCount = 0;
  if (Methode < GxTRAPEZ || Methode > GxNEWTONCOTES)
     return (1);                         /* correct quadrature method?*/
  if (N < 1)            return (2);      /* improper # of intervals   */
  if (RelError<1.e-12)  return (3);      /* rel. error too small ?    */
  if (N > Maximum)      return (4);      /* improper value of N       */

  switch (Methode)                       /* further checks            */
  {
    case GxGAUSSJORDAN:    if(MethPar<2||MethPar > 20) return(5); break;
    case GxCLENSHAWCURTIS: if(MethPar<2||MethPar%2!=0) return(6); break;
    case GxNEWTONCOTES:    if(MethPar<2||MethPar >  7) return(7); break;
  }

  /*********************************************************************
  * Allocate and initialize various aux arrays                         *
  *********************************************************************/
  Fields = (IntInfo **)vmalloc(vmblock, VVEKTOR,
                               Maximum, sizeof(*Fields));
  if (! vmcomplete(vmblock))
    return 8;

  St_Gew = (Methode == GxGAUSSJORDAN || Methode == GxCLENSHAWCURTIS);

  if (St_Gew)                          /* allocate storage as needed  */
  {
    StStellen = (REAL *)vmalloc(vmblock, VEKTOR, MethPar + 1, 0);
    Gewichte  = (REAL *)vmalloc(vmblock, VEKTOR, MethPar + 1, 0);
    if (! vmcomplete(vmblock))
      return 8;
  }
  switch (Methode)         /* compute nodes and weights as needed     */
  {
    case GxGAUSSJORDAN:
      if (AdaQuaStGew (MethPar, StStellen, 1, Gewichte) != 0)
         return 8;
      break;
    case GxCLENSHAWCURTIS:
      if (ClenCurtStGew (MethPar, StStellen, Gewichte) != 0)
         return 8;
      break;
  }

  /*********************************************************************
  * Compute intervals of all subintervals and sort by decreasing abs.  *
  *       error                                                        *
  *********************************************************************/
  QuadResult = AbsFehl = ZERO;

  for (i = 0; i < N; i++)
  {
    if ((Fields [FieldCount++] = CurrInt =
         (IntInfo *)vmalloc (vmblock, VVEKTOR,
          1, sizeof(**Fields))) == NULL)
      return 8;
    CurrInt->Links  = Zerlegung [i];
    CurrInt->Rechts = Zerlegung [i+1];
    if (GxQuad (CurrInt, Fkt, Methode, MethPar) != 0)
      return 9;
    QuadResult += CurrInt->QRes;
    AbsFehl    += FABS (CurrInt->FSch);
  }

  GxSort (Fields, FieldCount);              /* sort by error size     */

  if (AbsFehl == ZERO || QuadResult == ZERO)
    return 0;
  if (FABS(QuadResult) < AbsFehl)
    return 0;

  /*********************************************************************
  * Halve the interval with the largest error                          *
  *********************************************************************/
  while (AbsFehl > FABS (QuadResult * RelError)
         && FieldCount < Maximum)
  {
    new1 = Fields [0];                     /* new1 -> largest error   */
    if ((Fields[FieldCount++] = new2        /* new2 -> end of list    */
         = (IntInfo *)vmalloc (vmblock, VVEKTOR, 1, sizeof (**Fields))
       ) == NULL)
      return 8;
    QuadResult -= new1->QRes;            /* subtract partial integral */
    AbsFehl    -= FABS (new1->FSch);         /* ditto for error       */
    new2->Rechts = new1->Rechts;
    new2->Links  = new1->Rechts =                    /* new partition */
      0.5 * (new1->Links + new1->Rechts);

    if (GxQuad (new1, Fkt, Methode, MethPar) != 0)
      return 9;
    if (GxQuad (new2, Fkt, Methode, MethPar) != 0)
      return 9;

    QuadResult +=       new1->QRes  +       new2->QRes;        /* new */
    AbsFehl    += FABS (new1->FSch) + FABS (new2->FSch);

    GxSort (Fields, FieldCount);             /* sort by error size    */
  }

  if (FieldCount == Maximum && AbsFehl > FABS(QuadResult*RelError))
    RetVal = 10;
  return (RetVal);

  #undef FieldCount
  #undef Fields
  #undef QuadResult
  #undef AbsFehl

} /* End of GaxT */

/* ----------------------------------------------------------------- */

void GaxF (
/*.IX{GaxF}*/
           void *vmblock
          )
/***********************************************************************
* Call this routine when using GaxT instead of Gax in order to continue*
* work on partial results. (in Gax we also call GaxT)                  *
*                                                                      *
* The parameter T contains the address of the pointer being used for   *
* the IntInfo array, that was allocated in GaxT, the parameter N       *
* contains the number of current subintervals.                         *
***********************************************************************/
{
  vmfree(vmblock);
}

#if ! defined(__IBMC__) && ! defined(__IBMCPP__)
#define _Optlink
#endif
int _Optlink GxCmp (const void *a, const void *b)
/*.IX{GxCmp}*/
/***********************************************************************
* used for comparisons in GxSort ()                                    *
***********************************************************************/
{
  REAL diff = FABS ((*(const IntInfo **)a)->FSch) -
              FABS ((*(const IntInfo **)b)->FSch);
  return (diff == 0.0 ? 0 : diff > 0.0 ? -1 : 1);
}

void GxSort (IntInfo** F, int N)
/*.IX{GxSort}*/
/***********************************************************************
* sorts subintervals be decreasing absolute error                      *
***********************************************************************/
{
  qsort (F, N, sizeof (F[0]), GxCmp);
}

int GxQuad (
/*.IX{GxQuad}*/
            IntInfo* Info,
            REAL     Fkt (REAL),
            int      Methode,
            int      N
           )
/***********************************************************************
* Computes the value of the integral of Fkt(x) over the interval       *
* [a,b] using one of the quadratur formulas:                           *
*                                                                      *
* Methode = GxTRAPEZ        : Trapezoidal rule                         *
* Methode = GxGAUSSJORDAN   : N point Gauss formula (N = 2,...,20)     *
* Methode = GxCLENSHAWCURTIS: Clenshaw-Curtis formula with N+1 weights *
* Methode = GxROMBERG       : Romberg method                           *
* Methode = GxNEWTONCOTES   : Newton-Cotes with N (2,...,7) nodes      *
*                                                                      *
* Input parameters:                                                    *
*   IntInfo *Info    Data for the interval [a, b],                     *
*                            computed integral value and associated    *
*                            error                                     *
*   REAL    Fkt      Name of function                                  *
*   int     Methode  Quadrature formula used (see above)               *
*   int     N        number of nodes used for quadrature formula       *
*                                                                      *
* external subroutines used :                                          *
*   ClenCurt, QuaRom, QuaNeC, GxGauss                                  *
*                                                                      *
* global variables used depending on method :                          *
*   StStellen [], Gewichte []                                          *
*                                                                      *
* Author: Uli Eggermann, 10.3.1991                                     *
***********************************************************************/
{
  REAL Schrittweite;
  int  maxRomberg, Ordnung, Fehler = 0;

  switch (Methode)
  {
    case GxGAUSSJORDAN:        /* Gauss-Jordan quadrature of degree N */
           GxGauss (Fkt, N, Info);
           break;

    case GxCLENSHAWCURTIS:          /* Cleshaw-Curtis with N weights  */
           Fehler = GxClenCurt (Fkt, N, Info);
           break;

    case GxROMBERG:                              /* Romberg method   */
           Schrittweite =  ZERO;
           maxRomberg   = 10;
           Fehler = QuaRom (Info->Links, Info->Rechts, 1.0E-8,
                            &maxRomberg, &Schrittweite, Fkt,
                            &Info->QRes, &Info->FSch);
           break;

    case GxTRAPEZ:      /* Trapezoidal rule = Newton-Cotes with N = 1 */
           N = 1;

    case GxNEWTONCOTES:      /* Newton-Cotes with N subintervals      */
           Fehler = QuaNeC (Info->Links, Info->Rechts, 2, N, Fkt,
                            &Info->QRes, &Ordnung, &Info->FSch);

  }                                                         /* switch */
  return (Fehler);
}                                                           /* GxQuad */

void GxGauss (
/*.IX{GxGauss}*/
              REAL     Fkt (REAL),
              int      Anz,
              IntInfo* Info
             )
/***********************************************************************
* Computes the integral of Fkt over the interval [Info->Links,         *
* Info->Rechts] using the Gauss formula of degree Anz.                 *
*                                                                      *
* Parameters:                                                          *
*   REAL    Fkt (REAL)    Name of function                             *
*   int     Anz           number of nodes or weights                   *
*   IntInfo *Info         Interval information                         *
*                                                                      *
* global variables used :                                              *
*   StStellen [], Gewichte []                                          *
*                                                                      *
* subroutine  used :                                                   *
*   GxxGauss                                                           *
*                                                                      *
* Author: Uli Eggermann, 10.3.1991                                     *
***********************************************************************/
{
  REAL m, Res1, Res2;
  m = 0.5 * (Info->Links + Info->Rechts);

  Res1 = GxxGauss (Info->Links, Info->Rechts, Anz, Fkt);
  Res2 = GxxGauss (Info->Links, m,            Anz, Fkt) +
         GxxGauss (m,           Info->Rechts, Anz, Fkt);

  Info->QRes = Res2;
  Info->FSch = (Res1 - Res2) / ((1 << (2 * (Anz - 1) + 3)) - 1);
}                                                          /* GxGauss */

REAL GxxGauss (
/*.IX{GxxGauss}*/
               REAL a,
               REAL b,
               int  Anz,
               REAL Fkt(REAL)
              )
{
  int j, AnzHalbe = Anz / 2;                       /* half the number */
  REAL m, h, t, S = ZERO;
  m = 0.5 * (b + a);                                  /* center point */
  h = 0.5 * (b - a);                                  /* half size    */

  for (j = 0; j < AnzHalbe; j++)
  {
     t = h * StStellen [j];
     S += (Fkt (m - t) + Fkt (m + t)) * Gewichte [j];
  }
  if (Anz % 2)              /* number uneven: Interval center is node */
     S += Fkt (m) * Gewichte [AnzHalbe];
  return S * h;
}

int GxRiEp (
/*.IX{GxRiEp}*/
            REAL  qWert,
            REAL* L,
            REAL* Leps,
            REAL* EpsExe,
            REAL* EpsL,
            int   Num,
            REAL* Eps,
            REAL* AbsFehl,
            REAL* qExe
           )
/***********************************************************************
* Extrapolate current integral value using Richardson extrapolation    *
* and the Epsilon algorithm.                                           *
*                                                                      *
* Parameters :                                                         *
*   REAL    qWert           current integral value                     *
*   REAL   *Eps             relative accuracy                          *
*   REAL   *AbsFehl         absolute accuracy                          *
*   REAL    L   []          current Romberg scheme row                 *
*   REAL    Leps[]          current Epsilon row                        *
*   int     Num             number of calls                            *
*   REAL   *EpsExe          aux variable                               *
*   REAL   *EpsL            Hilfsvariable                              *
*   REAL   *qExe            neuer Quadraturwert                        *
*                                                                      *
* Return value :                                                       *
*   Num + 1                                                            *
***********************************************************************/
{
  REAL   LmSaved, LepsIndSav, Pot4, L1, Temp, Ls;
  int    i, Index;
  /*********************************************************************
  *                                        Richardson  extrapolation   *
  *********************************************************************/
  LmSaved = L [Num];
  Num ++;
  Pot4    = 1.0;
  L [Num] = 0.0;
  L1      = L [1];
  L [1]   = qWert;
  for  (i = 2; i <= Num; i++)
  {
    Pot4 *=  4;
    Temp  =  L [i];
    L [i] = (Pot4 * L [i-1] - L1) / (Pot4 - 1);
    L1    =  Temp;
  }

  /*********************************************************************
  * Extrapolation using the Epsilon algorithm                          *
  *********************************************************************/
  if (Num%2 == 0) { Index = Num-1; LepsIndSav = Leps [Index];   }
  else            { Index = Num;   LepsIndSav = Leps [Index-2]; }

  Leps [Num] = 0.0;
  Ls = L1    = Leps [1];
  Leps [0]   = qWert;
  for (i = 2; i <= Num; i++)
  {
     Temp = Leps [i];
     if (FABS (Leps [i-1] - L1) == 0.0)
       Leps [i] = Ls;
     else
     {
       Leps [i] = Ls + 1 / (Leps [i-1] - L1);
       Ls = L1;
     }
     L1   = Temp;
  }
  *AbsFehl = FABS (L[Num] - L[Num-1]) + FABS  (L[Num] - LmSaved);
  if (Num > 2)
  {
     *EpsExe = FABS (Leps [Index] - LepsIndSav)
             + FABS (Leps [Index] - Leps [Index - 2]);
     *qExe   = Leps [Index];
  }
  if (*AbsFehl < *EpsExe) *qExe = L [Num];
  *EpsL = *Eps * FABS (*qExe);
  if (*AbsFehl > *EpsExe) *AbsFehl  = *EpsExe;

  return (Num);

}                                                           /* GxRiEp */

int GxClenCurt (
/*.IX{GxClenCurt}*/
                REAL     f(REAL),
                int      n,
                IntInfo* Info
               )
/***********************************************************************
* Computes the integral of f over the interval [Info->Links,           *
* Info->Rechts] using a Clenshaw-Curtis formula with  n+1  nodes.      *
*                                                                      *
* Parameters :                                                         *
*   REAL     f (REAL x)     Name of function                           *
*   int      n              n + 1, number of nodes                     *
*   IntInfo *Info           Intervall information                      *
*                                                                      *
* global variables used :                                              *
*   StStellen [], Gewichte []                                          *
*                                                                      *
* Author: Uli Eggermann, 10.3.1991                                     *
***********************************************************************/
{
   REAL u[3], Res1, Res2;
   int i;
   u[0] = Info->Links; u[1] = u[2] = Info->Rechts;
   if ((i = ClenCurt (f,1,u,n,StStellen,Gewichte,&Res1)) != 0)
     return (i);
   u[1] = .5 * (u[0] + u[2]);
   if ((i = ClenCurt (f,2,u,n,StStellen,Gewichte,&Res2)) != 0)
     return (i);
   Info->QRes =  Res2;
   Info->FSch = (Res1 - Res2) / ((1 << (n + 2)) + 1);
   return (0);
}

/* ---------------------------- END gax.c --------------------------- */
