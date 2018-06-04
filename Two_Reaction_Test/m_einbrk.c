#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ----------------------- MODULE m_einbrk.c ------------------------ */

#include <basis.h>         /*  for  umleiten, fprintf, stderr, scanf, */
                           /*       printf, NULL, FABS, REAL, LZS,    */
                           /*       fehler_melden, LZP, TRUE          */
#include <vmblock.h>       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include <einb_rk.h>       /*  for  einb_rk                           */
#include <t_dgls.h>        /*  for  bsptyp, dgls_waehlen              */
#ifdef __TURBOC__
#ifdef __MSDOS__
#include <alloc.h>         /*  for  coreleft                          */
#endif
#endif
#ifndef VOLLTEST                            /* Test with input file?  */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
* Test program for the function einb_rk() from the module einb_rk      *
* to solve a first order ordinary system of differential equations     *
* using one of 23 embedding formulas.                                  *
*                                                                      *
* Scope of program :                                                   *
* ==================                                                   *
* The program reads the input data from the file stdin and writes      *
* output onto the file stdout.                                         *
* If the first command line entry exists, it is interpreted to be the  *
* input file and this is assigned to stdin.                            *
* The second command line parameter is treated analogously wrt. stdout.*
* Calls for input are directed to stderr which also collects error     *
* messages.                                                            *
*                                                                      *
* After reading the input, it is put out for control purposes; then    *
* the embedding formula is applied  and the computed results are put   *
* out.                                                                 *
*                                                                      *
* To solve a differential equation, please proceed as in example       *
* t_dgls.c .                                                           *
*                                                                      *
* Construction of input data files :                                   *
* ==================================                                   *
* bspnummer  Number of DE system in  t_dgls.c                          *
* x0         initial x-value                                           *
* xend       final desired x-value                                     *
* y0[0]   \  initial y-value at x0                                     *
* ...      >                                                           *
* y0[n-1] /                                                            *
* epsabs     desired absolute error bound                              *
* epsrel     desired relative error bound                              *
* fmax       maximal number of calls of right hand side                *
* neinb      Number of embedding formula  (0, 1,..., 14)               *
* hullstp    Step size control according to Hull (1) or using the      *
*            standard control (0)                                      *
*                                                                      *
* The number n of differential equations follows from the number of    *
* the DE system chosen; it is stored in conjunction with the right     *
* hand side function in t_dgls.c .                                     *
***********************************************************************/

{
  REAL   x0,             /* initial x-value                           */
         xend,           /* desired final x-value                     */
         *y0,            /* [0..n-1] vector: appr. initial y-value    */
         *yex,           /* [0..n-1] vector: exact y-value            */
         epsabs,         /* desired absolue error bound               */
         epsrel;         /* ditto for the relative error bound        */
  long   fmax,           /* maximal number of calls of the right      */
                         /* hand side in einb_rk()                    */
         aufrufe;        /* actual number of calls                    */
  int    bspnummer,      /* Number of the DE system from t_dgls.c     */
         n,              /* number of DEs                             */
         fehler,         /* error code from einb_rk()                 */
         i,              /* loop counter                              */
         neinb,          /* Number of the embedding formula           */
         hullstp;        /* Step size control according to Hull?      */
  bsptyp  *beispiel;     /* pointer to the structure describing the   */
                         /* DE system at hand                         */
  void    *vmblock;      /* List of dynamic allocations               */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard ones           */
    return fehler;                                          /* 1 or 2 */

  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "Example :                                       ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = dgls_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("non existing example ",
                  0, __FILE__, __LINE__);
    return 3;
  }

  n = beispiel->n;
  vmblock = vminit();                      /* initialize storage  */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  yex = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of storage", 0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "initial x-value x0:      ");
#endif
  scanf("%"LZS"f", &x0);

#ifdef INTERAKTIV
  fprintf(stderr, "desired x-value xend:   ");
#endif
  scanf("%"LZS"f", &xend);

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "Function value y0[%d] at x0:                  ",
            i);
#endif
    scanf("%"LZS"f", y0 + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "absolute error bound epsabs:                 ");
#endif
  scanf("%"LZS"f", &epsabs);

#ifdef INTERAKTIV
  fprintf(stderr, "relative error bound epsrel:                 ");
#endif
  scanf("%"LZS"f", &epsrel);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of calls for r.h. side: ");
#endif
  scanf("%ld", &fmax);

#ifdef INTERAKTIV
  fprintf(stderr, "Number of embedding formula (0...22):          ");
#endif
  scanf("%d", &neinb);

#ifdef INTERAKTIV
  fprintf(stderr, "general (0) or Hull step size control (1):     ");
#endif
  scanf("%d", &hullstp);


  /* ------------ put out input ------------------------------------ */

  printf("\n"
         "Solution of a first order ordinary DE system\n"
         "============================================\n"
         "using one of 23 embedding formulas:\n"
         "===================================\n"
         "\n\n"
         "System of DEs:\n"
         "--------------\n"
         "%s\n\n"
         "Input data:\n"
         "-----------\n"
         "Example  = %24d\n"
         "n        = %24d\n"
         "x0       = %24.15"LZP"e\n"
         "xend     = %24.15"LZP"e\n"
         "epsabs   = %24.15"LZP"e\n"
         "epsrel   = %24.15"LZP"e\n"
         "fmax     = %24ld\n"
         "neinb    = %24d\n"
         "hullstp  = %24d\n",
         (*beispiel->dgl_text)(), bspnummer, n, x0, xend, epsabs,
         epsrel, fmax, neinb, hullstp
        );

  for (i = 0; i < n; i++)
    printf("y0[%d]    = %24.15"LZP"e\n", i, y0[i]);


  /* ------------ Solve system of DEs ------------------------------- */
#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "\nfree before:  %u\n", coreleft());
#endif
#endif
#endif

  fehler = einb_rk(&x0, xend, n, beispiel->rechte_seite, y0, epsabs,
                   epsrel, neinb, hullstp, TRUE, FALSE, fmax, &aufrufe);

#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "free after: %u\n", coreleft());
#endif
#endif
#endif

  /* -------------------- put out results --------------------------- */

  printf("\n\n"
         "Output data:\n"
         "------------\n"
         "error code from einb_rk():         %24d\n"
         "number of calls of  dgl():         %24ld\n"
         "final x-value x =                  %24.15"LZP"e\n\n",
         fehler, aufrufe, x0
        );

  for (i = 0; i < n; i++)
    printf("approximate solution y%d(x) = %24.15"LZP"e\n",
           i + 1, y0[i]);

  if (beispiel->exakte_loesung != NULL)       /* exact solution       */
  {                                           /* available?           */
    (*beispiel->exakte_loesung)(x0, yex);
    printf("\n");
    for (i = 0; i < n; i++)
      printf("'exact' solution  y%d(x)    = %24.15"LZP"e\n",
             i + 1, yex[i]);
    printf("\nDifference:  approx sol. - 'exact' solution:\n");
    for (i = 0; i < n; i++)
      printf("%24.15"LZP"g\n", y0[i] - yex[i]);
  }

  if (fehler != 0)
  {
    fehler_melden("einb_rk()", fehler, __FILE__, __LINE__);
    return fehler;
  }


  return 0;
}

#else                                       /* Test from FNUM?       */



/* Test examples */

void dgl1(REAL x, REAL y[], REAL f[])
{
  f[0] = THREE * POW(FABS(y[0]), TWO / THREE) * COS(x);
}

void exakt1(REAL x, REAL y[])
{
  y[0] = POW(TWO + SIN(x), THREE);
}

void dgl2(REAL x, REAL y[], REAL f[])
{
  x    = x;
  f[0] = y[1];
  f[1] = -y[0] / POW(SQRT(y[0] * y[0] + y[2] * y[2]), THREE);
  f[2] = y[3];
  f[3] = -y[2] / POW(SQRT(y[0] * y[0] + y[2] * y[2]), THREE);
}

void exakt2(REAL x, REAL y[])
{
  y[0] = COS(x);
  y[1] = -SIN(x);
  y[2] = SIN(x);
  y[3] = COS(x);
}

int main(void)

/***********************************************************************
* Test program for the funktion einbrk()      (from FNUM)              *
*                                                                      *
* Test examples:                                                       *
* We solve two IVPs using the embedding formulas.                      *
*                                                                      *
* Test problem 1:                                                      *
*              Y'=COS(X)*3*Y^(2/3)           Y(0)=8    I=[0,20]        *
*                                                                      *
*              Every embedding formula is used to solve the problem    *
*              via step size control HULL and CIVPS.                   *
*                                                                      *
* Test problem 2:                                                      *
*              Y1'=Y2                        Y1(0)=1   I=[0,20]        *
*              Y2'=-Y1/(SQRT(Y1^2+Y3^2))^3   Y2(0)=0                   *
*              Y3'=Y4                        Y3(0)=0                   *
*              Y4'=-Y3/(SQRT(Y1^2+Y3^2))^3   Y4(0)=1                   *
*                                                                      *
*              The IVP is solved via the embedding formula RK8(7)13M   *
*              and step size control HULL. We compute values of the    *
*              solution at X = 5, 10, 15 and 20.                       *
*                                                                      *
* Here we keep integrating until we reach XEND. If `fehler==-2' (max.  *
* number of integration steps reached), we repeat the call to          *
* enb_rk().                                                            *
*                                                                      *
*                                                                      *
* Results:                                                             *
*                                                                      *
 C[                                                                  ]*
 C[ TEST 1:                                                          ]*
 C[                                                                  ]*
 C[ ABSERR= 1.0000000e-05    RELERR= 0.0000000e+00                   ]*
 C[                                                                  ]*
 C[ STEP CONTROL: 0                                                  ]*
 C[                                                                  ]*
 C[ METHOD | IFEHL | X     | APPROXIMATE VALUE | EXACT SOLUTION      ]*
 C[ -------+-------+-------+-------------------+------------------   ]*
 C[  0     I  -2   I 18.35 I 3.49022668846e+00 I 3.49022576600e+00   ]*
 C[  0     I   0   I 20.00 I 2.47170720660e+01 I 2.47170687870e+01   ]*
 C[  1     I   0   I 20.00 I 2.47170475491e+01 I 2.47170687870e+01   ]*
 C[  2     I   0   I 20.00 I 2.47174720111e+01 I 2.47170687870e+01   ]*
 C[  3     I   0   I 20.00 I 2.47171280422e+01 I 2.47170687870e+01   ]*
 C[  4     I   0   I 20.00 I 2.47170169782e+01 I 2.47170687870e+01   ]*
 C[  5     I   0   I 20.00 I 2.47172344214e+01 I 2.47170687870e+01   ]*
 C[  6     I   0   I 20.00 I 2.47177931998e+01 I 2.47170687870e+01   ]*
 C[  7     I   0   I 20.00 I 2.47170845846e+01 I 2.47170687870e+01   ]*
 C[  8     I   0   I 20.00 I 2.47176421879e+01 I 2.47170687870e+01   ]*
 C[  9     I   0   I 20.00 I 2.47172295381e+01 I 2.47170687870e+01   ]*
 C[ 10     I  -2   I  8.50 I 2.18662207854e+01 I 2.18659592493e+01   ]*
 C[ 10     I  -2   I 15.31 I 1.35733544141e+01 I 1.35729680506e+01   ]*
 C[ 10     I   0   I 20.00 I 2.47180039975e+01 I 2.47170687870e+01   ]*
 C[ 11     I   0   I 20.00 I 2.47170883243e+01 I 2.47170687870e+01   ]*
 C[ 12     I   0   I 20.00 I 2.47170736482e+01 I 2.47170687870e+01   ]*
 C[ 13     I   0   I 20.00 I 2.47173746805e+01 I 2.47170687870e+01   ]*
 C[ 14     I   0   I 20.00 I 2.47176148592e+01 I 2.47170687870e+01   ]*
 C[ 15     I   0   I 20.00 I 2.47170234046e+01 I 2.47170687870e+01   ]*
 C[ 16     I   0   I 20.00 I 2.47170701143e+01 I 2.47170687870e+01   ]*
 C[ 17     I   0   I 20.00 I 2.47170363752e+01 I 2.47170687870e+01   ]*
 C[ 18     I   0   I 20.00 I 2.47169824664e+01 I 2.47170687870e+01   ]*
 C[ 19     I   0   I 20.00 I 2.47170882811e+01 I 2.47170687870e+01   ]*
 C[ 20     I   0   I 20.00 I 2.47170863396e+01 I 2.47170687870e+01   ]*
 C[ 21     I   0   I 20.00 I 2.47170788593e+01 I 2.47170687870e+01   ]*
 C[ 22     I   0   I 20.00 I 2.47170961690e+01 I 2.47170687870e+01   ]*
 C[                                                                  ]*
 C[ STEP CONTROL: 1                                                  ]*
 C[                                                                  ]*
 C[ METHOD | IFEHL | X     | APPROXIMATE VALUE | EXACT SOLUTION      ]*
 C[ -------+-------+-------+-------------------+------------------   ]*
 C[  0     I   0   I 20.00 I 2.47171948788e+01 I 2.47170687870e+01   ]*
 C[  1     I   0   I 20.00 I 2.47168984960e+01 I 2.47170687870e+01   ]*
 C[  2     I   0   I 20.00 I 2.47178976577e+01 I 2.47170687870e+01   ]*
 C[  3     I   0   I 20.00 I 2.47172451874e+01 I 2.47170687870e+01   ]*
 C[  4     I   0   I 20.00 I 2.47167894645e+01 I 2.47170687870e+01   ]*
 C[  5     I   0   I 20.00 I 2.47174704023e+01 I 2.47170687870e+01   ]*
 C[  6     I   0   I 20.00 I 2.47187273348e+01 I 2.47170687870e+01   ]*
 C[  7     I   0   I 20.00 I 2.47170502642e+01 I 2.47170687870e+01   ]*
 C[  8     I   0   I 20.00 I 2.47180503733e+01 I 2.47170687870e+01   ]*
 C[  9     I   0   I 20.00 I 2.47172563102e+01 I 2.47170687870e+01   ]*
 C[ 10     I   0   I 20.00 I 2.47188606114e+01 I 2.47170687870e+01   ]*
 C[ 11     I   0   I 20.00 I 2.47170967600e+01 I 2.47170687870e+01   ]*
 C[ 12     I   0   I 20.00 I 2.47170766655e+01 I 2.47170687870e+01   ]*
 C[ 13     I   0   I 20.00 I 2.47175762901e+01 I 2.47170687870e+01   ]*
 C[ 14     I   0   I 20.00 I 2.47173916771e+01 I 2.47170687870e+01   ]*
 C[ 15     I   0   I 20.00 I 2.47170238160e+01 I 2.47170687870e+01   ]*
 C[ 16     I   0   I 20.00 I 2.47170748699e+01 I 2.47170687870e+01   ]*
 C[ 17     I   0   I 20.00 I 2.47169845422e+01 I 2.47170687870e+01   ]*
 C[ 18     I   0   I 20.00 I 2.47169723595e+01 I 2.47170687870e+01   ]*
 C[ 19     I   0   I 20.00 I 2.47171161729e+01 I 2.47170687870e+01   ]*
 C[ 20     I   0   I 20.00 I 2.47170651150e+01 I 2.47170687870e+01   ]*
 C[ 21     I   0   I 20.00 I 2.47170838513e+01 I 2.47170687870e+01   ]*
 C[ 22     I   0   I 20.00 I 2.47170647268e+01 I 2.47170687870e+01   ]*
 C[                                                                  ]*
 C[                                                                  ]*
 C[ TEST 2:                                                          ]*
 C[                                                                  ]*
 C[ ABSERR= 0.0000000e+00    RELERR= 1.0000000e-06                   ]*
 C[                                                                  ]*
 C[ METHOD: RK8(7)13M   STEP CONTROL: HULL                           ]*
 C[                                                                  ]*
 C[ IFEHL=  0     X=   5.000                                         ]*
 C[                                                                  ]*
 C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
 C[          ---------+----------------------+---------------------- ]*
 C[          Y(1)     I  2.8366239277040e-01 I  2.8366218546323e-01  ]*
 C[          Y(2)     I  9.5892436880233e-01 I  9.5892427466314e-01  ]*
 C[          Y(3)     I -9.5892376020970e-01 I -9.5892427466314e-01  ]*
 C[          Y(4)     I  2.8366282721367e-01 I  2.8366218546323e-01  ]*
 C[                                                                  ]*
 C[ IFEHL=  0     X=  10.000                                         ]*
 C[                                                                  ]*
 C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
 C[          ---------+----------------------+---------------------- ]*
 C[          Y(1)     I -8.3906907933749e-01 I -8.3907152907645e-01  ]*
 C[          Y(2)     I  5.4402454729505e-01 I  5.4402111088937e-01  ]*
 C[          Y(3)     I -5.4402376544152e-01 I -5.4402111088937e-01  ]*
 C[          Y(4)     I -8.3906962583324e-01 I -8.3907152907645e-01  ]*
 C[                                                                  ]*
 C[ IFEHL=  0     X=  15.000                                         ]*
 C[                                                                  ]*
 C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
 C[          ---------+----------------------+---------------------- ]*
 C[          Y(1)     I -7.5969379249209e-01 I -7.5968791285882e-01  ]*
 C[          Y(2)     I -6.5028028464851e-01 I -6.5028784015712e-01  ]*
 C[          Y(3)     I  6.5027969418122e-01 I  6.5028784015712e-01  ]*
 C[          Y(4)     I -7.5969479177045e-01 I -7.5968791285882e-01  ]*
 C[                                                                  ]*
 C[ IFEHL=  0     X=  20.000                                         ]*
 C[                                                                  ]*
 C[          SOLUTION I APPROXIMATE          I EXACT                 ]*
 C[          ---------+----------------------+---------------------- ]*
 C[          Y(1)     I  4.0806395049234e-01 I  4.0808206181339e-01  ]*
 C[          Y(2)     I -9.1295367373865e-01 I -9.1294525072763e-01  ]*
 C[          Y(3)     I  9.1295190308152e-01 I  9.1294525072763e-01  ]*
 C[          Y(4)     I  4.0806472045855e-01 I  4.0808206181339e-01  ]*
*                                                                      *
***********************************************************************/
{
  REAL    x0,            /* initial x-value                           */
          xend,          /* desired final x-value                     */
          y0[4],         /* [0..n-1] vector: appr. initial y-value    */
          yex[4],        /* [0..n-1] vector: exact solution           */
          epsabs,        /* desired absolue error bound               */
          epsrel;        /* ditto for the relative error bound        */
  long   fmax,           /* maximal number of calls of the right      */
                         /* hand side in einb_rk()                    */
         aufrufe;        /* actual number of calls                    */
  int    n,              /* number of DEs                             */
         fehler,         /* error code from einb_rk()                 */
         neinb,          /* Number of the embedding formula           */
         hullstp,        /* Step size control according to Hull?      */
         i;              /* loop counter                              */

  /* Test problem 1 */

  epsabs = (REAL)1.0e-5;
  epsrel = ZERO;
  fmax   = 32000l;
  n      = 1;
  printf(" C[%68s\n"
         " C[ TEST 1:%60s\n"
         " C[%68s\n"
         " C[ ABSERR=%14.7"LZP"e    RELERR=%14.7"LZP"e%21s\n",
         "]*", "]*", "]*", epsabs, epsrel, "]*"
        );

  /* Loop to choose step size control */

  for (i = 0; i <= 1; i++)
  {
    int j;
    hullstp = i;
    printf(" C[%68s\n"
           " C[ STEP CONTROL: %1d%52s\n"
           " C[%68s\n"
           " C[ METHOD | IFEHL | X     | APPROXIMATE VALUE |"
           " EXACT SOLUTION      ]*\n"
           " C[ -------+-------+-------+-------------------+"
           "------------------   ]*\n",
           "]*", hullstp, "]*", "]*"
          );

    /* Loop to choose embedding formula */

    for (j = 0; j <= 22; j++)
    {
      int neu;
      int save;
      neinb = j;
      neu   = TRUE;
      save  = TRUE;
      x0    = ZERO;
      xend  = (REAL)20.0;
      y0[0] = EIGHT;
      do
      {
        fehler = einb_rk(&x0, xend, n, dgl1, y0, epsabs, epsrel, neinb,
                         hullstp, neu, save, fmax, &aufrufe);
        exakt1(x0, yex);
        printf(" C[ %2d     I %3d   I%6.2"LZP"f I%18.11"LZP"e"
               " I%18.11"LZP"e   ]*\n",
               neinb, fehler, x0, y0[0], yex[0]
              );
        neu = FALSE;
      /* if maximal allowable number of integration steps has been   */
      /* reached, repeat call of RKTRB.                              */
      } while (fehler == -2);
    }
  }

  /* Test problem 2 */

  {
  int neu;
  int save;
  epsabs = ZERO;
  epsrel = (REAL)1.0e-6;
  n      = 4;
  neinb  = 19;
  neu    = TRUE;
  save   = TRUE;
  x0     = ZERO;
  y0[0]  = ONE;
  y0[1]  = ZERO;
  y0[2]  = ZERO;
  y0[3]  = ONE;

  printf(" C[%68s\n"
         " C[%68s\n"
         " C[ TEST 2:%60s\n"
         " C[%68s\n"
         " C[ ABSERR=%14.7"LZP"e    RELERR=%14.7"LZP"e%21s\n"
         " C[%68s\n"
         " C[ METHOD: RK8(7)13M   STEP CONTROL: HULL%29s\n",
         "]*", "]*", "]*", "]*", epsabs, epsrel, "]*", "]*", "]*"
        );

  /* Loop to compute value of the solution at             */
  /* X=5, 10, 15 and 20                                   */

  for (xend = FIVE; xend < (REAL)20.1; xend += FIVE)
  {
    do
    {
      fehler = einb_rk(&x0, xend, n, dgl2, y0, epsabs, epsrel, neinb,
                       hullstp, neu, save, fmax, &aufrufe);
      exakt2(x0, yex);
      printf(" C[%68s\n"
             " C[ IFEHL=%3d     X=%8.3"LZP"f%43s\n"
             " C[%68s\n"
             " C[          SOLUTION I APPROXIMATE          I EXACT"
             "                 ]*\n"
             " C[          ---------+----------------------+------"
             "---------------- ]*\n",
             "]*", fehler, x0, "]*", "]*"
            );
      for (i = 0; i < n; i++)
        printf(
               " C[          Y(%1d)     I %20.13"LZP"e I "
               "%20.13"LZP"e  ]*\n",
               i + 1, y0[i], yex[i]
              );
      neu = FALSE;
      /* if maximal allowable number of integration steps has been   */
      /* reached, repeat call of RKTRB.                              */
    } while (fehler == -2);
  }
  }

  return 0;
}
#endif

/* ------------------------- END  m_einbrk.c ------------------------ */
