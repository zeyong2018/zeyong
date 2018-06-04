#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE m_rwp.c ------------------------- */

#include <basis.h>       /*  for  umleiten, fprintf, stderr, printf,  */
                         /*       scanf, NULL, REAL, LZS, LZP,        */
                         /*       fehler_melden                       */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vmfree,        */
                         /*       vminit, VEKTOR                      */
#include <rwp.h>         /*  for  rwp                                 */
#include <t_rwp.h>       /*  for  bsptyp, rwp_waehlen                 */
#ifdef __TURBOC__
#ifdef __MSDOS__
#include <alloc.h>       /*  for  coreleft                            */
#endif
#endif
#ifndef VOLLTEST                            /* test with input file ? */



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
*                                                                      *
* Test program for the function rwp() from the module rwp for solving  *
* two point boundary value problems of first order using the shooting  *
* method.                                                              *
*                                                                      *
* Scope of program :                                                   *
* ==================                                                   *
* Input is read from stdin and output is gathered in stdout.           *
* If there is one file in the command line this file is assigned to the*
* standard input file stdin. Similar for the second file name.         *
* The input calls are registered in stderr together with the error     *
* messages.                                                            *
*                                                                      *
* After reading the input, the program print the input out, applies    *
* the shooting method to the data and gives the results.               *
*                                                                      *
* To solve other differential equations one can follow the examples    *
* in t_rwp.c.                                                          *
*                                                                      *
* Construction of input file:                                          *
* ===========================                                          *
* bspnummer       Number of differntial equation system from  t_rwp.c  *
* a               left end point                                       *
* b               right end point                                      *
* h               starting step size                                   *
* epsawp          desired absolute error bound                         *
* epsrb           desired relative error bound                         *
* fmax            maximal number of calls of right hand side of DE     *
* itmax           maximal number of Newton iterations                  *
* y_start[0]   \                                                       *
* ...           > Starting value for the initial value of the          *
* y_start[n-1] /  solution at a                                        *
* awpnumm         Number of the desired IVP solver                     *
*                                                                      *
* The size of the system of differential equations follows implicitly  *
* from the function declarations in t_rwp.c.                           *
***********************************************************************/

{
  REAL   a,
         b,
         h,
         epsrb,
         epsawp,
         *y_start;
  long   fmax;
  int    bspnummer,
         n,
         itmax,
         fehler,         /* error code for rwp()                      */
         i,              /* loop counter                              */
         aufrufe,        /* number of calls of right hand side of DE  */
                         /* in rwp()                                  */
         awpnumm;        /* Number of desired IVP solver              */
  bsptyp *beispiel;      /* pointer to stucture that describes the    */
                         /* boundary value problem                    */
  void   *vmblock;       /* List of dynamically allocated vectors     */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input -------------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "Example:                                          ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = rwp_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("improper number for the example",
                  0, __FILE__, __LINE__);
    return 3;
  }

  n = beispiel->n;
  vmblock = vminit();                 /* allocate storage   */
  y_start = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "left end point:                                   "
                  "     ");
#endif
  scanf("%"LZS"f", &a);

#ifdef INTERAKTIV
  fprintf(stderr, "right end point:                                  "
                  "     ");
#endif
  scanf("%"LZS"f", &b);

#ifdef INTERAKTIV
  fprintf(stderr, "initial step size h:                              "
                  "     ");
#endif
  scanf("%"LZS"f", &h);

#ifdef INTERAKTIV
  fprintf(stderr, "accuracy bound for initial value problem epsawp:  "
                  "     ");
#endif
  scanf("%"LZS"f", &epsawp);

#ifdef INTERAKTIV
  fprintf(stderr, "accuracy bound for boundary value problem epsrb:  "
                  "     ");
#endif
  scanf("%"LZS"f", &epsrb);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of calls of right hand side of DE: "
                  "     ");
#endif
  scanf("%ld", &fmax);

#ifdef INTERAKTIV
  fprintf(stderr, "maximal number of Newton iterations:              "
                  "     ");
#endif
  scanf("%d", &itmax);

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "initial value  y_start[%d] at a:"
                    "                   ",
                    i);
#endif
    scanf("%"LZS"f", y_start + i);
  }

#ifdef INTERAKTIV
  fprintf(stderr, "Number of initial value solver (1 - 5):            "
                  "           ");
#endif
  scanf("%d", &awpnumm);


  /* ------------ print out input data ----------------------------- */

  printf("\n"
         "Solution of a two point boundary value problem\n"
         "==============================================\n"
         "of first order via the shooting method\n"
         "======================================\n\n\n"
         "Given system of differential equations:\n"
         "---------------------------------------\n"
         "%s\n"
         "Boundary conditions:\n"
         "--------------------\n"
         "%s\n"
         "input data:\n"
         "-----------\n"
         "Example    = %24d\n"
         "n          = %24d\n"
         "a          = %24.15"LZP"e\n"
         "b          = %24.15"LZP"e\n"
         "h          = %24.15"LZP"e\n"
         "epsawp     = %24.15"LZP"e\n"
         "epsrb      = %24.15"LZP"e\n"
         "fmax       = %24ld\n"
         "itmax      = %24d\n"
         "awpnumm    = %24d\n",
         beispiel->dgl_text(), beispiel->randtext(),
         bspnummer, n, a, b, h, epsawp, epsrb,
         fmax, itmax, awpnumm);

  for (i = 0; i < n; i++)
    printf("y_start[%d] = %24.15"LZP"e\n", i, y_start[i]);


  /* ----------------- solve boundary value problem ----------------- */

#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "\nfree before:  %u\n", coreleft());
#endif
#endif
#endif
  fehler = rwp(a, b, h, y_start, n, beispiel->rechte_seite,
               beispiel->rand, awpnumm, epsawp, epsrb, fmax, itmax,
               &aufrufe);
#ifdef INTERAKTIV
#ifdef __TURBOC__
#ifdef __MSDOS__
  fprintf(stderr, "free after: %u\n", coreleft());
#endif
#endif
#endif

  if (fehler != 0)
  {
    fehler_melden("rwp()", 10 + fehler, __FILE__, __LINE__);
    return 10 + fehler;
  }


 /* -------------------- print results ---------------------------- */

  printf("\n\n"
         "Output data:\n"
         "-----------\n"
         "error code from rwp():          %24d\n"
         "number of Newton iterations:    %24d\n",
         fehler, aufrufe);

  for (i = 0; i < n; i++)
    printf("approximate solution y%d(a) =    %24.15"LZP"e\n",
           i + 1, y_start[i]);


  return 0;
}
#else



static void dgl(REAL x, REAL *y, REAL *f)
{
  x    = x;
  f[0] = y[1];
  f[1] = -y[0] * y[0] * y[0];
}

static void randbed(REAL *ya, REAL *yb, REAL *r)
{
  r[0] = ya[0];
  r[1] = yb[0];
}

int main(void)

/***********************************************************************
* Test program for rwp()               ( from FNUM )                   *
*                                                                      *
* Solve a boundary value problem for a first order DE system of the    *
* form:                                                                *
*                                                                      *
*                   Y1' = F1(X,Y1,Y2,...,YN)                           *
*                   Y2' = F2(X,Y1,Y2,...,YN)                           *
*                             ...                                      *
*                   YN' = FN(X,Y1,Y2,...,YN)                           *
*                                                                      *
* via the shooting method by determining an approximation for the      *
* initial value Y(A).                                                  *
*                                                                      *
* The given test data should produce the following results:            *
* ( Borland C++ 1.0 and OS/2 with `REAL=double'):                      *
*                                                                      *
 C[
 C[
 C[
 C[ TEST EXAMPLE (WITH IVP PROGRAM NO 1):
 C[ =============
 C[ ANALYZED SET OF DIFFERENTIAL EQUATIONS:
 C[ ---------------------------------------
 C[          Y'( 1) = Y(2)
 C[          Y'( 2) = -Y(1)**3
 C[ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)
 C[   AT POINT  A =  0.000E+00  : 0.000E+00
 C[   AT POINT  B =  1.000E+00  : 0.000E+00
 C[
 C[ AT LOCATION  A =  0.000E+00 THE FOLLOWING VALUES ARE PROVIDED:
 C[
 C[   Y( 1) =  0.000E+00
 C[   Y( 2) =  1.200E+01
 C[
 C[ REQUIRED PARAMETER:
 C[ ------------------:
 C[ -PRECISION FOR THE IVP =    1.0000000000E-07
 C[ -PRECISION FOR THE BVP =    1.0000000000E-08
 C[ -STEP WIDTH FOR THE IVP=    1.0000000000E-02
 C[ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000
 C[ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000
 C[
 C[ SOLUTION OF THE PROBLEM:
 C[ ------------------------
 C[ THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED
 C[ IN POINT A = 0.000E+00 AS FOLLOWS:
 C[
 C[   Y( 1) =     0.0000000000E+00
 C[   Y( 2) =     9.7229810242E+00
 C[
 C[ THIS REQUIRES     4 NEWTON-ITERATIONS!
 C[
 C[ DETERMINATION COMPLETED AS PLANNED!
 C[
 C[
 C[
 C[ TEST EXAMPLE (WITH IVP PROGRAM NO 2):
 C[ =============
 C[ ANALYZED SET OF DIFFERENTIAL EQUATIONS:
 C[ ---------------------------------------
 C[          Y'( 1) = Y(2)
 C[          Y'( 2) = -Y(1)**3
 C[ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)
 C[   AT POINT  A =  0.000E+00  : 0.000E+00
 C[   AT POINT  B =  1.000E+00  : 0.000E+00
 C[
 C[ AT LOCATION  A =  0.000E+00 THE FOLLOWING VALUES ARE PROVIDED:
 C[
 C[   Y( 1) =  0.000E+00
 C[   Y( 2) =  1.200E+01
 C[
 C[ REQUIRED PARAMETER:
 C[ ------------------:
 C[ -PRECISION FOR THE IVP =    1.0000000000E-07
 C[ -PRECISION FOR THE BVP =    1.0000000000E-08
 C[ -STEP WIDTH FOR THE IVP=    1.0000000000E-02
 C[ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000
 C[ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000
 C[
 C[ SOLUTION OF THE PROBLEM:
 C[ ------------------------
 C[ THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED
 C[ IN POINT A = 0.000E+00 AS FOLLOWS:
 C[
 C[   Y( 1) =     0.0000000000E+00
 C[   Y( 2) =     9.7229809702E+00
 C[
 C[ THIS REQUIRES     3 NEWTON-ITERATIONS!
 C[
 C[ DETERMINATION COMPLETED AS PLANNED!
 C[
 C[
 C[
 C[ TEST EXAMPLE (WITH IVP PROGRAM NO 3):
 C[ =============
 C[ ANALYZED SET OF DIFFERENTIAL EQUATIONS:
 C[ ---------------------------------------
 C[          Y'( 1) = Y(2)
 C[          Y'( 2) = -Y(1)**3
 C[ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)
 C[   AT POINT  A =  0.000E+00  : 0.000E+00
 C[   AT POINT  B =  1.000E+00  : 0.000E+00
 C[
 C[ AT LOCATION  A =  0.000E+00 THE FOLLOWING VALUES ARE PROVIDED:
 C[
 C[   Y( 1) =  0.000E+00
 C[   Y( 2) =  1.200E+01
 C[
 C[ REQUIRED PARAMETER:
 C[ ------------------:
 C[ -PRECISION FOR THE IVP =    1.0000000000E-07
 C[ -PRECISION FOR THE BVP =    1.0000000000E-08
 C[ -STEP WIDTH FOR THE IVP=    1.0000000000E-02
 C[ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000
 C[ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000
 C[
 C[ SOLUTION OF THE PROBLEM:
 C[ ------------------------
 C[ THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED
 C[ IN POINT A = 0.000E+00 AS FOLLOWS:
 C[
 C[   Y( 1) =     4.0541851510E-23
 C[   Y( 2) =     9.7229804225E+00
 C[
 C[ THIS REQUIRES     4 NEWTON-ITERATIONS!
 C[
 C[ DETERMINATION COMPLETED AS PLANNED!
 C[
 C[
 C[
 C[ TEST EXAMPLE (WITH IVP PROGRAM NO 4):
 C[ =============
 C[ ANALYZED SET OF DIFFERENTIAL EQUATIONS:
 C[ ---------------------------------------
 C[          Y'( 1) = Y(2)
 C[          Y'( 2) = -Y(1)**3
 C[ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)
 C[   AT POINT  A =  0.000E+00  : 0.000E+00
 C[   AT POINT  B =  1.000E+00  : 0.000E+00
 C[
 C[ AT LOCATION  A =  0.000E+00 THE FOLLOWING VALUES ARE PROVIDED:
 C[
 C[   Y( 1) =  0.000E+00
 C[   Y( 2) =  1.200E+01
 C[
 C[ REQUIRED PARAMETER:
 C[ ------------------:
 C[ -PRECISION FOR THE IVP =    1.0000000000E-07
 C[ -PRECISION FOR THE BVP =    1.0000000000E-08
 C[ -STEP WIDTH FOR THE IVP=    1.0000000000E-02
 C[ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000
 C[ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000
 C[
 C[ SOLUTION OF THE PROBLEM:
 C[ ------------------------
 C[ THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED
 C[ IN POINT A = 0.000E+00 AS FOLLOWS:
 C[
 C[   Y( 1) =    -4.3143378054E-21
 C[   Y( 2) =     9.7229810232E+00
 C[
 C[ THIS REQUIRES     3 NEWTON-ITERATIONS!
 C[
 C[ DETERMINATION COMPLETED AS PLANNED!
 C[
 C[
 C[
 C[ TEST EXAMPLE (WITH IVP PROGRAM NO 5):
 C[ =============
 C[ ANALYZED SET OF DIFFERENTIAL EQUATIONS:
 C[ ---------------------------------------
 C[          Y'( 1) = Y(2)
 C[          Y'( 2) = -Y(1)**3
 C[ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)
 C[   AT POINT  A =  0.000E+00  : 0.000E+00
 C[   AT POINT  B =  1.000E+00  : 0.000E+00
 C[
 C[ AT LOCATION  A =  0.000E+00 THE FOLLOWING VALUES ARE PROVIDED:
 C[
 C[   Y( 1) =  0.000E+00
 C[   Y( 2) =  1.200E+01
 C[
 C[ REQUIRED PARAMETER:
 C[ ------------------:
 C[ -PRECISION FOR THE IVP =    1.0000000000E-07
 C[ -PRECISION FOR THE BVP =    1.0000000000E-08
 C[ -STEP WIDTH FOR THE IVP=    1.0000000000E-02
 C[ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000
 C[ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000
 C[
 C[ SOLUTION OF THE PROBLEM:
 C[ ------------------------
 C[ THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED
 C[ IN POINT A = 0.000E+00 AS FOLLOWS:
 C[
 C[   Y( 1) =     0.0000000000E+00
 C[   Y( 2) =     9.7229810104E+00
 C[
 C[ THIS REQUIRES     3 NEWTON-ITERATIONS!
 C[
 C[ DETERMINATION COMPLETED AS PLANNED!
*                                                                      *
* (    Je nach Startwertvorgabewert fuer  Y(2)   ( Y(1) = 0          ) *
* (    vorgegeben! )  ist bei dem betrachteten Beispiel einer        ) *
* (    der folgenden Werte fuer Y(2) im Punkt A naeherungsweise      ) *
* (    zu erreichen:                                                 ) *
* (  0 ; +- 9.7229810... ; +- 38.8919241... ; +- 87.5068292... ; ... ) *
*                                                                      *
***********************************************************************/

{
  /* Initialize data; test examples can be exchanged by adjusting     */
  /* the functions dgl(), randbed() and the macro N                   */

  #define N  2

  static REAL yanf0[N]     = { ZERO,
                               (REAL)12.0 };
  static char ystri[N][80] = { "Y(2)",           /* alpha-numeric form*/
                               "-Y(1)**3" };     /* of the DE         */
  REAL        a            = ZERO;
  REAL        b            = ONE;
  REAL        h            = (REAL)0.01;
  REAL        wa           = ZERO;
  REAL        wb           = ZERO;
  long        ifmax        = 20000;
  int         itmax        = 1000;
  REAL        epsawp       = (REAL)1.0e-7;
  REAL        epsrb        = (REAL)1.0e-8;
  REAL        yanf[N];
  int         fehler;
  int         iter;
  int         verfahren;
  int         j;

  #define print(text)  printf(" C[ %s\n", text);

  /* output of the test examples for varying methods        */

  for (verfahren = 1; verfahren <= 5; verfahren++)
  {
    copy_vector(yanf, yanf0, N);                      /* yanf = yanf0 */
    print("");
    print("");
    print("");
    printf(" C[ TEST EXAMPLE (WITH IVP PROGRAM NO %1d):\n",
           verfahren);
    print("=============");
    print("ANALYZED SET OF DIFFERENTIAL EQUATIONS:");
    print("---------------------------------------");
    for (j = 0; j < N; j++)
      printf(" C[          Y'(%2d) = %s\n", j + 1, ystri[j]);
    print("WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)");
    printf(" C[   AT POINT  A = %10.3"LZP"E  :%10.3"LZP"E\n"
           " C[   AT POINT  B = %10.3"LZP"E  :%10.3"LZP"E\n",
           a, wa, b, wb);
    print("");
    printf(" C[ AT LOCATION  A = %10.3"LZP"E THE FOLLOWING "
           "VALUES ARE PROVIDED:\n", a);
    print("");
    for (j = 0; j < N; j++)
      printf(" C[   Y(%2d) = %10.3"LZP"E\n",
             j + 1, yanf[j]);
    print("");
    print("REQUIRED PARAMETER:");
    print("------------------:");
    printf(" C[ -PRECISION FOR THE IVP =%20.10"LZP"E\n"
           " C[ -PRECISION FOR THE BVP =%20.10"LZP"E\n"
           " C[ -STEP WIDTH FOR THE IVP=%20.10"LZP"E\n"
           " C[ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE "
           "IVP'S = %6ld\n"
           " C[ -MAX. NUMBER OF NEWTON-ITERATION STEPS            "
           "      = %6d\n", epsawp, epsrb, h, ifmax, itmax
          );

    fehler = rwp(a, b, h, yanf, N, dgl, randbed, verfahren, epsawp,
                 epsrb, ifmax, itmax, &iter);

    /* report result */

    if (fehler != 0)
    {
      print("");
      switch (fehler)
      {
        case 1:  print("ERRORBOUND(S) TOO SMALL");
                 break;
        case 2:  print("ERROR: B <= A");
                 break;
        case 3:  print("ERROR: STEP SIZE H <= 0");
                 break;
        case 4:  print("ERROR: NUMBER N OF DEQ'S NOT CORRECT");
                 break;
        case 5:  print("ERROR: WRONG IVP-PROGRAM CHOICE");
                 break;
        case 6:  print("IFMAX TO SMALL FOR SOLUTION OF THE "
                       "IVP-PROGRAM");
                 break;
        case 7:  print("AFTER ITMAX STEPS PRECISION WAS NOT REACHED!");
                 break;
        case 8:  print("ERROR: JACOBI-MATRIX IS SINGULAR!");
                 break;
        case 9:  print("LACK of MEMORY!");
                 break;
        default: print("UNKNOWN ERROR");
                 printf(" C[ ERROR NUMBER  = %2d\n", fehler);
                 break;
      }
    }

    else                                              /* no error? */
    {
      /* Output of solution */

      print("");
      print("SOLUTION OF THE PROBLEM:");
      print("------------------------");
      print("THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED");
      printf(" C[ IN POINT A =%10.3"LZP"E AS FOLLOWS:\n", a);
      print("");
      for (j = 0; j < N; j++)
        printf(" C[   Y(%2d) = %20.10"LZP"E\n", j + 1, yanf[j]);
      print("");
      printf(" C[ THIS REQUIRES %5d NEWTON-ITERATIONS!\n", iter);
      print("");
      print("DETERMINATION COMPLETED AS PLANNED!");
    }
  }

  return 0;
}
#endif

/* --------------------------- END m_rwp.c -------------------------- */
