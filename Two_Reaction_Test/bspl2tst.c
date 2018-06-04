#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>      /*  for  umleiten, printf, scanf, fprintf,    */
                        /*       stderr, REAL, LZS, fehler_melden,    */
                        /*       LZP, size_t, readln                  */
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vminit, PMATRIX */
#include <bspline.h>    /*  for  bspline2                             */
#include <zeigkrv2.h>   /*  for  zeigflaeche                          */
#ifdef __TURBOC__
#if defined(__MSDOS__)
#include <alloc.h>      /*  for  coreleft                             */
#if defined(__TINY__) || defined(__SMALL__) || defined(__MEDIUM__)
#define memtyp  unsigned int
#define MT      "u"
#else                   /* memory models Compact, Large or Huge?      */
#define memtyp  unsigned long
#define MT      "lu"
#endif
#elif defined(__TOS__)
#include <ext.h>        /*  for  coreleft                             */
#define memtyp  unsigned long
#define MT      "lu"
#endif
#endif



/*--------------------------------------------------------------------*/

int main
        (
         int  argc,
         char *argv[]
        )

/***********************************************************************
* This is a test program for the function bspline2() from the module   *
* bspline which computes a desired number of points on an open or      *
* closed uniform B spline surface of specified order.                  *
*                                                                      *
* Scope of program:                                                    *
* =================                                                    *
* The program reads and writes from and to stdin and stdout.           *
* The first and second entries of the command line are associated with *
* these files. Error messages and read requests are written to stderr. *
* If the macro INTERAKTIV is not defined before compilation, the input *
* prompts are suppressed.                                              *
*                                                                      *
* Reading of the input is followed by its output. Then bspline2()      *
* computes the desired points on the curve and the results are put     *
* out. For Turbo C or QuickC compilers we also plot the                *
* B spline surface together with the underlying de Boor points.        *
* curve. This plot will not be generated if the third entry in the     *
* command line is "n".                                                 *
*                                                                      *
* Several test data sets are included in the files bspl2tst.ei*. The   *
* correct output should mimick the entries in bsplitst.au* .           *
*                                                                      *
* Structure of input files:                                            *
* =========================                                            *
* k                         order of B spline surface (3<=k<=m,n)      *
* nv                        number of v curve points to be computed    *
* nw                        number of w curve points to be computed    *
* voffen                    open (1) or closed v curves (0)            *
* woffen                    open (1) or closed w curves (0)            *
* m                         number of de Boor points in v direction    *
* n                         number of de Boor points in w direction    *
* d[0][0][0]    ...d[0][0][2]     \                                    *
* d[0][1][0]    ...d[0][1][2]      \  1st de Boor-Vektor               *
* ...                              /                                   *
* d[0][n-1][0]  ...d[0][n-1][2]   /                                    *
* d[1][0][0]    ...d[1][0][2]     \                                    *
* d[1][1][0]    ...d[1][1][2]      \  2nd de Boor-Vektor               *
* ...                              /                                   *
* d[1][n-1][0]  ...d[1][n-1][2]   /                                    *
*     .                                                                *
*     .                                                                *
*     .                                                                *
* d[m-1][0][0]  ...d[m-1][0][2]   \                                    *
* d[m-1][1][0]  ...d[m-1][1][2]    \  mth de Boor-Vektor               *
* ...                              /                                   *
* d[m-1][n-1][0]...d[m-1][n-1][2] /                                    *
***********************************************************************/

{
  REAL ***d,      /* [0..m-1,0..n-1,0..2] matrix with given m*n       */
                  /* de Boor points                                   */
       ***c;      /* [0..nv-1,0..nw-1,0..2] matrix with the computed  */
                  /* points on the B spline surface                   */
  int  voffen,    /* open v curves?                                   */
       woffen,    /* open w curves?                                   */
       k,         /* order of the B spline surface (3<=k<=max(m,n))   */
       nv,        /* number of v curve points to be computed          */
       nw,        /* number of w curve points to be computed          */
       m,         /* number of de Boor points in v direction          */
       n,         /* number of de Boor points in w direction          */
       fehler,    /* error code from umleiten(), bspline2(),          */
                  /* zeigflaeche()                                    */
       i,         /* loop variable                                    */
       j,         /* loop variable                                    */
       l;         /* loop variable                                    */
  void *vmblock;  /* list of dynamically allocated vectors/matrices   */
#if defined(__TURBOC__) && (defined(__MSDOS__) || defined(__TOS__))
  memtyp SpeicherVorher,      /* free heap memory at several places   */
         SpeicherVorGrafik,   /* of main() (necessary to make finding */
         SpeicherNachGrafik,  /* errors in dynamic array easier)      */
         SpeicherNachher;
#endif


#ifdef __TURBOC__
#if defined(__MSDOS__) || defined(__TOS__)
  SpeicherVorher = coreleft();
#endif
#endif

  if ((fehler = umleiten(argc, argv))    /* assign input/output files */
      != 0)                              /* to standard ones          */
    return fehler;  /* 1 or 2 */


  /* ------------------------ read input data ----------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
          "order of B spline surface:                             ");
#endif
  scanf("%d", &k);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr,
          "number of desired v curve points:                      ");
#endif
  scanf("%d", &nv);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr,
          "number of desired w curve points:                      ");
#endif
  scanf("%d", &nw);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr,
          "open (1) or closed (0) v curves?:                      ");
#endif
  scanf("%d", &voffen);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr,
          "open (1) or closed (0) w curves?:                      ");
#endif
  scanf("%d", &woffen);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr,
          "number of de Boor points in v direction:               ");
#endif
  scanf("%d", &m);
  readln();
#ifdef INTERAKTIV
  fprintf(stderr,
          "number of de Boor points in w direction:               ");
#endif
  scanf("%d", &n);
  readln();

  /* -------------- allocate memory for point matrices -------------- */

  vmblock = vminit();                      /* initialize memory block */
  d = (REAL ***)vmalloc(vmblock, PMATRIX, m,  n);         /* allocate */
  c = (REAL ***)vmalloc(vmblock, PMATRIX, nv, nw);        /* memory   */
  if (! vmcomplete(vmblock))                 /* allocations partially */
  {                                          /* unsuccessful?         */
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
#ifdef INTERAKTIV
      fprintf(stderr, "de Boor point (%2d,%2d) "
                      "(x, y and z coordinate):         ", i, j);
#endif
      for (l = 0; l < 3; l++)
        scanf("%"LZS"f", &d[i][j][l]);
      readln();
    }

  /* ----------------------- print input data ----------------------- */

  printf("\n\n"
         "Uniform B spline surface\n"
         "========================\n\n"
         "number m of de Boor points in v direction:    %6d\n"
         "number n of de Boor points in w direction:    %6d\n"
         "order k of B spline surface:                  %6d\n"
         "flag voffen (open v curves?):                 %6d\n"
         "flag woffen (open w curves?):                 %6d\n"
         "number nv of desired v curve points:          %6d\n"
         "number nw of desired w curve points:          %6d\n\n"
         "de Boor points:\n\n"
         "  i  j    d[i,j,0]    d[i,j,1]    d[i,j,2]\n"
         "------------------------------------------\n",
         m, n, k, voffen, woffen, nv, nw
        );

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      printf("%3d%3d", i, j);
      for (l = 0; l < 3; l++)
        printf("%12.5"LZP"f", d[i][j][l]);
      printf("\n");
    }

  /* ------------------- compute B spline surface ------------------- */

  fehler = bspline2(d, m, n, k, voffen, woffen, c, &nv, &nw);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("bspline2(): wrong value for m, k, or nv",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("bspline2(): wrong value for n, k, or nw",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("bspline2(): unknow error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler)
    return 10 + fehler;

  /* --------------- print computed points on surface --------------- */

  printf("\n"
         "computed points on the B spline surface:\n\n");
  printf("  i  j    c[i,j,0]    c[i,j,1]    c[i,j,2]\n"
         "------------------------------------------\n");
  for (i = 0; i < nv; i++)
    for (j = 0; j < nw; j++)
    {
      printf("%3d%3d", i, j);
      for (l = 0; l < 3; l++)
        printf("%12.5"LZP"f", c[i][j][l]);
      printf("\n");
    }

  /* ---------------- plot spline surface if desired ---------------- */

#ifdef __TURBOC__
#if defined(__MSDOS__) || defined(__TOS__)
  SpeicherVorGrafik = coreleft();
#endif
#endif

  fehler = 0;
  if (argc <= 3 || *argv[3] != 'n')                  /* plot desired? */
    fehler = zeigflaeche(c, nv, nw, d, m, n, voffen, woffen, 1);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("zeigflaeche(): nv too small",
                    30 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("zeigflaeche(): nw too small",
                   30 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("zeigflaeche(): m too small",
                    30 + fehler, __FILE__, __LINE__);
      break;
    case 4:
      fehler_melden("zeigflaeche(): n too small",
                    30 + fehler, __FILE__, __LINE__);
      break;
    case 5:
      fehler_melden("zeigflaeche(): nw or n too large",
                    30 + fehler, __FILE__, __LINE__);
      break;
    case 6:
      fehler_melden("zeigflaeche(): graphics error",
                    30 + fehler, __FILE__, __LINE__);
      break;
    case 7:
      fehler_melden("zeigflaeche(): lack of memory",
                    30 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("zeigflaeche(): other error",
                   30 + fehler, __FILE__, __LINE__);
  }


#ifdef __TURBOC__
#if defined(__MSDOS__) || defined(__TOS__)
  SpeicherNachGrafik = coreleft();
#endif
#endif
  vmfree(vmblock);
#ifdef __TURBOC__
#if defined(__MSDOS__) || defined(__TOS__)
  SpeicherNachher = coreleft();

  fprintf(stderr, "\n"
                  "free memory at program start:       %"MT"\n"
                  "free memory before graphics:        %"MT"\n"
                  "free memory after graphics:         %"MT"\n"
                  "free memory at program end:         %"MT"\n"
                  "difference:                         %"MT"\n",
                  SpeicherVorher, SpeicherVorGrafik, SpeicherNachGrafik,
                  SpeicherNachher, SpeicherVorher - SpeicherNachher
         );
#endif
#endif


  return 0;
}
