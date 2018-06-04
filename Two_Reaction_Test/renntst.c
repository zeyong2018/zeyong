#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>      /*  for  umleiten, abl_mat1, printf, scanf,   */
                        /*       fprintf, stderr, NULL, REAL, LZS,    */
                        /*       ZERO, fehler_melden, LZP             */
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vminit, VEKTOR, */
                        /*       MATRIX                               */
#include <subsplin.h>   /*  for  renner, rennwert, renntab            */
#include <zeigkrv2.h>   /*  for  zeigkrv2                             */



/* ------------------------------------------------------------------ */

int main
        (
         int  argc,
         char *argv[]
        )

/***********************************************************************
* Test program for renner() from the module subsplin to compute        *
* interpolating parametric Renner subsplines in R2 or R3               *
*                                                                      *
* Scope of program:                                                    *
* =================                                                    *
* The program reads from stdin and writes to stdout.                   *
* The first and second entries of the command line are associated with *
* these files. Error messages and input prompts are written to stderr. *
* If the macro INTERAKTIV is not defined before compilation, the input *
* prompts are suppressed.                                              *
*                                                                      *
* Reading of the input is followed by its output. Then renner()        *
* computes the spline coefficients, and finally rennwert() computes    *
* some functional values and derivatives of the subspline function,    *
* and renntab() tabulates a part of it. Each of these three sections - *
* computation of coefficients, evaluation and tabulation - is followed *
* by the output of the results, if no error occurred.                  *
* For Turbo C or QuickC compilers we also plot the spline curve        *
* together with the underlying nodes (in case of a curve in R3,        *
* however, without the z component). This plot will not be generated   *
* if the third entry in the command line is "n".                       *
*                                                                      *
* Several test data sets are provided in renntst.ei*. The output       *
* should be identical to the contents of renntst.au* .                 *
*                                                                      *
* Form of input file:                                                  *
* ===================                                                  *
* rund                         Parameter for rounding:                 *
*                                0 < rund < 1: round corners,          *
*                                for other values of rund: no rounding *
* dim                          dimension of spaces where the nodes lie *
*                                (2 <= dim <= 3)                       *
* maxtab                       maximal number of curve points to be    *
*                                tabulated                             *
* n                            number of nodes                         *
* P[0][0],  ...,P[0][dim-1]    1st node                                *
* P[1][0],  ...,P[1][dim-1]    2nd note                                *
*   ...           ...                                                  *
* P[n-1][0],...,P[n-1][dim-1]  nth node                                *
***********************************************************************/

{
  int      n,         /* number of nodes                              */
           n0,        /* n + [n / 2]                                  */
           nm1,       /* n - 1                                        */
           dim,       /* space dimension (2 or 3)                     */
           fehler,    /* error code from umleiten(), renner(),        */
                      /* renntab() resp. zeigkrv2()                   */
           maxtab,    /* maximal number of curve points to be         */
                      /* tabulated                                    */
           lentab,    /* actual length of table of curve points       */
           i,         /* loop variable                                */
           j;         /* loop variable                                */
                      /* The vectors P, T, boglg, b, c, d must offer  */
                      /* room for maximally n+[n/2] entries,          */
                      /* comprised of the original n nodes and at     */
                      /* most additional [n/2] points from rounding   */
                      /* corners.                                     */
  REAL     **sptab,   /* [0..maxtab-1,0..dim-1] vector                */
           **P,       /* [0..n0-1,0..dim-1] vector                    */
           *T,        /* [0..n0-2] vector                             */
           *boglg,    /* [0..n0-1] vector                             */
           **b,       /* [0..n0-2,0..dim-1] vector                    */
           **c,       /* [0..n0-2,0..dim-1] vector                    */
           **d,       /* [0..n0-2,0..dim-1] vector                    */
           rund,      /* degree of rounding corners                   */
           b0,        /* loop variable when calling rennwert()        */
           s[3],      /* curve deliverd by rennwert()                 */
           banf,      /* left border of tabulation interval           */
           bend,      /* right border of tabulation interval          */
           bstep,     /* step size for tabulation                     */
           ds[4][3];  /* 0th - 3rd derivative of the spline at b0     */
  void     *vmblock;  /* List of dynamic allocations                  */


  if ((fehler = umleiten(argc, argv))   /* assign in/output files to  */
      != 0)                             /* standard files             */
    return fehler;  /* 1 or 2 */


  /* ---------------------- read input   ---------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "rounding corners (0 < rund < 1): ");
#endif
  scanf("%"LZS"f", &rund);
#ifdef INTERAKTIV
  fprintf(stderr, "space dimension of curve:        ");
#endif
  scanf("%d", &dim);
#ifdef INTERAKTIV
  fprintf(stderr, "size of value table:             ");
#endif
  scanf("%d", &maxtab);
#ifdef INTERAKTIV
  fprintf(stderr, "Number of nodes:                 ");
#endif
  scanf("%d", &n);

  n0  = n + n / 2;
  nm1 = n - 1;

  vmblock = vminit();                     /* initialize storage block */
  P     = (REAL **)vmalloc(vmblock, MATRIX, n0,     dim);
  T     = (REAL *) vmalloc(vmblock, VEKTOR, n0 - 1, 0);
  boglg = (REAL *) vmalloc(vmblock, VEKTOR, n0,     0);
  b     = (REAL **)vmalloc(vmblock, MATRIX, n0 - 1, dim);
  c     = (REAL **)vmalloc(vmblock, MATRIX, n0 - 1, dim);
  d     = (REAL **)vmalloc(vmblock, MATRIX, n0 - 1, dim);
  sptab = (REAL **)vmalloc(vmblock, MATRIX, maxtab, dim);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "node        (x[%2d],y[%2d]):       ",
                    i, i);
#endif
    for (j = 0; j < dim; j++)
      scanf("%"LZS"f", &P[i][j]);
  }


  /* --------------------- print out input data  -------------------- */

  printf("interpolating Renner subspline\n"
         "==============================\n\n\n"
         "nodes:\n\n"
         "  i     x[i]        y[i]%s\n"
         "---------------------------%s\n",
         (dim == 3) ? "        z[i]" : "",
         (dim == 3) ? "------------" : ""
        );
  for (i = 0; i < n; i++)
  {
    printf("%3d", i);
    for (j = 0; j < dim; j++)
      printf("%12.5"LZP"f", P[i][j]);
    printf("\n");
  }
  printf("\n");
  if (rund == ZERO)
    printf("Corners are not rounded.\n");
  else
    printf("rounding parameter:    %7.3"LZP"f\n", rund);


  /* ------------------ compute subspline --------------------------- */

  fehler = renner(&nm1, n0 - 1, dim, P, rund, T, b, c, d);
  n = nm1 + 1;

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("renner(): less than 5 nodes",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("renner(): two successive nodes coincide",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 3:
      fehler_melden("renner(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 4:
      fehler_melden("renner(): lack of room for rounding nodes",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 5:
      fehler_melden("renner(): wrong space dimension",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("renner(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* --------------- accumulate parameter values      --------------- */

  for (boglg[0] = ZERO, i = 1; i < n; i++)
    boglg[i] = boglg[i - 1] + T[i - 1];


  /* ---------------- print subspline coefficients    --------------- */

  if (rund > ZERO)
  {
    printf("\n\nnodes with corner rounding nodes included:\n\n"
            "  i     x[i]        y[i]%s\n"
            "---------------------------%s\n",
            (dim == 3) ? "        z[i]" : "",
            (dim == 3) ? "------------" : ""
          );
    for (i = 0; i < n; i++)
    {
      printf("%3d", i);
      for (j = 0; j < dim; j++)
        printf("%12.5"LZP"f", P[i][j]);
      printf("\n");
    }
  }

  for (j = 0; j < dim; j++)
  {
    printf("\n\n"
           "Coefficients of the subspline component S%c:\n\n"
           "  i   T[i](accumul.)  a%c[i]          b%c[i]"
           "          c%c[i]          d%c[i]\n"
           "------------------------------------------"
           "---------------------------------\n",
           'X' + j, 'x' + j, 'x' + j, 'x' + j, 'x' + j
          );

    for (i = 0; i < n - 1; i++)
      printf("%3d%12.5"LZP"f%15.7"LZP"f%15.7"LZP"f%15.7"LZP"f"
             "%15.7"LZP"f\n",
             i, boglg[i], P[i][j], b[i][j], c[i][j], d[i][j]);
  }


  /* -------- compute several function and derivative values -------- */

  printf("\n\n"
         "Example using pspwert():\n"
         "Compute several function and derivative values"
         "\n\n"
         "       "
        );
  if (dim == 2)
    printf("                        "
           "         .           .  "
           "         ..          .. "
           "         ...         ...");
  else
    printf("                                    "
           "         .           .           .  "
           "         ..          ..          .. "
           "         ...         ...         ...");
  printf("\n"
         "T (accumul.)    ");
  if (dim == 2)
    printf("SX(t)       SY(t)       "
           "SX(t)       SY(t)       "
           "SX(t)       SY(t)       "
           "SX(t)       SY(t)\n"
           "----------------"
           "------------------------"
           "------------------------"
           "------------------------"
           "------------------\n");
  else
    printf("SX(t)       SY(t)       SZ(t)       "
           "SX(t)       SY(t)       SZ(t)       "
           "SX(t)       SY(t)       SZ(t)       "
           "SX(t)       SY(t)       SZ(t)\n"
           "----------------"
           "------------------------------------"
           "------------------------------------"
           "------------------------------------"
           "------------------------------\n");

  bstep = (boglg[n - 1] - boglg[0]) / (REAL)15.0;
  for (b0 = boglg[0]; b0 <= boglg[n - 1]; b0 += bstep)
  {
    rennwert(n - 1, dim, b0, boglg, P, b, c, d, s, ds);
    printf("%10.5"LZP"f", b0);
    for (i = 0; i < 4; i++)
      for (j = 0; j < dim; j++)
        printf("%12.5"LZP"f", ds[i][j]);
    printf("\n");
  }


  /* - make a table of values for the spline from  boglg[0]        -  */
  /* - to boglg[n-1] with minimal step size                        -  */

  bstep = boglg[n - 1] - boglg[0];
  banf  = boglg[0]     - (REAL)0.000 * bstep;
  bend  = boglg[n - 1] + (REAL)0.000 * bstep;
  bstep = (bend - banf) / (maxtab - 1 - n);

  fehler = renntab(n - 1, banf, bend, bstep, maxtab - 1,
                   boglg, dim, P, b, c, d, sptab, &lentab);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("renntab(): tanf > tend",
                    20 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("renntab(): step size <= 0",
                    20 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("renntab(): unknown error",
                    20 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 20 + fehler;

  printf("\n\n"
         "Example using renntab():\n"
         "Table of values from arclength =%10.5"LZP"f to %10.5"LZP"f, "
         "Step size =%10.5"LZP"f\n\n",
         banf, bend, bstep);

  printf("  i        xtab[i]        ytab[i]%s\n"
         "---------------------------------%s\n",
         (dim == 3) ? "        ztab[i]" : "",
         (dim == 3) ? "---------------" : ""
        );
  for (i = 0; i <= lentab; i++)
  {
    printf("%3d", i);
    for (j = 0; j < dim; j++)
      printf("%15.7"LZP"f", sptab[i][j]);
    printf("\n");
  }


  /* ------------ plot spline if desired (attention:      ----------- */
  /* ------------ curves in R3 do not show their z part.) ----------- */

  if (argc <= 3 || *argv[3] != 'n')     /* plot not suppressed?       */
  {
    fehler = zeigkrv2(lentab + 1, sptab, n, P);
    switch (fehler)
    {
      case 0:
        break;
      case 3:
        fehler_melden("zeigkrv2(): lack of memory",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 4:
        fehler_melden("zeigkrv2(): graphics error",
                      30 + fehler, __FILE__, __LINE__);
        break;
      default:
        fehler_melden("zeigkrv2(): other error",
                      30 + fehler, __FILE__, __LINE__);
    }
    if (fehler != 0)
      return 30 + fehler;
  }


  return 0;
}
