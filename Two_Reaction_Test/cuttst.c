#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, fprintf, stderr, scanf,   */
                         /*       printf, sprintf, REAL, LZS, LZP,    */
                         /*       ZERO, fehler_melden                 */
#include <vmblock.h>     /*  for  vminit, vmalloc, vmcomplete, VEKTOR */
#include <cuthill.h>     /*  for  cutgaucho                           */



/* ------------------------------------------------------------------ */

int main
        (
         int  argc,
         char *argv[]
        )

/***********************************************************************
* Test program for the function cutgaucho().                           *
*                                                                      *
* Scope:                                                               *
* ======                                                               *
* input a sparse matrix, produce several right hand sides, reduce the  *
* matrix band width via Cuthill-McKee and solve the linear systems     *
* using the Gauss or Cholesky method.                                  *
*                                                                      *
* Form of the input file:                                              *
* =======================                                              *
* gauss        Flag, that asks for Gaussian elimination (1) or the     *
*              Cholesky method (0)                                     *
* n            size of matrix                                          *
* <row  1> \   for every nonzero matrix element, one pair              *
* ...       >  <column index, matrix element> is stored. A zero in the *
* <row  n> /   index position denotes the end of a row. Note that the  *
*              column indices start with 1, which is unusual in C.     *
*              The column indices of one row must be monotonically     *
*              increasing.                                             *
***********************************************************************/

{
#define MAX_NV  185  /* maximal number of tupels allowed in input     */

#define NRS       5  /* number of right hand sides                    */

  int     gauss;      /* flag: Gauss (TRUE) or Cholesky (FALSE)       */
  REAL    v[MAX_NV],  /* Vector for read in matrix elements           */
          *rs,        /* [0..n*NRS-1] vector for NRS r. h. sides      */
          *x;         /* [0..n*NRS-1] vector for NRS solutions        */
  int     n,          /* size of matrix, as read in                   */
          nv,         /* number of <index, element> tupels read in    */
          ic[MAX_NV], /* vector of column indices in v                */
          m,          /* Half band width of the condensed matrix      */
          fehler,     /* return value of umleiten(), scanf() and      */
                      /* cutgaucho()                                  */
          i,          /* Loop variables                               */
          j,          /*                                              */
          k;          /*                                              */
  void    *vmblock;   /* List of dynamically allocated vectors and    */
                      /* matrices                                     */
  char    meldung[200];             /* aux buffer for error messages  */


  if ((fehler = umleiten(argc, argv))   /* assign given in/output     */
      != 0)                             /* files to standard in/output*/
    return fehler;  /* 1 or 2 */        /* files                      */


  /* --------------------- Set up input ----------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\nMethod of solution (0 = Cholesky, 1 = Gauss): ");
#endif
  fehler = scanf("%d", (int *)&gauss);
  if (fehler != 1)                                    /* Input error? */
  {
    fehler_melden("wrong character or premature end of read in",
                  0, __FILE__, __LINE__);
    return 3;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "\nsize of matrix:  ");
#endif
  fehler = scanf("%d", &n);
  if (fehler != 1)                                  /* Input error?   */
  {
    fehler_melden("wrong character or premature end of input",
                  0, __FILE__, __LINE__);
    return 3;
  }


#ifdef INTERAKTIV
  fprintf(stderr, "put in nonzero matrix elements row wise"
                  "\n(Indexing from 1, zero index: end of row):"
                  "\n");
#endif

  for (i = nv = 0; i < n; )                 /* loop over matrix rows  */
  {
    if (nv >= MAX_NV)                       /* too many input pairs?  */
    {                                       /* interrupt program      */
      fehler_melden("too many input pairs", 0, __FILE__, __LINE__);
      return 3;
    }
#ifdef INTERAKTIV
    fprintf(stderr, "row %4d: column:  ", i + 1);
#endif
    fehler = scanf("%d", ic + nv);
    if (fehler != 1)                                /* input error?   */
    {
      fehler_melden("wrong character or premature end of data",
                    0, __FILE__, __LINE__);
      return 3;
    }
#ifdef INTERAKTIV
    fprintf(stderr, "            Element: ");
#endif
    fehler = scanf("%"LZS"f", v + nv);
    if (fehler != 1)                                /* input error ?  */
    {
      fehler_melden("wrong character or premture end in input",
                    0, __FILE__, __LINE__);
      return 3;
    }
    if (ic[nv] == 0)                       /* End of row?             */
      i++;                                 /* next row                */
    nv++;                                  /* next pair               */
  }


  /* ------------ dynamical allocation of vectors and matrices ------ */

  vmblock = vminit();
  rs = (REAL *) vmalloc(vmblock, VEKTOR, n * NRS, 0);
  x  = (REAL *) vmalloc(vmblock, VEKTOR, n * NRS, 0);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }


  /* -------------- Print out input for control purposes ------------ */

  printf("\n"
         "Solve a linear system with a sparse symmetric matrix\n"
         "====================================================\n"
         "using Cuthill-KcKee and %s\n"
         "========================%s\n"
         "\n\n\n",
         gauss ? "Gauss" : "Cholesky", gauss ? "=====" : "========"
        );

  printf("Input data:\n"
         "-----------\n\n"
         "size of matrix = %d\n\n"
         "read in matrix elements:\n", n);

  for (i = 0; i < nv; i++)                   /* Matrix as read in     */
    printf("%3i%12.2"LZP"e\n", ic[i], v[i]);

#ifdef DEBUG
#define FORMAT  "%3.2"LZP"g"
#else
#define FORMAT  "%12.2"LZP"e"
#endif
  printf("\nMatrix:\n");
  for (i = k = 0; i < n; i++, k++)    /* Matrix in standard form      */
  {
    for (j = 1; j <= n; j++)
      if (j == ic[k])
      {
        printf(FORMAT, v[k]);
        k++;
      }
      else
        printf(FORMAT, ZERO);
    printf("\n");
  }


  /* --------------- construct NRS right hand sides ----------------- */
  /* ------------- (kth solution vector = (k+1,...,k+n) ------------- */

  if (n > 0)
    for (k = 0; k < NRS; k++, rs += n)
      for (j = i = 0, rs[0] = ZERO; j < nv; j++)
      {
        if (ic[j] == 0)
        {
          i++;
          if (i < n)
            rs[i] = ZERO;
        }
        else
          rs[i] += v[j] * (ic[j] + k);
      }
  rs -= NRS * n;                /* reset pointer rs to original value */


  /* ---------------- solve linear system --------------------------- */

  fehler = cutgaucho((boolean)gauss, n, nv, ic, v, NRS, rs, x, &m);


  /* -------------------- print results ----------------------------- */

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      sprintf(meldung, "cutgaucho(): n <= 0 or nv <= 0 or "
                       "nrs <= 0 or m < 0");
      break;
    case 2:
      sprintf(meldung, "cutgaucho(): lack of memory");
      break;
    case 3:
      sprintf(meldung, "cutgaucho(): Matrix numerically not%s regular",
                       gauss ? "" : " strongly");
      break;
    case 5:
      sprintf(meldung, "cutgaucho(): n doea not match the zero indices "
                       "in ic");
      break;
    case 6:
      sprintf(meldung, "cutgaucho(): "
                       "read in column index invalid");
      break;
    case 7:
      sprintf(meldung, "cutgaucho(): Matrix not %ssymmetric",
                       gauss ? "weakly " : "");
      break;
    default:
      sprintf(meldung, "cutgaucho(): unknown error");
  }

  if (fehler)
  {
    fehler_melden(meldung, 4 + fehler, __FILE__, __LINE__);
    return 4 + fehler;
  }

  printf("\n\n\n"
         "Output data:\n"
         "------------\n\n\n"
         "Half band width after Cuthill-McKee = %d\n", m);

  for (k = 0; k < NRS; k++) /* print computed and compare with exact  */
  {                         /* solutions                              */
    printf("\n\nsolution %3d:\n\n"
           "  i    right hand side     x[i]"
           "      exact solution  difference\n"
           "-------------------------------"
           "--------------------------------\n", k);
    for (i = 0; i < n; i++)
      printf("%3d%16.8"LZP"f%16.8"LZP"f%16.8"LZP"f%11.3"LZP"f\n",
             i, rs[i], x[i], (REAL)(i + k + 1), (i + k + 1) - x[i]);
    rs += n;                                 /* next right hand side  */
    x  += n;                                 /* next solution         */
  }


  return 0;
}
