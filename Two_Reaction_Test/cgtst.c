#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
* This is a test program for the function cg_verfahren() from the      *
* module cg to solve a linear system                                   *
*                       A * X  =  Y                                    *
* for a symmetric positive definite matrix A using the conjugate       *
* gradient method.                                                     *
*                                                                      *
* Scope of the program :                                               *
* ======================                                               *
* The first parameter on the command line is treated as the input file *
* the second one is the output file.                                   *
* If theis second one is void, all output is played to the screen.     *
* If the first file is void also, the input must be made from the      *
* keyboard.                                                            *
* The correct solution for all inputs is a vector of all ones.         *
* The input file has the format :                                      *
*                                                                      *
* 4                     size N of matrix                               *
* 5.0 -1.0 -1.0 -1.0    the upper triangle (with diagonal)             *
*      5.0 -1.0 -1.0    of the positive definite matrix A              *
*           5.0 -1.0                                                   *
*                5.0                                                   *
*                                                                      *
* The text to the right of the numerical entries is ignored.           *
* Some test matrices are available in this format in cgtst.ei*.        *
************************************************************************
* Programing language: ANSI C                                          *
* Compiler:            Turbo C 2.0                                     *
* Computer:            IBM PS/2 70 with 80387                          *
* Author:              Juergen Dietel, Computer Center, RWTH Aachen    *
* Date:                7.31.1992                                       *
***********************************************************************/

#include <basis.h>     /*  for  REAL, scanf, readln, fprintf, stderr, */
                       /*       ZERO, LZS, printf, LZP                */
#include <vmblock.h>   /*  for  vminit, vmalloc, vmcomplete, vmfree,  */
                       /*       VEKTOR, MATRIX                        */
#include <cg.h>        /*  for  cg_verfahren                          */



int main(int argc, char *argv[])

{
  REAL **a,       /* the upper riangle of a positive definite real    */
                  /* matrix.                                          */
       *y,        /* right hand side of the linear system             */
       *x,        /* solution vector of the system                    */
       sum;       /* aux sum for right hand side                      */
  int  n,         /* size of matrix A                                 */
       i,         /* row index                                        */
       j,         /* column index                                     */
       fehler;    /* return value of cg_verfahren()                   */
  void *vmblock;  /* List of dynamically allocated vectors and        */
                  /* matrices                                         */


  if ((fehler = umleiten(argc, argv))   /* assign a given input or    */
      != 0)                             /* output file to the         */
    return fehler;  /* 1 or 2 */        /* standard in/output file    */

  /* ----------------------------- Set inputs ----------------------- */

  scanf("%d", &n);      /* read size of A from standard input file    */
  readln();

  vmblock = vminit();
  a = (REAL **)vmalloc(vmblock, MATRIX, n, n);  /* allocate buffers   */
  y = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);  /* for matrx, right   */
  x = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);  /* hand side and      */
                                                /* solution vector    */

  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

  for (i = 0; i < n; i++)              /* assign test case values to  */
  {                                    /* A and right hand side y     */
    sum = ZERO;
    for (j = i; j < n; j++)            /* We set up the right hand    */
    {                                  /* side y as y = A * vector of */
      scanf("%"LZS"f", &a[i][j]);      /* ones.                       */
      sum += a[i][j];
    }
    readln();
    for (j = 0; j < i; j++)
      sum += a[j][i];
    y[i] = sum;
  }


  /* -------------- Print input data as a safeguard ---------------- */

  printf("\nTest data (Matrix and right hand side):\n");
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < i; j++)
      printf(" %9.5"LZP"f", a[j][i]);
    for (j = i; j < n; j++)
      printf(" %9.5"LZP"f", a[i][j]);
    printf("   %9.5"LZP"f\n", y[i]);
  }


  /* ---------------- solve linear system --------------------------- */

  fehler = cg_verfahren(n, a, y, x);           /* perform  CG method  */


  /* -------------------- Print output results ---------------------- */

  if (fehler)                                  /* Error in CG method? */
  {
    fehler_melden("cg_verfahren()", 4 + fehler, __FILE__, __LINE__);
    return 4 + fehler;
  }

  printf("\nSolution for the linear system:\n");
  for (i = 0; i < n; i++)
    printf(" %"LZP"f", x[i]);

  return 0;
}
