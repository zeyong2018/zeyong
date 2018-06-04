#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*--------------------------------------------------------------------*/
/* Test program for the adaptive SOR method                           */
/*--------------------------------------------------------------------*/

#include <basis.h>      /*  for  REAL, umleiten, ZERO, ONE, MACH_EPS  */
#include <vmblock.h>    /*  for  vminit, vmalloc, VEKTOR, MATRIX,     */
                        /*       vmcomplete, WriteHead,               */
                        /*       LogError, ReadMat, WriteMat,         */
                        /*       CopyMat, SetVec, WriteVec, WriteEnd  */
#include <fadsor.h>     /*  for  adsor                                */



int main (int  argc,
          char *argv[]
         )

{
  REAL **a,
       **c,
       *b,
       *x,
       *residu,
       omega;
  int  n,
       j,
       rc,
       krit,
       iter,
       fehler;
  void *vmblock;


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to the standard in/output  */
    return fehler;  /* 1 or 2 */        /* files                      */


  WriteHead("adaptive SOR method");

  if (scanf("%d ", &n) <= 0)
  {
    LogError("error reading with scanf()", 0,  __FILE__, __LINE__);
    return 3;
  }

  if (n < 1)
  {
    LogError("n < 1", 0,  __FILE__, __LINE__);
    return 4;
  }

  vmblock = vminit();
  a      = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  c      = (REAL **)vmalloc(vmblock, MATRIX, n, n);
  b      = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
  x      = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
  residu = (REAL *) vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    LogError("lack of memory", 0, __FILE__, __LINE__);
    return 5;
  }

  if (ReadMat(n, n, a) != 0)
  {
    LogError("error reading with ReadMat()", 0,  __FILE__, __LINE__);
    return 6;
  }

  printf("size of matrix: %d\n", n);
  printf("Input matrix:\n");

  WriteMat(n, n, a);
  CopyMat(n, n, a, c);

  printf("\ntransposed inverse matrix:\n\n");

  omega   = (REAL)1.5;
  krit    = 0;

  for (j = 0; j < n; j++)
  {
    SetVec(n, x, ZERO);
    SetVec(n, b, ZERO);
    b[j] = ONE;
    CopyMat(n, n, c, a);

    rc = adsor(krit, n, a, b, &omega, x, residu, &iter,
               5, (REAL)2048.0 * MACH_EPS, 300, 0);
    if (rc != 0)
    {
      LogError("adsor()", rc,  __FILE__, __LINE__);
      return 7;
    }
    printf("%5d Iterations:  ", iter);
    WriteVec(n, x);
  }

  printf("\nResidual vector:\n");
  WriteVec(n, residu);

  WriteEnd();


  return 0;
}
