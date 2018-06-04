#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*--------------------------------------------------------------------*/
/* Testing environment for cyclic tridiagonal symmetric matrices      */
/* (zytrdsy())                                                        */
/*--------------------------------------------------------------------*/

#include <basis.h>     /*  for  REAL, freopen, stdout, NULL, fprintf, */
                       /*       stderr, ONE, FOUR, printf, ZERO,      */
                       /*       FORMAT_126LF, FORMAT_LE, WriteHead,   */
                       /*       LogError, SetVec,                     */
                       /*       WriteVec, WriteEnd                    */
#include <vmblock.h>   /*  for  vminit, vmalloc, VEKTOR, vmcomplete   */
#include <fzytrdsy.h>  /*  for  zytrdsy                               */



#define DIM 4                               /* size of test matrix    */



/* ------------------------------------------------------------------ */

int main(int  argc,
         char *argv[]
        )

{
  int  n,
       modus,
       fehler,
       i,
       j;
  REAL *diag,
       *oben,
       *rechts,
       *rs,
       determ;
  void *vmblock;


  /* --- assign the input file to standard input -------------------- */

  if (argc >= 2)                           /* at least one entry ?    */
    if (freopen(argv[1], "w", stdout) == NULL)        /* open output  */
    {                                                 /* file         */
      fprintf(stderr, "error in opening file %s!\n", argv[1]);
      return 1;
    }


  WriteHead("cyclic tridiagonal symmetric matrix (zytrdsy())");

  n = DIM;

  vmblock = vminit();
  diag   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  oben   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  rechts = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  rs     = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    LogError("lack of storage", 0, __FILE__, __LINE__);
    return 1;
  }

  for (i = 0; i < n; i++)
    oben[i] = -ONE,
    diag[i] = FOUR;

  printf("Condensed input matrix:\n\n");
  for (i = 0; i < n; i++)
    printf(FORMAT_126LF FORMAT_126LF"\n", diag[i], oben[i]);

  printf("\ninverse matrix:\n\n");
  modus = 0;
  for (j = 0; j < n; j++)
  {
    SetVec(n, rs, ZERO);
    rs[j] = ONE;
    fehler = zytrdsy(modus, n, diag, oben, rechts, rs);
    if (fehler)
    {
      LogError("zytrdsy()", fehler,  __FILE__, __LINE__);
      return 1;
    }

    WriteVec(n, rs);
    modus = 2;
  }

  for (determ = ONE, i = 0; i < n; i++)
    determ *= diag[i];

  printf("\nDeterminant = "FORMAT_LE, determ);

  WriteEnd();


  return 0;
}
