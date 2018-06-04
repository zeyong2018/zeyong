#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*--------------------------------------------------------------------*/
/* Test environment for symmetric tridiagonal matrices trdiasy()      */
/*--------------------------------------------------------------------*/

#include <basis.h>     /*  for  REAL, freopen, stdout, NULL, fprintf, */
                       /*       stderr, ONE, TEN, printf, ZERO,       */
                       /*       FORMAT_126LF, FORMAT_LE, WriteHead,   */
                       /*       LogError, SetVec,                     */
                       /*       WriteVec, WriteEnd                    */
#include <vmblock.h>   /*  for  vminit, vmalloc, VEKTOR, vmcomplete   */
#include <ftrdiasy.h>  /*  for  trdiasy                               */



#define DIM  10                             /* size of test matrix    */



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
       *rs,
       determ;
  void *vmblock;


  /* --- assign the input file to the standard input file ----------- */

  if (argc >= 2)                           /* at least one entry ?    */
    if (freopen(argv[1], "w", stdout) == NULL)        /* open input   */
    {                                                 /* file         */
      fprintf(stderr, "error in opening input file %s!\n", argv[1]);
      return 1;
    }


  WriteHead("tridiagonal symmetric matrix (using trdiasy())");

  n = DIM;

  vmblock = vminit();
  diag = (REAL *)vmalloc(vmblock, VEKTOR, n,     0);
  oben = (REAL *)vmalloc(vmblock, VEKTOR, n - 1, 0);
  rs   = (REAL *)vmalloc(vmblock, VEKTOR, n,     0);
  if (! vmcomplete(vmblock))
  {
    LogError("lack of memory space", 0, __FILE__, __LINE__);
    return 1;
  }

  for (i = 0; i < n; i++)
    oben[i] = ONE,
    diag[i] = TEN;

  printf("Condensed input matrix:\n\n");
  for (i = 0; i < n; i++)
    printf(FORMAT_126LF FORMAT_126LF"\n", diag[i], oben[i]);

  printf("\ninverse matrix:\n\n");
  modus = 0;
  for (j = 0; j < n; j++)
  {
    SetVec(n, rs, ZERO);
    rs[j]  = ONE;
    fehler = trdiasy(modus, n, diag, oben, rs);
    if (fehler)
    {
      LogError("trdiasy()", fehler,  __FILE__, __LINE__);
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
