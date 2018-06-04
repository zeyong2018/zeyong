#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>       /*  for  umleiten, printf, scanf, fprintf,   */
                         /*       stderr, REAL, LZS, fehler_melden,   */
                         /*       LZP                                 */
#include <vmblock.h>     /*  for  vmalloc, vmcomplete, vminit, MATRIX */
#include <bspline.h>     /*  for  bspline                             */
#include <zeigkrv2.h>    /*  for  zeigkrv2                            */



/*--------------------------------------------------------------------*/

int main(int  argc,
         char *argv[]
        )

/***********************************************************************
* This is a test program for the function bspline() from the module    *
* bspline which computes a a desired number of points on a uniform     *
* open or closed B spline curve of specified order.                    *
*                                                                      *
* Scope of program:                                                    *
* ==================                                                   *
* The program reads and writes from and to stdin and stdout.           *
* The first and second entries of the command line are associated with *
* these files. Error messages and read requests are stored in stderr.  *
* If the macro INTERAKTIV is not defined before compilation, the input *
* prompts are suppressed.                                              *
*                                                                      *
* Reading of the input is followed by its output. Then bspline()       *
* computes the desired points on the curve and the results are put     *
* out. For Turbo C or QuickC compilers we also plot the B spline curve *
* together with the underlying de Boor points in case of a planar      *
* curve. This plot will not be generated if the third entry in the     *
* command line is "n".                                                 *
*                                                                      *
* Several test data sets are included in the files bsplitst.ei*. The   *
* correct output should mimick the entries in bsplitst.au* .           *
*                                                                      *
* Structure of input files:                                            *
* =========================                                            *
* k                          Order of  B spline curve   (3 <= k <= n)  *
* m                          Dimension of space for de Boor points     *
*                            (m >= 2)                                  *
* nc                         number of computed points on the curve    *
* n                          number of de Boor points                  *
* d[0][0],  ...,d[0][m-1]    1st de Boor point                         *
* d[1][0],  ...,d[1][m-1]    2nd de Boor point                         *
*   ...           ...        ...                                       *
* d[n-1][0],...,d[n-1][m-1]  nth de Boor point                         *
***********************************************************************/

{
  REAL **d,        /* [0..n-1,0..m-1] matrix with de Boor points      */
       **c;        /* [0..nc-1,0..m-1] matrix with points on curve    */
  int  k,          /* Order of B spline curve (3 <= k <= n)           */
       m,          /* Dimension of space for de Boor points           */
                   /* (m >= 2)                                        */
       nc,         /* number of desired points on curve               */
       n,          /* number of de Boor points                        */
       fehler,     /* error code from bspline()                       */
       i,          /* loop counter                                    */
       j;          /* loop counter                                    */
  void *vmblock;   /* List of dynamically allocated vectors/matrices  */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard ones           */
    return fehler;  /* 1 or 2 */

  /* ---------------------- read input data ------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n"
                  "Order of B spline curve:             ");
#endif
  scanf("%d", &k);
#ifdef INTERAKTIV
  fprintf(stderr, "space dimension for points of curve: ");
#endif
  scanf("%d", &m);
#ifdef INTERAKTIV
  fprintf(stderr, "number of desired points:            ");
#endif
  scanf("%d", &nc);
#ifdef INTERAKTIV
  fprintf(stderr, "number of de Boor points:            ");
#endif
  scanf("%d", &n);

  vmblock = vminit();                 /* initialize storage buffers   */
  d = (REAL **)vmalloc(vmblock, MATRIX, n,  m);
  c = (REAL **)vmalloc(vmblock, MATRIX, nc, m);
  if (! vmcomplete(vmblock))
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 3;
  }

  for (i = 0; i < n; i++)
  {
#ifdef INTERAKTIV
    fprintf(stderr, "de Boor point %2d (%1d coordinates):        ",
                    i + 1, m);
#endif
    for (j = 0; j < m; j++)
      scanf("%"LZS"f", &d[i][j]);
  }


  /* --------------------- print input data ------------------------- */

  printf("\n\n"
         "Open uniform B spline curve\n"
         "===========================\n\n\n"
         "Order k of the curve:                 %4d\n"
         "Space dimension m of the points:      %4d\n"
         "Number nc of desired points on curve: %4d\n"
         "Number n of de Boor points:           %4d\n\n"
         "de Boor points:\n\n"
         "  i     d[i][0] ... d[i][m-1]\n"
         "---------------------------------------\n",
         k, m, nc, n
        );
  for (i = 0; i < n; i++)
  {
    printf("%3d", i);
    for (j = 0; j < m; j++)
      printf("%12.5"LZP"f", d[i][j]);
    printf("\n");
  }


  /* ------------------- Compute  B spline curve  ------------------- */

  fehler = bspline(d, n, k, m, 1, c, &nc);

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("bspline(): n < 3  or  k < 3  or  k > n  "
                    "or  nc < 2*(n-k+1)+1)",
                    10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("bspline(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
      break;
    default:
      fehler_melden("bspline(): unknown error",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ------------- print out desired points ------------------------ */

  printf("\n\n"
         "Computed points on the curve:\n\n"
         "  i     c[i][0] ... c[i][m-1]\n"
         "---------------------------------------\n"
        );
  for (i = 0; i < nc; i++)
  {
    printf("%3d", i);
    for (j = 0; j < m; j++)
      printf("%12.5"LZP"f", c[i][j]);
    printf("\n");
  }


  /* ------------ plot spline if desired and if curve is planar  ---- */

  if (m == 2 &&                         /* planar curve ?             */
      (argc <= 3 || *argv[3] != 'n')    /* plot desired ?             */
     )
  {
    fehler = zeigkrv2(nc, c, n, d);
    switch (fehler)
    {
      case 0:
        break;
      case 3:
        fehler_melden("zeigkrv2(): lack of memory",
                      30 + fehler, __FILE__, __LINE__);
        break;
      case 4:
        fehler_melden("zeigkrv2(): Graphics error",
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
