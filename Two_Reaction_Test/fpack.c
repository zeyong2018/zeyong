#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"



/*.FE{C 4.13}{Linear Systems with Band Matrices}
             {Linear Systems with Band Matrices}*/

/* ------------------------- MODULE fpack.c ------------------------- */

#include <basis.h>
#include <u_proto.h>


int pack                /* condense a row ............................*/
/*.IX{pack}*/
         (
          int     n,              /* size of matrix ..................*/
          int     ld,             /* number of lower co-diagonals ....*/
          int     ud,             /* number of upper co-diagonals ....*/
          int     no,             /* row index .......................*/
          REAL    row[],          /* original row ....................*/
          REAL    prow[]          /* condensed row ...................*/
         )
/*====================================================================*
 *                                                                    *
 *  pack converts the no-th row stored in row of an n x n matrix A    *
 *  with ld lower co-diagonals and ud upper co-diagonals into the row *
 *  prow of length ld+ud+1 as needed for the routines band and bando. *
 *  The original row called row thus becomes the no-th row of the     *
 *  input matrix  packmat for the functions band and bando.           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 2 )                                    *
 *               row dimension of the matrix                          *
 *      ld       int ld; ( ld >= 0 )                                  *
 *               number of lower co-diagonals                         *
 *      ud       int ud; ( ud >= 0 )                                  *
 *               number of upper co-diagonals                         *
 *      no       int no;                                              *
 *               row index of the matrix to be converted              *
 *      row      REAL   row[];                                        *
 *               no-th row of length n of the original matrix         *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      prow     REAL   prow[];                                       *
 *               no-th row in condensed form of length ld+ud+1.       *
 *                                                                    *
 *====================================================================*/
{
  register k, m, j;

  m = 0;
  k = ld - no;
  while ( k-- > 0 ) prow[m++] = ZERO;
  for (j = 0; j < n; j++)
    if ( (no - j  <= ld) && (j - no <= ud) ) prow[m++] = row[j];

  k = ld + ud + 1 - m;
  while ( k-- > 0 ) prow[m++] = ZERO;
  return (0);
}


int unpack              /* uncondense row ............................*/
/*.IX{unpack}*/
         (
          int     n,              /* size of matrix ..................*/
          int     ld,             /* number of lower co-diagonals ....*/
          int     ud,             /* number of upper co-diagonals ....*/
          int     no,             /* row index .......................*/
          REAL    row[],          /* uncondensed row .................*/
          REAL    prow[]          /* condensed row ...................*/
         )
/*====================================================================*
 *                                                                    *
 *  unpack converts the no-th row stored in prow of the (n x ld+ud+1) *
 *  band matrix in condensed form  with ld lower co-diagonals and ud  *
 *  upper co-diagonals into the no-th row of length n of the original *
 *  matrix A. This is the inverse operation to pack.                  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 2 )                                    *
 *               row dimension of the matrix                          *
 *      ld       int ld; ( ld >= 0 )                                  *
 *               number of lower co-diagonals                         *
 *      ud       int ud; ( ud >= 0 )                                  *
 *               number of upper co-diagonals                         *
 *      no       int no;                                              *
 *               row index of the row to be converted                 *
 *      prow     REAL   row[];                                        *
 *               no-th row of length ld+ud+1 of the condensed matrix  *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      row      REAL   row[];                                        *
 *               no-th row in uncondensed form of length n.           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Macros in use :  min, max                                        *
 *   ===============                                                  *
 *                                                                    *
 *====================================================================*/
{
  register  m, k;
  int       i, j;

  m = 0;
  k = i = no - ld;
  while ( k-- > 0 ) row[m++] = ZERO;

  k = min (n - m + ld, ud + ld);
  for (j = max (0, -i); j <= k; j++)
    row[m++] = prow[j];

  k = n - m;
  while ( k-- > 0 ) row[m++] = ZERO;
  return (0);
}

/* -------------------------- END fpack.c --------------------------- */
