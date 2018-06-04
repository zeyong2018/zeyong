#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE fzyfsy.c ------------------------ */

#include <basis.h>                    /*  for  MACH_EPS, REAL,        */
#include <glsp.h>                     /*  for  fzyfsy, fzyfsz, fzyfsl */

#define MACH4_EPS (FOUR*MACH_EPS)
/*.BA*/

int fzyfsy (int   n,
/*.IX{fzyfsy}*/
            REAL* md,     REAL* ud1,    REAL* ud2,
            REAL* rs,
            REAL* x,
            REAL* cmd,
            REAL* cld_1,  REAL* cld_2,
            REAL* cld_l2, REAL* cld_l1,
            REAL* bud_1,  REAL* bud_2,
            REAL* brs_2,  REAL* brs_1
           )
/***********************************************************************
*                                                                      *
*  fzyfsy  solves the linear system  A * X = RS  for a symmetric,      *
*  almost cyclic five diagonal system matrix A.                        *
.BE*)
*                                                                      *
*  The matrix A is given by its diagonals  md, ud1 and ud2 .           *
*  The linear system has the form :                                    *
*                                                                      *
*  md[1]*x[1] + ud1[1]*x[2] + ud2[1]*x[3]                              *
*                           + ud2[n-1]*x[n-1] + ud1[n]*x[n]  = rs[1]   *
*  ud1[1]*x[1] + md[2]*x[2] + ud1[2]*x[3]                              *
*                           + ud2[2]*x[4] + ud2[n]*x[n]      = rs[2]   *
*                                                                      *
*  ud2[i-2]*x[i-2] + ud1[i-11]*x[i-1] + md[i]*x[i]                     *
*                           + ud1[i]*x[i+1] + ud2[i]*x[i+2]  = rs[i]   *
*                           for i=3, ..., n-2                          *
*                                                                      *
*  ud2[n-1]*x[1] + ud2[n-3]*x[n-3] + ud1[n-2]*x[n-2]                   *
*                           + md[n-1]*x[n-1] + ud1[n-1]*x[n] = rs[n-1] *
*  ud1[n]*x[1] + ud2[n]*x[2] + ud2[n-2]*x[n-2]                         *
*                           + ud1[n-1]*x[n-1] + md[n]*x[n]   = rs[n]   *
*                                                                      *
* ==================================================================== *
*                                                                      *
*  Input parameters:                                                   *
*  -----------------                                                   *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*                                                                      *
*   n       int/--              number of equations  (n > 5)           *
*   md      REAL  /[n+1]        }  diagonal and co-diagonal entries    *
*   ud1     REAL  /[n+1]        }  of system matrix (the 0th entry is  *
*   ud2     REAL  /[n+1]        }   not used )                         *
*   rs      REAL  /[n+1]        right hand side (0th entry not used)   *
*                                                                      *
*  Output parameters:                                                  *
*  ------------------                                                  *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*   cmd     REAL  /[n+1]        }  co-digonals of the lower triangular *
*   cld_1   REAL  /[n+1]        }  matrix C and of the upper triangular*
*   cld_2   REAL  /[n+1]        }  B with  A = C * B ;                 *
*   cld_l1  REAL  /[n-3]        }  see fzyfsz                          *
*   cld_l2  REAL  /[n-2]        }  (the 0th entries again are not used)*
*                                                                      *
*   bud_l   REAL  /[n+1]        }                                      *
*   bud_2   REAL  /[n+1]        }                                      *
*   brs_2   REAL  /[n-3]        }                                      *
*   brs_1   REAL  /[n-2]        }                                      *
*                                                                      *
*   x       REAL  /[n+1]        solution vector (position 0 not used)  *
*                                                                      *
*                                                                      *
*  Return value :                                                      *
*  --------------                                                      *
*                                                                      *
*   = 0 : no error                                                     *
*   = 1 : error in  fzyfsy                                             *
*   = 2 : n < 6                                                        *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   subprograms used :        fzyfsz, fzyfsl                           *
*   ------------------                                                 *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int error;

 /* int i; */
  if (n < 6)     return (2);
 /*
     Factor system matrix into C * B for triangular matrices C, B
 */
  error = fzyfsz (n, md, ud1, ud2, cmd, cld_1, cld_2, cld_l2, cld_l1,
          bud_1, bud_2, brs_2, brs_1);

  if (! error)                        /* factorization without error */
    fzyfsl (n, rs, x, cmd, cld_1, cld_2, cld_l2, cld_l1,
            bud_1, bud_2, brs_2, brs_1);

  return (0);
}

/* ------------------------------------------------------------------ */
/*.BA*/

int fzyfsz (int   n,
/*.IX{fzyfsz}*/
            REAL* md,     REAL* ud1, REAL* ud2,
            REAL* cmd,
            REAL* cld_1,  REAL* cld_2,
            REAL* cld_l2, REAL* cld_l1,
            REAL* bud_1,  REAL* bud_2,
            REAL* brs_2,  REAL* brs_1
           )
/***********************************************************************
*                                                                      *
*  fzyfsz  factors a symmetric, almost cyclic five diagonal matrix A   *
*  into the product of a lower triangular matrix C and an upper        *
*  triangular matrix B.                                                *
*  The matrix A is given by the vectors of its (co)diagonals md, ud1   *
*  and ud2 . (refer to  fzyfsy).                                       *
.BE*)
*                                                                      *
* ==================================================================== *
*                                                                      *
*  Input parameters:                                                   *
*  -----------------                                                   *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*   n       int/--              number of equations  (n > 5)           *
*   md      REAL  /[n+1]        }  Matrix elements in vector form      *
*   ud1     REAL  /[n+1]        }  (0th entry not used)                *
*   ud2     REAL  /[n+1]        }                                      *
*                                                                      *
*                                                                      *
*  Output parameters:                                                  *
*  ------------------                                                  *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*   cmd     REAL  /[n+1]        }  Elements of the lower and upper     *
*   cld_1   REAL  /[n+1]        }  triangular factors C and B of A :   *
*   cld_2   REAL  /[n+1]        }         A = C * B                    *
*   cld_l1  REAL  /[n-3]        }  ( the 0th entries are not used)     *
*   cld_l2  REAL  /[n-2]        }                                      *
*                                                                      *
*   bud_l   REAL  /[n+1]        }                                      *
*   bud_2   REAL  /[n+1]        }                                      *
*   brs_2   REAL  /[n-3]        }                                      *
*   brs_1   REAL  /[n-2]        }                                      *
*                                                                      *
*   DETAILS:                                                           *
*   cmd   : main diagonal of C; entries  1 to n used                   *
*   cld_1 : lower co-diagonal of C; entries  2 to n used               *
*   cld_2 : second lower co-diagonal of C; entries  3 to n used        *
*   cld_l2: second to last row of C; elements  1 to n-4 used           *
*   cld_l1: last row of C; entries  1 to n-3 used                      *
*   bud_1 : upper co-diagonal of B; elements  1 to n-1 used            *
*   bud_2 : second upper co-diagonal of B; entries  1 to n-2 used      *
*   brs_2 : (n-1)th column of B; elements  1 to n-4 used               *
*   brs_1 : nth column of B; elements  1 to n-3 used                   *
*                                                                      *
*                                                                      *
*                                                                      *
*  Return value :                                                      *
*  --------------                                                      *
*                                                                      *
*   = 0 : no error                                                     *
*   = 1 : at least one elemet on the main diagonal of C has magnitude  *
*         less than 4 * MACH_EPS, resulting in division near overflow. *
*                                                                      *
* ==================================================================== *
*                                                                      *
*   constants used :       CH_EPS, MACH4_EPS                           *
*   ----------------                                                   *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int i, j, k;
  REAL   h_var_1, h_var_2, h_var_3;
 /*
     Factor A = C * B.
     Stop if an element on the main diagonal of C decomes nearly zero.
 */
  cmd [1] = md [1];
  if (cmd[1] <= MACH4_EPS)     return (1);
  bud_1 [1] = ud1[1] / cmd[1];
  brs_2 [1] = ud2[n-1] / cmd[1];
  brs_1 [1] = ud1[n] / cmd[1];
  cld_1 [2] = ud1[1];
  cmd [2]   = md[2] - cld_1[2] * bud_1[1];
  if (cmd[2] <= MACH4_EPS)     return (1);
  brs_2 [2] = - (brs_2[1] * cld_1[2]) / cmd[2];
  brs_1 [2] = (ud2[n] - cld_1[2] * brs_1[1]) / cmd[2];

  for (i=3; i<=n-2; ++i)
  {
    j = i - 2;
    k = i - 1;
    cld_2 [i] = ud2[i-2];
    bud_2 [j] = ud2[j] / cmd[j];
    bud_1 [k] = (ud1[k] - cld_1[k] * bud_2[j]) / cmd[k];
    cld_1 [i] = ud1[i-1] - cld_2[i] * bud_1[j];
    cmd [i]   = md[i] - cld_1[i] * bud_1[k] - cld_2[i] * bud_2[j];
    if (cmd[i] <= MACH4_EPS)     return (1);
  }

  for (i=3; i<=n-4; ++i)
    brs_2 [i] = - (cld_2[i] * brs_2[i-2]
                   + cld_1[i] * brs_2[i-1]) / cmd[i];

  for (i=3; i<=n-3; ++i)
    brs_1 [i] = - (cld_2[i] * brs_1[i-2]
                   + cld_1[i] * brs_1[i-1]) / cmd[i];

  bud_2 [n-3] = (ud2[n-3] - cld_1[n-3] * brs_2[n-4] -
                                cld_2[n-3] * brs_2[n-5]) / cmd[n-3];
  bud_2 [n-2] = (ud2[n-2] - cld_1[n-2] * brs_1[n-3] -
                                cld_2[n-2] * brs_1[n-4]) / cmd[n-2];
  bud_1 [n-2] = (ud1[n-2] - cld_1[n-2] * bud_2[n-3] -
                                cld_2[n-2] * brs_2[n-4]) / cmd[n-2];

  cld_l2 [1] = ud2[n-1];
  cld_l2 [2] = -cld_l2[1] * bud_1[1];
  for (i=3; i<=n-4; ++i)
    cld_l2 [i] = -(cld_l2[i-2] * bud_2[i-2] + cld_l2[i-1] * bud_1[i-1]);

  cld_l1 [1] = ud1[n];
  cld_l1 [2] = ud2[n] - cld_l1[1] * bud_1[1];
  for (i=3; i<=n-3; ++i)
    cld_l1 [i] = -(cld_l1[i-2] * bud_2[i-2] + cld_l1[i-1] * bud_1[i-1]);

  cld_2 [n-1] = ud2[n-3] - (cld_l2[n-5] * bud_2[n-5]
                            + cld_l2[n-4] * bud_1[n-4]);
  cld_2 [n]   = ud2[n-2] - (cld_l1[n-4] * bud_2[n-4]
                            + cld_l1[n-3] * bud_1[n-3]);
  cld_1 [n-1] = ud1[n-2] - (cld_l2[n-4] * bud_2[n-4]
                            + cld_2[n-1] * bud_1[n-3]);

  h_var_1 = h_var_2 = h_var_3 = 0.0;
  for (i=1; i<=n-4; ++i)
  {
    h_var_1 += cld_l1[i] * brs_2[i];
    h_var_2 += cld_l2[i] * brs_2[i];
    h_var_3 += cld_l2[i] * brs_1[i];
  }

  cld_1 [n] = ud1[n-1] - h_var_1 - cld_l1[n-3] * bud_2[n-3]
              - cld_2[n] * bud_1[n-2];

  cmd [n-1] = md[n-1] - h_var_2 - cld_2[n-1] * bud_2[n-3]
              - cld_1[n-1] * bud_1[n-2];
  if (cmd[n-1] <= MACH4_EPS)     return (1);

  bud_1[n-1] = (ud1[n-1] - h_var_3 - cld_2[n-1] * brs_1[n-3]
                                   - cld_1[n-1] * bud_2[n-2]) /
                                     cmd[n-1];

  for (h_var_1=ZERO, i=1; i<=n-3; ++i)
    h_var_1 += cld_l1[i] * brs_1[i];
  cmd [n] = md[n] - h_var_1 - cld_2[n] * bud_2[n-2]
                                   - cld_1[n] * bud_1[n-1];
  if (cmd[n] <= MACH4_EPS)     return (1);

  return (0);
}

/* ------------------------------------------------------------------ */
/*.BA*/

int fzyfsl (int   n,
/*.IX{fzyfsl}*/
            REAL* rs,     REAL* x,
            REAL* cmd,
            REAL* cld_1,  REAL* cld_2,
            REAL* cld_l2, REAL* cld_l1,
            REAL* bud_1,  REAL* bud_2,
            REAL* brs_2,  REAL* brs_1
           )
/***********************************************************************
*                                                                      *
*  fzyfsl solves the linear system  A * X = RS  for a symmetric, almost*
*  cyclic five diagonal syatem matrix A which is given in factored form*
*  A = C * B  from   fzyfsy or fzyfsz.                                 *
.BE*)
*                                                                      *
* ==================================================================== *
*                                                                      *
*  Input parameters:                                                   *
*  -----------------                                                   *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*   n       int/--              number of equations  (n > 5)           *
*   rs      REAL  /[n+1]        right hand side (0th entry not used)   *
*                                                                      *
*   cmd     REAL  /[n+1]        }  Elements of the triangular factors  *
*   cld_1   REAL  /[n+1]        }  C and B of A                        *
*   cld_2   REAL  /[n+1]        }  ( 0th elements not used)            *
*   cld_l1  REAL  /[n-3]        }  ( for details see  fzyfsl)          *
*   cld_l2  REAL  /[n-2]        }                                      *
*                                                                      *
*   bud_l   REAL  /[n+1]        }                                      *
*   bud_2   REAL  /[n+1]        }                                      *
*   brs_2   REAL  /[n-3]        }                                      *
*   brs_1   REAL  /[n-2]        }                                      *
*                                                                      *
*                                                                      *
*  Output parameters:                                                  *
*  ------------------                                                  *
*                                                                      *
*   Name    Type/size           Meaning                                *
*  ------------------------------------------------------------------- *
*                                                                      *
*   x       REAL  /[n+1]        solution vector (0th entry not used)   *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/
{
  int i;
  REAL   h_var_1;
 /*
     Solve system by updating right hand side and backsubstitution
 */
  x[1] = rs[1] / cmd[1];                   /*  update right hand side */
  x[2] = (rs[2] - x[1] * cld_1[2]) / cmd[2];
  for (i=3; i<=n-2; ++i)
    x [i] = (rs[i] - x[i-2] * cld_2[i] - x[i-1] * cld_1[i]) / cmd[i];

  for (h_var_1=0.0, i=1; i<=n-4; ++i)
    h_var_1 += x[i] * cld_l2[i];
  x[n-1] = (rs[n-1] - h_var_1 - x[n-3] * cld_2[n-1]
            - x[n-2] * cld_1[n-1]) / cmd[n-1];

  for (h_var_1=0.0, i=1; i<=n-3; ++i)
    h_var_1 += x[i] * cld_l1[i];
  x[n] = (rs[n] - h_var_1 - x[n-2] * cld_2[n]
                          - x[n-1] * cld_1[n]) / cmd[n];

  x[n-1] -= bud_1[n-1] * x[n];                   /*  back substitute  */
  x[n-2] -= (bud_1[n-2] * x[n-1] + bud_2[n-2] * x[n]);
  x[n-3] -= (bud_1[n-3] * x[n-2] + bud_2[n-3] * x[n-1]
                                 + brs_1[n-3] * x[n]);
  for (i=n-4; i>=1; --i)
    x [i] -= (bud_1[i] * x[i+1] + bud_2[i] * x[i+2]
                                + brs_2 [i] * x[n-1] + brs_1[i] * x[n]);
  return 0;
}

/* -------------------------- END fzyfsy.c -------------------------- */
