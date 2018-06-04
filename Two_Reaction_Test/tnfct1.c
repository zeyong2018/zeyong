#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*--------------------------------------------------------------------*/
/* Test functions for fnewt                                           */
/*                                                                    */
/*       Newton's method for several dimensions                       */
/*                                                                    */
/*                                                                    */
/*--------------------------------------------------------------------*/

#include <basis.h>
#include <u_proto.h>
#include <tnfct1.h>


int fx1 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   1                         ***/
/**********************************************************************/
{
  if (n != 3) return -1;
  if (x[2] == (REAL) 0.0) return -1;

  f[0] = x[0] - x[1] + SIN (x[0]) + COS (x[1]) - PI / (REAL) 8.0;
  f[1] = SQR (x[0]) - (REAL)4.0 * SQR (x[1]) - LOG (SQR (x[2]));
  f[2] = x[2] - (REAL)1.0 + (REAL)2.0 * x[0] - (REAL)4.0 * x[1];

  return (0);
}

int dfx1 (int n, REAL x[], REAL *df[])
{
  if (n != 3) return -1;
  if (x[2] == (REAL)0.0) return -1;

  df[0][0] =  (REAL)1.0 + COS (x[0]);
  df[0][1] = -(REAL)1.0 - SIN (x[1]);
  df[0][2] =  (REAL)0.0;
  df[1][0] =  (REAL)2.0 * x[0];
  df[1][1] = -(REAL)8.0 * x[1];
  df[1][2] = -(REAL)2.0 / x[2];
  df[2][0] =  (REAL)2.0;
  df[2][1] = -(REAL)4.0;
  df[2][2] =  (REAL)1.0;

  return (0);
}

int fx2 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   2                         ***/
/**********************************************************************/
{
  if (n != 2) return -1;
  if (x[0] <= (REAL)0.0) return -1;

  f[0] = SQR (x[0]) + x[0] + SQR (x[1]) - (REAL)2.0;
  f[1] = SQR (x[0]) + x[1] - SQR (x[1]) - (REAL)1.0 + LOG (x[0]);

  return (0);
}


int dfx2 (int n, REAL x[], REAL *df[])
{
  if (n != 2) return -1;
  if (x[0] == (REAL)0.0) return -1;

  df[0][0] = x[0] + x[0] + (REAL)1.0;
  df[0][1] = x[1] + x[1];
  df[1][0] = x[0] + x[0] + (REAL)1.0 / x[0];
  df[1][1] = (REAL)1.0 - x[1] + x[1];

  return (0);
}

int fx3 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   3                         ***/
/**********************************************************************/
{
  if (n != 4) return -1;
  if (x[0] < MACH_EPS) return -1;
  if ((x[1] + x[2]) < MACH_EPS) return -1;

  f[0] = (REAL)10.0 * x[0] + x[1] + x[2] + x[3] - (REAL)20.0 +
         SQR (SIN (x[0])) + SQR (COS (x[1]));

  f[1] = x[0] + (REAL)20.0 * x[1] + x[2] + x[3]
         - (REAL)48.0 + (REAL)1.0 / POW (x[0],6.0);

  f[2] = SQR (x[0] + x[1]) + (REAL)30.0 * x[2] + x[3]
         - (REAL)97.0 + LOG (x[0]) + LOG (x[1] + x[2]);

  f[3] = x[0] + x[1] + x[2] + (REAL)40.0 * x[3]
         - (REAL)166.0 + SQR (x[0]);

  return (0);
}

int dfx3 (int n, REAL x[], REAL *df[])
{
  if (n != 4) return -1;
  if (x[0] < MACH_EPS) return -1;
  if (ABS (x[1] + x[2]) < MACH_EPS) return -1;

  df[0][0] = (REAL)10.0 + SIN (x[0] + x[0]);
  df[0][1] = (REAL)1.0 - SIN (x[1] + x[1]);
  df[0][2] = (REAL)1.0;
  df[0][3] = (REAL)1.0;
  df[1][0] = (REAL)1.0 - (REAL)6.0 / POW(x[0],7.0);
  df[1][1] = (REAL)20.0;
  df[1][2] = (REAL)1.0;
  df[1][3] = (REAL)1.0;
  df[2][0] = (REAL)2.0 * (x[0] + (REAL)2.0 * x[1]) + (REAL)1.0 / x[0];
  df[2][1] = (REAL)4.0 * (x[0] + (REAL)2.0 * x[1]) +
             (REAL)1.0 / (x[1] + x[2]);
  df[2][2] = (REAL)30.0 + (REAL)1.0 / (x[1] + x[2]);
  df[2][3] = (REAL)1.0;
  df[3][0] = (REAL)1.0 + x[0] + x[0];
  df[3][1] = (REAL)1.0;
  df[3][2] = (REAL)1.0;
  df[3][3] = 40.0;

  return (0);
}

int fx4 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   4                         ***/
/**********************************************************************/
{
  if (n != 2) return -1;

  f[0] = x[0] - EXP (SIN (x[1]));
  f[1] = PI * EXP (SIN (PI * x[0])) + x[1];
  return (0);
}

int dfx4 (int n, REAL x[], REAL *df[])
{
  if (n != 2) return -1;

  df[0][0] = (REAL)1.0;
  df[0][1] = - EXP (SIN (x[1])) * COS (x[1]);
  df[1][0] = SQR(PI) * EXP (SIN (PI*x[0])) * COS (PI*x[0]);
  df[1][1] = (REAL)1.0;

  return (0);
}


int fx5 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   5                         ***/
/**********************************************************************/
{
  if (n != 2) return -1;

  f[0] = x[0] * x[0] * x[0] + (REAL)1.0;
  f[1] = x[1] * x[1] * x[1] + (REAL)1.0;

  return (0);
}

int dfx5 (int n, REAL x[], REAL *df[])
{
  if (n != 2) return -1;

  df[0][0] = (REAL)3.0 * SQR (x[0]);
  df[0][1] = (REAL)0.0;
  df[1][1] = (REAL)3.0 * SQR (x[1]);
  df[1][0] = (REAL)0.0;
  return (0);
}


int fx6 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   6                         ***/
/**********************************************************************/
{
  if (n != 2) return -1;

  f[0] = SIN (x[0]) - x[1];
  f[1] = x[0] - COS (x[1]);

  return (0);
}

int dfx6 (int n, REAL x[], REAL *df[])
{
  if (n != 2) return -1;
  df[0][0] = COS (x[0]);
  df[0][1] = -(REAL)1.0;
  df[1][0] =  (REAL)1.0;
  df[1][1] = SIN (x[1]);

  return(0);
}


int fx7 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   7                         ***/
/**********************************************************************/
{
  register int i;
  REAL fi;

  for (fi = (REAL)0.0, i = 0; i < n; i++)
    fi += x[i];

  for (i = 0; i < n; i++)
    f[i] = fi - x[i] - i - (REAL)1.0 + SQR (x[i]);

  return (0);
}

int dfx7 (int n, REAL x[], REAL *df[])
{
  register int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      df[i][j] = (REAL)1.0;

  for (i = 0; i < n; i++)
    df[i][i] = x[i] + x[i];

  return (0);
}


int fx8 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   8                         ***/
/**********************************************************************/
{
  if (n != 2) return -1;

  f[0] = SQR (x[0]) + SQR (x[1]) - (REAL)5.0;
  f[1] = x[0] + x[0] - x[1];

  return (0);
}


int dfx8 (int n, REAL x[], REAL *df[])
{
  if (n != 2) return -1;

  df[0][0] = x[0] + x[0];
  df[0][1] = x[1] + x[1];
  df[1][0] = (REAL)2.0;
  df[1][1] = -(REAL)1.0;

  return (0);
}


int fx9 (int n, REAL x[], REAL f[])
/**********************************************************************/
/***                    f u n c t i o n   9                         ***/
/**********************************************************************/
{
  if (n != 2) return -1;

  f[0] = (REAL)0.5 * SQR (x[0] - (REAL)1.0) +
         (REAL)0.25 * SQR (x[1]) - (REAL)1.0;
  f[1] = x[1] - SIN (x[0]);

  return (0);
}


int dfx9 (int n, REAL x[], REAL *df[])
{
  if (n != 2) return -1;

  df[0][0] = x[0] - (REAL)1.0;
  df[0][1] = (REAL)0.5 * x[1];
  df[1][0] = - COS (x[0]);
  df[1][1] = (REAL)1.0;

  return (0);
}
