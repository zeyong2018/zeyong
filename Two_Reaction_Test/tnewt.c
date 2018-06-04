#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*--------------------------------------------------------------------*/
/* Test program for newton method                                     */
/*--------------------------------------------------------------------*/


#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>
#include <tnfct1.h>


typedef struct FUNCTION_ENTRY
        {
          int      n;           /* Dimension of the system ...........*/
          FNFCT    fct;         /* Non linear function ...............*/
          JACOFCT  jaco;        /* Jacobi matrix .....................*/
          REAL     x0[4];       /* starting vector ...................*/
        } FUNCTION_ENTRY;

static FUNCTION_ENTRY ftable[] =
       {
         { 3, fx1, dfx1, { (REAL)1.5,  THREE,   ONE,  ZERO } },
         { 2, fx2, dfx2, { TWO,        TWO,     ZERO, ZERO } },
         { 4, fx3, dfx3, { TWO,        TWO,     TWO,  TWO }  },
         { 2, fx4, dfx4, { ONE,  -(REAL)2.5,    ZERO, ZERO } },
         { 2, fx5, dfx5, { ONE,        ONE,     ZERO, ZERO } },
         { 2, fx6, dfx6, { TWO,       -ONE,     ZERO, ZERO } },
         { -1,fx7, dfx7, { ONE,       -ONE,     THREE, (REAL)4.0 } },
         { 2, fx8, dfx8, { TWO,        TWO,     ZERO, ZERO } },
         { 2, fx9, dfx9, { TWO,        TWO,     ZERO, ZERO } },
         { 0, NULL, NULL, { ZERO,      ZERO,    ZERO, ZERO } }
       };

int main (void)
{
  REAL *fval, *x;
  int i, k, j, rc, kmax, prim, iter, n;
  char * proto = NULL;
  JACOFCT df;
  void *vmblock;

  WriteHead ("Newton's method (multi-dimensional)");

  prim = 0;
  kmax = 10;
  n = 20;       /* maximum dimension .................................*/

  vmblock = vminit();
  fval = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  x    = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
  {
    LogError ("No Memory", 0, __FILE__, __LINE__);
    return 1;
  }

  for (j = 0; j < 2; j++)
  {
    if (j == 0)
      WriteHead ("Use Jacobi Matrix");
    else
      WriteHead ("Use Approximation of the Jacobi Matrix");

    for (k = 0; ; k++)    /* for all functions in table ..............*/
    {
      n = ftable[k].n;
      if (n == 0) break;

      if (n < 0) n = 2;     /* set up arbritarily ....................*/

      if (j == 0)
        df = ftable[k].jaco;
      else
        df = NULL;

      for (i = 0; i < n; i++)
        x[i] = ftable[k].x0[i];

      rc = newt (n, x, ftable[k].fct, df, kmax, prim, proto,
                 fval, &iter, (REAL)(TWO * MACH_EPS));

      printf ("\nFunction number = %d\n", k + 1);
      printf ("Return Code     = %d\n", rc);
      printf ("Iterations      = %d\n\n", iter);
      printf ("Result:\n");

      if (rc <= 0)
      {
        for (i = 0; i < n; i++)
        {
          printf ("x[%2d] = ", i);
          printf (FORMAT_LE, x[i]);
          printf ("    f[%2d] = ", i);
          printf (FORMAT_LE, fval[i]);
          printf ("\n");
        }
      }
      else
      {
        LogError ("newt", rc,  __FILE__, __LINE__);
      }
    }   /* end of k */

    printf ("\n");
  }   /* end of j  */

  WriteEnd ();

  return (0);
}
