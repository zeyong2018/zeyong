#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/
/*.KA{C 6}{Systems of Nonlinear Equations}
          {Systems of Nonlinear Equations}*/
/*.FE{C 6.2.1.2}{Damped Newton Method for Systems}
                {Damped Newton Method for Systems}*/

/*.BE*/
/* -------------------------- MODULE fnewt.c ------------------------ */

#include <basis.h>
#include <vmblock.h>
#include <u_proto.h>


#define ITERMAX 300                    /* Maximal number of iterations*/

static FILE *fp;                    /* File pointer for protocol file */


static int japprox        /* approximate the Jacobi matrix ...........*/
                   (int     n,
                    REAL    x[],
                    REAL *  jmat[],
                    FNFCT   fct,
                    REAL    f0[],
                    REAL    tmpvec[]
                   );

static REAL l2norm          /* L2 Vector norm ........................*/
                    (int   n,
                      REAL  x[]
                    );

static int protopen        /* open protocol file .....................*/
                    (int     n,
                     REAL    x[],
                     int     kmax,
                     int     prim,
                     char *  pfile
                    );

static int protwrite       /* write into protocol file ...............*/
                     (int   iter,
                      REAL  fxnorm,
                      REAL  dnorm,
                      int   k
                     );

static void protclose       /* close protocol file ...................*/
                      (int   n,
                       int   iter,
                       REAL  x[],
                       REAL  fvalue[],
                       int   rc
                      );


/* Macro to shorten 'return' in newt .................................*/

#define RETURN(code)                                   \
    if (flag) protclose (n, *iter, x, fvalue, code);   \
    vmfree(vmblock);                                   \
    return (code)

/*.BA*/

int newt                 /* Multidimensional Newton method ...........*/
/*.IX{newt}*/
         (
          int       n,            /* size of system ..................*/
          REAL      x[],          /* Starting/solution vector ........*/
          FNFCT     fct,          /* Function ........................*/
          JACOFCT   jaco,         /* Function for Jacobi matrix ......*/
          int       kmax,         /* Maximal number of damped steps ..*/
          int       prim,         /* Maximal number of basic steps ...*/
          char *    pfile,        /* Name of the protocol file .......*/
          REAL      fvalue[],     /* Function value at solution ......*/
          int *     iter,         /* number of iteration steps .......*/
          REAL      eps           /* error bound .....................*/
         )
/*====================================================================*
 *                                                                    *
 *    newt determines a solution of the non-linear system of equations*
 *                                                                    *
 *                  f0 (x[0],...,x[n-1])     = 0                      *
 *                  f1 (x[0],...,x[n-1])     = 0                      *
 *                  :                                                 *
 *                  f(n-1) (x[0],...,x[n-1]) = 0                      *
 *                                                                    *
 *    using the damped Newton method.                                 *
.BE*)
 *    Here the function fct must be externally given, the Jacobi      *
 *    matrix of the partial derivatives of fct can either be given    *
 *    as another function. If this is not possible or not desired,    *
 *    the Jacobi matrix must be preassigned to be the NULL pointer;   *
 *    and the Jacobi matrix will be internally approximated by the    *
 *    formard difference matrix.                                      *
 *                                                                    *
 *    The method converges quadratically for suitably chosen starting *
 *     vectors if the problem is solvable.                            *
 *                                                                    *
 *    If the nonlinear system is solvable, then the Newton iteration  *
 *    converges to a solution depending on :                          *
 *                                                                    *
 *      1.  Starting vector of the iteration                          *
 *      2.  number of basic iteration steps                           *
 *      3.  Maximal allowed number of damped steps.                   *
 *                                                                    *
 *    We use three break-off criteria :                               *
 *                                                                    *
 *      1.  The L2 norm of the differenz deltax between the current   *
 *          x-value and the previous one is less than eps * norm(x),  *
 *      2.  the L2 norm of the functional value at the new x is less  *
 *          eps, or                                                   *
 *      3.  the maximal number of iterations has been performed.      *
 *                                                                    *
 *    Here  eps denotes the desired accuracy,  eps >= MACH_EPS.       *
 *                                                                    *
 *    Intermediate results are collected in a protocol file if the    *
 *    input parameter pfile is different from " ", a blank.           *
 *    If the protocol file name with the letter a or A, the protocol  *
 *    file is extensive, otherwise a shortened file is kept only.     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Applications:                                                     *
 *  =============                                                     *
 *      To solve nonlinear systems of equations with n functions and  *
 *      n variables.                                                  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Literature:                                                       *
 *  ===========                                                       *
 *      Conte, S.D., de Boor, C.: Elementary Numerical Analysis, an   *
 *      algorithmic approch. New York - Sidney - Toronto,             *
 *      3. ed. 1980.  [CONT80]                                        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  =================                                                 *
 *      n        int n; (n > 1).                                      *
 *               number of equations and number of unknown            *
 *      x        REAL   x[];                                          *
 *               Starting vector for the iteration.                   *
 *      fct      int fct();                                           *
 *               Function, that evaluates the n function values       *
 *               f0, ..., f(n-1)                                      *
 *               The defining program for fct has the form :          *
 *                                                                    *
 *                  int fct (int n, REAL x[], REAL fvalue[])          *
 *                  {                                                 *
 *                    fval[0] = ......       ;                        *
 *                    :                                               *
 *                    fval[n-1] = ....       ;                        *
 *                    return (0);                                     *
 *                  }                                                 *
 *                                                                    *
 *               fct return an error code which is 0 if the function  *
 *               was evaluated successfully.                          *
 *               fval denotes the functional value of of fct at x for *
 *               x = (x[i]), i = 0, ..., n-1.                         *
 *                                                                    *
 *      jaco     int jaco();                                          *
 *               user supplied function for the Jacobi matrix of fct  *
 *               jaco has the form :                                  *
 *                                                                    *
 *                  int jaco (int n, REAL x[], REAL *mem[])           *
 *                  {                                                 *
 *                   REAL **df;                                       *
 *                   df = mem;                                        *
 *                   for (i = 0; i < n; i++)                          *
 *                      for (j = 0; j < n; j++)                       *
 *                         df[i][j] = ...;                            *
 *                   return (0);                                      *
 *                  }                                                 *
 *                                                                    *
 *               Here for each i, j  .....  denotes the explicit      *
 *               partial derivative of the ith function fi with       *
 *               respect to x(j).                                     *
 *               mem is the memory location which contains the Jacobi *
 *               matrix upon execution. The return value must equal 0.*
 *               Alternatively jaco can be set equal to the zero      *
 *               matrix NULL. In this case the Jacobi matrix is       *
 *               approximated by the forward difference quotient      *
 *               matrix.                                              *
 *      kmax     int kmax;  ( 0 <= kmax <= 10 )                       *
 *               Maximal number of damped steps.                      *
 *               kmax = 0  ==> ordinary Newton method;                *
 *               kmax = 4  is a good initial value to start testing   *
 *      prim     int prim;                                            *
 *               number of basic Newtin steps, i.e., after prim       *
 *               iterations the Jacobi matrix is reevaluated.         *
 *               Experiment with prim =  0, 1, 2, 3.                  *
 *      pfile    char *pfile;                                         *
 *               Name of the protocol file.                           *
 *               = NULL or                                            *
 *               = " "    : no intermediate output.                   *
 *               = "A ..." or                                         *
 *               = "a ...": detailed report of intermediate data      *
 *               = other  : Starting value and end result             *
 *      eps      desired accuracy;                                    *
 *               If eps <= 0.0, we set eps = 4 * MACH_EPS.            *
 *                                                                    *
 *  Output parameters:                                                *
 *  ==================                                                *
 *      x        solution                                             *
 *      fvalue   REAL   fvalue[];                                     *
 *               Functional value at the solution                     *
 *      iter     int *iter;                                           *
 *               number of iterations performed.                      *
 *                                                                    *
 *  Return value :                                                    *
 *  =============                                                     *
 *      = -1     Warning: solution with L2 norm of fvalue > 128 * eps *
 *      =  0     solution found with L2 norm of fvalue <= eps         *
 *      =  1     wrong input : n < 2 or kmax < 0 or prim < 0          *
 *      =  2     lack of memory                                       *
 *      =  3     Jacobi matrix is singular                            *
 *      =  4     Iteration maximum exceeded                           *
 *      =  5     Protocol file cannot be opened                       *
 *      =  6     Protocol file does not accept input                  *
 *      =  7     user supplied function fct cannot be evaluated       *
 *      =  8     error in computeing Jacobi matrix                    *
 *      =  9     Norm of the function value > MAXROOT -> Divergence   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   functions in use  :                                              *
 *   ===================                                              *
 *                                                                    *
 *             int  gauss ():    solves the linear system             *
 *                                 jaco * deltax = f.                 *
 *   static   REAL  l2norm():    computes the L2 norm of a vector     *
 *   static    int  protopen (): Open protocol file                   *
 *   static    int  protwrite(): Write onto protocol file             *
 *   static   void  protclose(): close protocol file                  *
 *            void *vmalloc():   allocate vector or matrix            *
 *            void vmfree():     free list of vectors and matrices    *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used   :   NULL, ITERMAX, MACH_EPS                     *
 *   ==================                                               *
 *                                                                    *
 *   Macros: RETURN                                                   *
 *   =======                                                          *
 *                                                                    *
.BA*)
 *====================================================================*/
/*.BE*/
{
  int  i, k, rc, count, vordet, cas;
  int  flag, flag1, approx, *perm = NULL;
  void *vmblock = NULL;

  REAL **jmat = NULL,               /* For  LU decomposition in gauss */
       *deltax = NULL, *xtemp, *fvalue0, *tmpvec,
       fxnorm, fxnorm1, dnorm, omega;


  if (n < 2 || kmax < 0 || prim < 0)                   /* wrong input */
    return (1);

  if (x == NULL || fct == NULL) return (1);

  if (eps < MACH_EPS)
    eps = (REAL) ((REAL)4.0 * MACH_EPS);
                                          /* open protocol file if    */
  flag = flag1 = 0;                       /* desired                  */
  if (pfile)
  {
    flag = (pfile[0] != ' ');
    flag1 = (pfile[0] == 'a') || (pfile[0] == 'A');
  }

  rc = 0;                                 /* initialize return code   */

  if (flag) rc = protopen (n, x, kmax, prim, pfile);

  if (rc)                     /* Protocol file cannot be opened       */
  {
    flag = 0;
    RETURN (5);
  }

  if (jaco == NULL)           /* Approximate Jacobi matrix by         */
    approx = 1;               /* forward difference matrix            */
  else
    approx = 0;

  *iter = 0;                  /* initialize iteration counter         */
  count = prim;               /* prim steps for fixed Jacobi matrix   */
                              /* im 1. iteration compute the J matrix */

  rc = (*fct) (n, x, fvalue); /* Functional value at starting vector  */
  if (rc)
  {
    RETURN (7);
  }
                                       /* allocate storage to vectors */
                                       /* deltax, xtemp, fvalue0,     */
  vmblock = vminit();                  /* tmpvec                      */
  deltax  = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  xtemp   = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  fvalue0 = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  tmpvec  = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  perm    = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));
                                   /* Speicher fuer Jacobi Matrix     */
  jmat    = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  if (! vmcomplete(vmblock))
  {
    if (flag)
      protclose (n, *iter, x, fvalue, 2);   \
    vmfree(vmblock);
    return 2;
  }

  fxnorm = l2norm (n, fvalue);   /* L2 norm of function at x          */

  if (fxnorm <= eps)             /* Starting vector is solution?      */
  {
    RETURN (0);
  }

  do   /* Newton iteration */
  {
    (*iter)++;                   /* up iterations counter             */

    if (count < prim)            /* if counter < prim, use basic      */
    {                            /* Newton method                     */
      count++;
      cas = 2;
    }
    else                        /* otherwise: recompute Jacobi matrix */
    {
      count = 0;
      cas = 0;

      if (approx)                          /* use approximate J matrix*/
        rc = japprox (n, x, jmat, fct, fvalue, tmpvec);
      else
        rc = (*jaco) (n, x, jmat);         /* or : use user function  */

      if (rc)                       /* Jacobi matrix does not compute */
      {
        RETURN (8);
      }
    }

    /* solve linear system  jmat * deltax = fvalue  for delmax;       */
    /* jmat noe contains the LU factorization which for cas = 2 need  */
    /* not be recomputed.                                             */

    rc = gauss (cas, n, jmat, jmat, perm, fvalue, deltax, &vordet);

    if (rc)                          /* Jacobi matrix singular, or    */
    {                                /*  gauss  fails                 */
      RETURN (3);
    }

    omega = TWO;                   /* initialize damping factor       */
    k = -1;

    do    /*  Daempfung  */
    {
      k++;
      omega *= 0.5;                             /* omega = 2^{-k}     */
      for (i = 0; i < n; i++)
        xtemp[i] = x[i] - omega * deltax[i];

      rc = (*fct) (n, xtemp, fvalue);         /* fct does not compute */
      if (rc)
      {
        RETURN (7);
      }

      fxnorm1 = l2norm (n, fvalue);
      if (kmax == 0) break;             /* no damping wanted          */
      if (k == 0)
        for (i = 0; i < n; i++)         /* store function values in   */
          fvalue0[i] = fvalue[i];       /* fvalue0                    */
    }
                                              /* use damping until    */
    while (fxnorm < fxnorm1 && k <= kmax);    /* k = kmax or an x with*/
                                              /* smaller fct value has*/
                                              /* been found.          */

    if ( (0 < k  && k <= kmax) || kmax == 0 ) /* in case of damping   */
    {                                         /* or kmax = 0          */

      for (i = 0; i < n; i++)                 /* use current values   */
        x[i] = xtemp[i];
      fxnorm = fxnorm1;
      dnorm = omega * l2norm (n, deltax);
    }
    else                                      /* otherwise set        */
    {                                         /*  x = x - deltax      */
      for (i = 0; i < n; i++)
      {
        x[i] -= deltax[i];
        fvalue[i] = fvalue0[i];
      }
      fxnorm = l2norm (n, fvalue);
      dnorm = l2norm (n, deltax);
    }

    if (flag1)                            /* keep protocol if desired */
    {
      rc = protwrite (*iter, fxnorm, dnorm, k);
      if (rc)
      {
        RETURN (6);
      }
    }

  }                                    /* as long as:                 */
  while (dnorm > eps * l2norm (n, x)   /* Norm(deltax) > eps*Norm(x), */
         && fxnorm > eps               /* Norm(fx) > eps,             */
         && fxnorm < MAXROOT           /* Norm(fx) < MAXROOT,         */
         && *iter < ITERMAX);          /* iter < Iterationsmax:       */

  if (*iter >= ITERMAX) rc = 4;       /* Iteration max exceeded       */
  else
    if (fxnorm >= MAXROOT)            /* Functional value too large   */
      rc = 9;                         /* suspect divergence           */
    else
      if (fxnorm > (REAL)128.0 * eps) /* Warning: bad approximation   */
        rc = -1;

  RETURN (rc);
}


static REAL l2norm          /* L2 Vector norm .. .....................*/
/*.IX{l2norm}*/
                   (
                    int  n,
                    REAL x[]
                   )
/*====================================================================*
 *                                                                    *
 *  l2norm computes the euclidean or L2 norm of a vector              *
 *  x = (x[0],x[1],...,x[n-1]). It avoides underflow.                 *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constant used : EPSQUAD                                          *
 *   ===============                                                  *
 *                                                                    *
 *   Macros: ABS, SQRT                                                *
 *   ======                                                           *
 *====================================================================*/
{
  int  i, j;
  REAL scale, sum, tmp, xiabs;

  if (n <= 0) return (ZERO);                   /* n <= 0 ==> Norm = 0 */

  for (i = 0; i < n; i++)
    if (x[i] != ZERO) break;

  if (i == n) return (ZERO);                   /* zero vector         */

  scale = ABS (x[i]);
  if (i == n - 1) return (scale);        /* only one component  != 0  */

  j = i + 1;
  for (sum = ONE, i = j; i < n; i++)
  {
    xiabs = ABS (x[i]);
    if (xiabs <= scale)                  /* scale = previous max      */
    {                                    /* of  ABS(x[i])             */
      tmp = xiabs / scale;
      if (tmp > EPSQUAD)
        sum += tmp * tmp;                /* sum = sum + temp*temp     */
    }
    else
    {
      tmp = scale / xiabs;
      if (tmp <= EPSQUAD) tmp = ZERO;
      sum *= tmp * tmp;
      sum += ONE;                     /* sum = sum * temp * temp + 1  */
      scale = xiabs;
    }
  }

  return (scale * SQRT (sum));
}


static int japprox        /* approximate Jacobi matrix ...............*/
/*.IX{japprox}*/
                   (
                    int     n,
                    REAL    x[],
                    REAL *  jmat[],
                    FNFCT   fct,
                    REAL    f0[],
                    REAL    tmpvec[]
                   )
/*====================================================================*
 *                                                                    *
 *  japprox approximates the Jacobi matrix of a vector valued function*
 *  by forming forward difference quotients. This function is an      *
 *  alternative to the explicit  Jacobi matrix in the Newton algorithm*
 *  newt.                                                             *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  =================                                                 *
 *      n        int n;                                               *
 *               number of equations and number of unknown            *
 *      x        REAL   x[];                                          *
 *               Vector, where the Jacobi matrix is to be determined  *
 *      jmat     REAL   * jmat[];                                     *
 *               pointer for output storage                           *
 *      fct      REAL   *fct();                                       *
 *               Function, which computes the n values in  f0, ...,   *
 *               f(n-1) (same as in newt)                             *
 *      f0       REAL   f0[];                                         *
 *               Vector of function value at x                        *
 *      tmpvec   REAL   tmpvec[];                                     *
 *               aux vector                                           *
 *                                                                    *
 *  Output parameter:                                                 *
 *  =================                                                 *
 *      jmat     REAL   *jmat[];                                      *
 *               forward difference approximation for the Jacobi      *
 *               matrix.                                              *
 *                                                                    *
 *  Return value :                                                    *
 *  =============                                                     *
 *              = 0: Jacobi matrix approximation found                *
 *              other : error                                         *
 *====================================================================*
 *                                                                    *
 *   Constants in use  :  NULL, EPSROOT                               *
 *   ===================                                              *
 *                                                                    *
 *====================================================================*/
{
  REAL xj, *f1, h, denom;
  int  i, j, rc = 0;

  if (jmat == NULL) return (-1);

  f1 = tmpvec;
  if (f1 == NULL) return (-1);

  for (j = 0; j < n; j++)
  {
    xj = x[j];

    h = (REAL)(EPSROOT * HALF);  /* search for least h >= EPSROOT     */
    do                           /* with xj != x[j]                   */
    {
      h += h;
    } while (xj + h == x[j]);

    x[j] += h;                   /* temporary: x[j] = x[j] + h        */
    denom = ONE / h;

    rc = (*fct) (n, x, f1);      /* f1 is the function value at the   */
                                 /* new x-value                       */

    x[j] = xj;                   /* reassign old x component          */

    if (rc) return (rc);

    for (i = 0; i < n; i++)      /* forward difference quotient       */
      jmat[i][j] = ( f1[i] - f0[i] ) * denom;
  }

  return (0);
}


static int protopen        /* open protokol file .....................*/
/*.IX{protopen}*/
                    (int     n,
                     REAL    x[],
                     int     kmax,
                     int     prim,
                     char *  pfile
                    )
/*====================================================================*
 *                                                                    *
 * popen opens the protokol file pfile im append-mode and writes the  *
 * input parameters for Newton onto this file.                        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input             : values of n, x, kmax, prim, pfile from newt.  *
 *  Return value : =0 : Protocol file has been opened                 *
 *                 =1 : file not opened                               *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used    :  NULL                                        *
 *   ===================                                              *
 *                                                                    *
 *====================================================================*/
{
  int i;

  fp = fopen (pfile, "a");
  if (fp == NULL) return (1);

  fprintf (fp, "Damped Newton method \n");
  fprintf (fp, "--------------------- \n\n");
  fprintf (fp, "size of system                  n   : %3d \n",n);
  fprintf (fp, "number of basic steps           prim: %3d \n",prim);
  fprintf (fp, "Maximal number of damped steps  kmax: %3d\n\n",kmax);

  fprintf (fp, "Starting vector x :\n");
  for (i = 0; i < n; i++)
  {
    fprintf (fp, "\t x[%2d] = ", i);
    fprintf (fp, FORMAT_LE, x[i]);
    fprintf (fp, "\n");
  }
  if ( pfile[0] == 'a' || pfile[0] == 'A' )
    fprintf (fp, "\n Iter\t Norm(f)\t Norm(deltax)\tk \n\n");
  else fprintf (fp, "\n");

  fflush (fp);
  return (0);
}


static int protwrite       /* write onto protocol file ...............*/
/*.IX{protwrite}*/
                     (int     iter,
                      REAL    fxnorm,
                      REAL    dnorm,
                      int     k
                     )
/*====================================================================*
 *                                                                    *
 *  protwrite wries onto pfile :                                      *
 *     -  the actual iteration number,                                *
 *     -  L2 norm of the current functional value,                    *
 *     -  L2 norm of the actual step size,                            *
 *     -  Damping factor.                                             *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input : values of  iter, fxnorm, dnorm, k from newt.              *
 *                                                                    *
 *====================================================================*
 *  Return value : =  0: ok                                           *
 *                 != 0: no input into file possible                  *
 *                                                                    *
 *====================================================================*/
{
  int rc;

  rc = fprintf (fp, " %3d\t", iter);
  rc = fprintf (fp, FORMAT_LE, fxnorm);
  rc = fprintf (fp, "\t");
  rc = fprintf (fp, FORMAT_LE, dnorm);
  rc = fprintf (fp, "\t%2d \n", k);
  if (rc <= 0) return (rc);
  fflush (fp);
  return (0);
}


static void protclose       /* close protocol file ...................*/
/*.IX{protclose}*/
                      (int   n,
                       int   iter,
                       REAL  x[],
                       REAL  fvalue[],
                       int   rc
                      )
/*====================================================================*
 *                                                                    *
 *  pclose writes onto the protocol file :                            *
 *     -  size n of problem,                                          *
 *     -  total number of iterations,                                 *
 *     -  solution x of the nonlinear system,                         *
 *     -  functional value fvalue at x,                               *
 *     -  Return value of newt;                                       *
 *  and then closes the file.                                         *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input : values of  n, iter, x, fvalue, rc from newt.              *
 *                                                                    *
 *====================================================================*/
{
  int i;

  fprintf (fp,"\n\nNumber of Newton iterations  iter: %3d \n",iter);
  fprintf (fp,"Return value                      rc : %3d \n\n",rc);
  fprintf (fp,"Approximate solution x:\t\t Function value:\n\n");
  for (i = 0; i < n; i++)
  {
    fprintf (fp,"x[%2d] = ", i);
    fprintf (fp, FORMAT_LE, x[i]);
    fprintf (fp, "\t\t f(%2d) = ", i);
    fprintf (fp, FORMAT_LE, fvalue[i]);
    fprintf (fp, "\n");
  }
  fprintf (fp,
     "\n--------------------------------------------------------\n");
  fclose (fp);
}

/* --------------------------- END fnewt.c -------------------------- */
