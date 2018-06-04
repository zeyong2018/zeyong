#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE brown.c ------------------------- */

/***********************************************************************
*                                                                      *
* Brown's method for nonlinear systems                                 *
* -------------------------------------                                *
*                                                                      *
* Programing language: ANSI C                                          *
* Compiler:            Borland C++ 2.0                                 *
* Computer:            IBM PS/2 70 mit 80387                           *
* Literature:          Brown, K. M.:                                   *
*                      A quadratically convergent Newton-like method   *
*                      based upon Gaussian elimination                 *
*                      SIAM J. Numer. Anal. Vol. 6 , (1969), p. 560    *
* Source:              QuickBASIC program                              *
* Author:              Johannes Karfusehr (FORTRAN)                    *
* Adaption:            Juergen Dietel, Computer Center, RWTH Aachen    *
* Date:                8.17.1992                                       *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  REAL, ZERO, printf, LZP, FALSE, TRUE, */
                       /*       MACH_EPS, copy_vector, boolean, FABS, */
                       /*       EIGHT, NULL                           */
#include <vmblock.h>   /*  for  egen vmalloc, vmcomplete, vmfree,     */
                       /*       vminit, VEKTOR, MATRIX, IMATRIX       */
#include <brown.h>     /*  for  nlgls, brown                          */



/* ------------------------------------------------------------------ */

static void subst
/*.IX{subst}*/
                 (
                  int  n,
                  int  k,
                  int  *ihf[],
                  REAL *hf[],
                  REAL rslin[],
                  REAL x1[]
                 )

/***********************************************************************
* Solve a linear system of equations                                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n      size of system                                                *
* k      Index for coordinates: 0,...,n-1                              *
* ihf    [0..n-1,0..n] matrix, register of row and column interchanges *
* hf     [0..n-1,0..n-1] matrix, the system matrix                     *
* rslin  [0..n-1] right nahd side vector                               *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* x1     [0..n-1] solution vector                                      *
*                                                                      *
* Constants in use :  REAL, ZERO                                       *
* ==================                                                   *
*                                                                      *
***********************************************************************/

{
  REAL sum;              /* aux variable for finding x1[kmax]        */
  int  km,               /* Loop counter                             */
       kmax,             /* original row index                       */
       jsub,             /* original column index                    */
       j;                /* Loop counter                             */


  for (km = k; km > 0; km--)
  {
    kmax = ihf[km - 1][n];

    for (sum = ZERO, j = km; j < n; j++)
      jsub =  ihf[km][j],
      sum  += hf[km - 1][jsub] * x1[jsub];

    x1[kmax] = sum + rslin[km - 1];
  }
}



/* ------------------------------------------------------------------ */

/***********************************************************************
* aux vectors and matrices needed in iter4() and allocated in brown(). *
* These are defined globally for the module now so that they are       *
* available for iter4() (for actual computations) and for brown() (for *
* memory allocation).                                                  *
* Local definition only in iter4() is possible but would involve many  *
* unnecessary allocation and free calls later.                         *
***********************************************************************/

static
  int  **ihf;   /* [0..n-1,0..n] matrix eith row and column swaps     */
static
  REAL *rslin,  /* [0..n-1] right hand side vector                    */
       *dquot,  /* [0..n-1] vector for difference quotients           */
       **hf;    /* [0..n-1,0..n-1] array for the system matrix        */

/* ------------------------------------------------------------------ */

static int iter4
/*.IX{iter4}*/
                (
                 nlgls   fkt,
                 int     n,
                 REAL    epsm,
                 REAL    xalt[],
                 REAL    x1[],
                 boolean *sing
                )

/***********************************************************************
* Compute one approximation via  Brown's method.                       *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* fkt   Function, which evaluates the right hand side. It has          *
*       components from  0 to n - 1.                                   *
* n     number of equations                                            *
* epsm  machine constant                                               *
* xalt  previous iterate                                               *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x1    [0..n-1] vector, the new iterate for the zero                  *
* sing  error code, matrix is singular                                 *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok                                                          *
* = 1: error in evaluating  fkt()                                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ihf, rslin, dquot, hf, nlgls, REAL, boolean, subst, FABS, FALSE,     *
* TRUE, ZERO, copy_vector                                              *
***********************************************************************/

{
  int  i,       /* Loop counters                                      */
       j,       /*                                                    */
       k,       /* index of current function                          */
       anzahl,  /* counter for those difference quotients that are 0  */
       temp = 0,/* aux variable                                       */
       kmax = 0,/* Index of largest difference quotient               */
       jsub;    /* aux variable for an index in ihf                   */
  REAL hold,    /* aux variable for x1[temp]                          */
       h,       /* Step size for difference quotient in direction of  */
                /*  x1[temp]                                          */
       faktor,  /* quotient  h / hold                                 */
       dermax,  /* maximal difference quotient magnitude              */
       sum,     /* sum in rslin[k]                                    */
       f,       /* value of function k at  x1                         */
       fplus;   /* value of function k at                             */
                /* (x1[0],..., x1[temp]+h,...,x1[n-1])                */


  /* --------------------- initialize variables --------------------- */

  for (j = 0; j < n; j++)
    ihf[0][j] = j,
  copy_vector(x1, xalt, n);



  /* ----------- linearize the kth coordinate function -------------- */

  for (k = 0; k < n; k++)
  {

    anzahl = 0;
    faktor = (REAL)0.001;

    for (j = 0; j < 3; j++)
    {
      if (k > 0)
        subst(n, k, ihf, hf, rslin, x1);
      if ((*fkt)(k, x1, &f))
        return 1;


      /* --- find ith diskretization size and difference quotient --- */

      for (i = k; i < n; i++)
      {
        temp = ihf[k][i];
        hold = x1[temp];
        h    = faktor * hold;
        if (FABS(h) <= epsm)
          h = (REAL)0.001;
        x1[temp] = hold + h;
        if (k > 0)
          subst(n, k, ihf, hf, rslin, x1);
        if ((*fkt)(k, x1, &fplus))
          return 1;
        x1[temp] = hold;

        dquot[temp] = (fplus - f) / h;

        if (FABS(dquot[temp]) <= epsm)
          anzahl++;
        else if (FABS(f / dquot[temp]) >= (REAL)1.0e20)
          anzahl++;
      }


      if (anzahl < n - k)
      {
        *sing = FALSE;
        break;
      }
      else
      {
        *sing  =  TRUE;
        faktor *= TEN;
        anzahl =  0;
      }
    }


    if (*sing)
      break;

    else if (k < n - 1)
    {
      kmax = ihf[k][k];

      /* --- find largest magnitude difference quotient ------------ */

      for (dermax = FABS(dquot[kmax]), i = k + 1; i < n; i++)
      {
        jsub = ihf[k][i];
        if (FABS(dquot[jsub]) < dermax)
          ihf[k + 1][i] = jsub;
        else
          ihf[k + 1][i] = kmax,
          kmax          = jsub;
      }

      if (FABS(dquot[kmax]) <= epsm)
        *sing = TRUE;

      ihf[k][n] = kmax;


      if (*sing)
        break;

      else
      {

        /* ---------- solve kth equation for  xmax ------------------ */

        for (sum = ZERO, j = k + 1; j < n; j++)
          jsub        =  ihf[k + 1][j],
          hf[k][jsub] =  -dquot[jsub] / dquot[kmax],
          sum         += dquot[jsub] * x1[jsub];
        rslin[k] = (sum - f) / dquot[kmax] + x1[kmax];

      }
    }

    else
    {

      /* ----- solve (n-1)th coordinate function via discrete ------- */
      /* ----- Newton method for one variable ----------------------- */

      if (FABS(dquot[temp]) <= epsm)
        *sing = TRUE;
      else
        kmax     = temp,
        rslin[k] = -f / dquot[kmax] + x1[kmax];

    }
  }



  if (! *sing)
  {
    x1[kmax] = rslin[n - 1];                 /* compute approximation */
    if (n > 1)                               /* by back substitution  */
      subst(n, n - 1, ihf, hf, rslin, x1);
  }


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int brown   /* Brown's method for nonlinear systems of equations .....*/
/*.IX{brown}*/
         (
          nlgls   fkt,      /* Function ..............................*/
          int     n,        /* number of equations ...................*/
          REAL    x0[],     /* Starting value for iteration ..........*/
          REAL    eps,      /* error bound ...........................*/
          int     prot,     /* Protokol switch .......................*/
          int     maxit,    /* maximal number of steps ...............*/
          REAL    x1[],     /* solution ..............................*/
          int     *itanz    /* actual steps performed ................*/
         )                  /* error code ............................*/

/***********************************************************************
* find a zero of a nonlinear system of n equations in n unknown using  *
* Brown's method                                                       *
.BE*)
*                                                                      *
* One step of Brown's method is performed by calling iter4().          *
* The iteration is continued until the admissable number of iterations *
* is reached or one of the following break-off criteria is satisfied : *
* 1. the relative change between two successive iterates is less than  *
*    eps.                                                              *
* 2. The function value has mgnitude less than MACH_EPS at the new x   *
* 3. The desired accuracy has been reached.                            *
*                                                                      *
* If specified in prot, intermediate results are tabulated, see prot   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* fkt    Function, whichdescribes the right hand side of the system;   *
*        its components are numbered 0 to n-1.                         *
* n      size of system                                                *
* x0     [0..n-1] starting vector for the iteration                    *
* eps    desired accuracy; if specified as <= 0, we set                *
*        eps = 0.8 * MACH_EPS.                                         *
* prot   Protocol flag. If TRUE we tabulate the differenz of successive*
*        iterates, the last iterate and the functional value there     *
*        after each iteration in the standard output file. If FALSE,   *
*        there is no output.                                           *
* maxit  Maximal number of iterations                                  *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x1     [0..n-1] vector, the approximate solution for the zero        *
* itanz  number of iterations executed                                 *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* error code.                                                          *
* = 0: successsful iteration                                           *
* = 1: desired accuracy not achieved after maxit iterations            *
* = 2: system matrix singular                                          *
* = 3: lack of memory                                                  *
* = 4: wrong input : fkt = NULL or n < 1 or maxit < 1                  *
* = 5: error in evaluating  fkt()                                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ihf, rslin, dquot, hf, nlgls, REAL, boolean, MACH_EPS, FALSE,        *
* copy_vector, iter4, subst, vminit, vmalloc, vmcomplete, vmfree,      *
* VEKTOR, IMATRIX, MATRIX, printf, LZP, FABS, EIGHT, NULL              *
.BA*)
***********************************************************************/
/*.BE*/

{

  int     i,       /* Loop variable                                   */
          m,       /* iterations counter                              */
          krit;    /* list of causes for stop of program:             */
                   /* = 0: too many steps                             */
                   /* = 1: step size too small                        */
                   /* = 2: Function value sufficiently small          */
                   /* = 3: desired accuracy reached                   */
  boolean sing;    /* Flag, indicating whether the linearized system  */
                   /* is singular                                     */
  REAL    relf = ZERO,  /* aux variables                              */
          fwert,   /*                                                 */
          delta0,  /* maximum norm of previous vector                 */
          delta1,  /* maximum norm of current vector                  */
          epsm,    /* machine constant                                */
          *xalt;   /* [0..n-1] vector used to store old iterate       */
  void    *vmblock;     /* dynamically allocated vectors and matrices */

  /* --------------- catch wrong input ------------------------------ */

  if (fkt == NULL || n < 1 || maxit < 1)
    return 4;

  if (eps < MACH_EPS)            /* desired accuracy rediculously     */
    eps = EIGHT * MACH_EPS;      /* small ? correct                   */


  /* ------- allocate aux vectors, matrices ------------------------- */

  vmblock = vminit();
  xalt  = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  ihf   = (int **) vmalloc(vmblock, IMATRIX, n, n + 1);
  hf    = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  rslin = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  dquot = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  if (! vmcomplete(vmblock))                       /* lack of memory? */
  {
    vmfree(vmblock);
    return 3;
  }


  epsm   = MACH_EPS;                         /* initialize variables  */
  sing   = FALSE;
  krit   = 0;
  delta0 = (REAL)0.01;
  copy_vector(xalt, x0, n);


  if (prot)
    printf("\n");


  for (m = 1; m <= maxit; m++)                    /* start iteration  */
  {

    if (iter4(fkt, n, epsm, xalt, x1, &sing))
    {
      vmfree(vmblock);
      return 5;
    }


    /* ---------- if desired document each step --------------------- */

    if (prot)
    {
      printf("%3d. Iteration step\n", m);
      if (! sing)
      {
        printf("      Difference         Component       "
               "Approximation             Function value\n");
        for (i = 0; i < n; i++)
        {
          if ((*fkt)(i, x1, &fwert))
          {
            vmfree(vmblock);
            return 5;
          }
          printf("%22"LZP"e  %4d  %22"LZP"e  %22"LZP"e\n",
                 x1[i] - xalt[i], i, x1[i], fwert);
        }
      }
      else
        printf("Jacobi matrix singular!\n");
    }


    if (! sing)
    {

      /* ---------------- test break-off criteria ------------------- */

      for (i = 0; i < n; i++)  /* test relative change in new iterate */
      {
        relf = (x1[i] - xalt[i]) / (xalt[i] + eps);
        if (FABS(relf) >= eps)
          break;
      }
      if (FABS(relf) < eps)
      {
        krit = 1;
        break;
      }

      for (i = 0; i < n; i++)            /* check function value     */
      {
        if ((*fkt)(i, x1, &fwert))
        {
          vmfree(vmblock);
          return 5;
        }
        if (FABS(fwert) > epsm)
          break;
      }
      if (FABS(fwert) <= epsm)
      {
        krit = 2;
        break;
      }

      delta1 = FABS(x1[0] - xalt[0]);         /* compare with desired */
      for (i = 1; i < n; i++)                 /* accuracy             */
        if (delta1 < FABS(x1[i] - xalt[i]))
          delta1 = FABS(x1[i] - xalt[i]);
      if (delta1 <= (REAL)0.001)
        if (delta0 <= delta1)
        {
          krit = 3;
          break;
        }
      delta0 = delta1;


      if (m < maxit)
        copy_vector(xalt, x1, n);
    }

    else
      break;

  }


  *itanz = m;


  vmfree(vmblock);
  if (sing)
    return 2;
  else if (krit == 0)
    return 1;
  else
    return 0;
}

/* --------------------------- END brown.c -------------------------- */
