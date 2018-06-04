#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE implruku.c ----------------------- */

/***********************************************************************
*                                                                      *
* Solve a first order system of ordinary differential equations using  *
* -------------------------------------------------------------------  *
* implicit Runge-Kutta methods                                         *
* ----------------------------                                         *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Thomas Eul                                     *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Literature:           D. Sommer: Neue implizite Runge-Kutta-Formeln  *
*                       und deren Anwendungsmoeglichkeiten, Ph. D.     *
*                       thesis, Aachen 1967,   [SOMM67]                *
* Date:                 6.2.1992 ; 11,3 1995                           *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  SQRT, POW, MACH_EPS, sqr, dglsysfnk, */
                        /*       FALSE, TRUE, SEEK_END, REAL, LZS,    */
                        /*       ZERO, HALF, ONE, TWO, FIVE, TEN, LZP */
#include <vmblock.h>    /*  for  vmalloc, vmcomplete, vmfree, vminit, */
                        /*       VEKTOR, MATRIX                       */
#include <legendre.h>   /*  for  legendre                             */
#include <implruku.h>   /*  for  implruku                             */

/* put out a real vector of length `n' (and its name)                 */
#define zeig(v, n)                   \
  {                                  \
    int i;                           \
    printf("%-9s", #v": ");          \
    for (i = 0; i < n; i++)          \
      printf("%30.19"LZP"e", v[i]);  \
    printf("\n");                    \
  }



/* ------------------------------------------------------------------ */

#define FLTSIZE     (sizeof(REAL))                /* for abbreviation */
/*.IX{FLTSIZE}*/



/* ------------------------------------------------------------------ */

static int kopf_erzeugen       /* Generate headings ..................*/
/*.IX{kopf\unt erzeugen}*/
                        (
                         int  n,           /* number of DEs ..........*/
                         int  mmax,        /* max. order of IRKMs  ...*/
                         REAL xend,        /* desired final x-value ..*/
                         REAL x0,          /* initial x-value ........*/
                         REAL y0[],        /* initial y-value at x0 ..*/
                         REAL g[],         /* weights for  y .........*/
                         REAL eps_rel,     /* relative error bound ...*/
                         char file_out[],  /* name of output file ....*/
                         char file_pro[],  /* name of protocol file ..*/
                         FILE **ausdat,    /* output file ............*/
                         FILE **prodat     /* protocol file ..........*/
                        )                  /* error code .............*/

/***********************************************************************
* Generate headings for output and protocol files, if those are wanted.*
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* (see implruku())                                                     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ausdat  potential output file                                        *
* prodat  potential protocol file                                      *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 4: error opening protocol  or output file                          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, fopen, fprintf, NULL, FILE, LZS, LZP                           *
***********************************************************************/

{
  int i;


  /* ------- prepare heading for output file ------------------------ */

  if (*file_out != '\0')                 /* does name of file exist?  */
  {
    if ((*ausdat = fopen(file_out, "w")) == NULL)
      return 4;

    fprintf(*ausdat, " Initial values:\n\n"
                     "          x0            comp.            y0\n"
                     "%23.15"LZP"e   1   %23.15"LZP"e", x0, y0[0]);

    for (i = 1; i < n; i++)
      fprintf(*ausdat, "\n                       %4d   %23.15"LZP"e",
                       i, y0[i]);

    fprintf(*ausdat, "\n\n"
                     " Final x-value:    %23.15"LZP"e\n"
                     " Desired accuracy: %23.15"LZP"e\n"
                     " Maximal order:    %2d\n\n"
                     " Name of protocol file: %s\n\n"
                     " step  final x-value          comp."
                     "    approximate solution    error estimate",
                     xend, eps_rel, mmax, file_pro);
  }
  else
    *ausdat = NULL;


  /* ------ generate heading for protocol file ---------------------- */

  if (*file_pro != '\0')                  /* does file name exist?    */
  {
    if ((*prodat = fopen(file_pro, "w")) == NULL)
      return 4;

    fprintf(*prodat, " Initial values :\n\n"
                     "          x0            comp.            y0\n"
                     "%23.15"LZP"e   1   %23.15"LZP"e", x0, y0[0]);

    for (i = 1; i < n; i++)
      fprintf(*prodat, "\n                       %4d   %23.15"LZP"e",
                       i, y0[i]);

    fprintf(*prodat, "\n\n"
                     " Final x-value:    %23.15"LZP"e\n"
                     " Desired accuracy: %23.15"LZP"e\n"
                     " Maximal order:    %2d\n\n"
                     " Reason for decreasing step size:\n"
                     "     0  no decrease of step size\n"
                     "     1  e_rel >= eps_rel\n"
                     "     2  delta_g >= dk\n"
                     "     3  dg_rel >= eps_rel",
                     xend, eps_rel, mmax);

    fprintf(*prodat, "\n\n Weights G:       comp.            G\n");
    for (i = 0; i < n; i++)
      fprintf(*prodat, "                  %3d   %23.15"LZP"e\n",
                       i, g[i]);

    fprintf(*prodat, "\n\n"
                     " Name of output file: %s\n\n"
                     " step ord-  step size    upper bound      "
                     "approxim.  error      # fct.  #    cau"
                     "\n      er          h      of interval  "
                     "       Y       estimate   calls   iter se",
                     file_out);
  }
  else
    *prodat = NULL;

  return 0;
}



/* ------------------------------------------------------------------ */

static int machstuetz   /* Create file with implicit R-K coefficients */
/*.IX{machstuetz}*/
                     (
                      int  mmax,        /* max. order of coefficients */
                      char file_st[]    /* name of file ..............*/
                     )

/***********************************************************************
* This subroutine computes the coefficients for implicit Runge-Kutta   *
* Methods (IRM) from order 1 to mmax, if a node file does not exist    *
* or if it does not contain enough nodes.                              *
* The result is written unformatted onto the file file_st. This data   *
* can be read and used in implruku().                                  *
* For each order m we first compute the Gauss-Legendre nodes alpha[j], *
* j=1, ..., m, which are the roots of the Legendre polynomials         *
* transformed to the interval [0,1].                                   *
* The coefficients beta[i][j] and A[j], i,j=1,...,m, are given by the  *
* solution of m*(m+1) linear systems of equations. These solutions can *
* be found by multiplying the Lagrange polynomials.                    *
* The node file has the form:                                          *
*   1         coefficient   alpha  \                                   *
*   1         coefficient   beta    > of order 1                       *
*   1         coefficient   A      /                                   *
*   2         coefficients  alpha  \                                   *
*   4         coefficients  beta    > of order 2                       *
*   2         coefficients  A      /                                   *
*   ...                                                                *
*   mmax      coefficients  alpha  \                                   *
*   mmax*mmax coefficients  beta    > of order  mmax                   *
*   mmax      coefficients  A      /                                   *
* The matrix beta is filled rowwise.                                   *
* The file contains   3+8+...+mmax*(mmax+2)  REAL numbers.             *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* mmax     maximal order up to which we want to produce implicit R-K   *
*          coefficients                                                *
* file_st  Name of file with nodes for implruku()                      *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 3: lack of available memory                                        *
* = 5: error opening node file                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, vminit, vmalloc, vmcomplete, vmfree, VEKTOR, MATRIX, legendre, *
* fopen, NULL, fseek, SEEK_END, ftell, fclose, fwrite, ZERO, HALF,     *
* ONE                                                                  *
***********************************************************************/

{
  void *vmblock;   /* List of dynamic allocations                     */
  REAL *c,         /* [0..mmax-1] vector that stores the coefficients */
                   /* of the Lagrange polynomials                     */
       *alpha,     /* [0..mmax-1] vector            \  arrays for the */
       **beta,     /* [0..mmax-1,0..mmax-1] matrix   > coefficients of*/
       *A,         /* [0..mmax-1] vector            /  the impl. R-K  */
       zj,         /* numerator of the jth Lagrange polynomial        */
       betajk,     /* aux variable for finding beta[j][k] by multi-   */
                   /* plying the Lagrange polynomials                 */
       alphak,     /* aux storage for alpha[k]                        */
       Aj;         /* aux variable for A[j]                           */
  int  m,          /* current order of coefficients being computed    */
       i, j, k,    /* Loop variables                                  */
       ng;         /* counter for number of factors                   */
                   /* (alpha[k]-alpha[j]), k=1,...,j-1,j+1,...,m-1,   */
                   /* which are multiplied                            */
  FILE *stuetzdat; /* node file                                       */


  if ((stuetzdat = fopen(file_st, "r")) != NULL)     /* does the node */
  {                                                  /* file exist?   */
    fseek(stuetzdat, 0l, SEEK_END);             /* go to end of file  */
    if (ftell(stuetzdat) >=                     /* file large enough? */
        (long)((mmax * (mmax + 1) *
                (2 * mmax + 7) / 6) * FLTSIZE)
       )
    {
      fclose(stuetzdat);                    /* use current node file  */
      return 0;
    }
    fclose(stuetzdat);
  }


  /* ------- allocate storage --------------------------------------- */

  vmblock = vminit();                           /* initialize storage */
  c     = (REAL *) vmalloc(vmblock, VEKTOR, mmax, 0);
  A     = (REAL *) vmalloc(vmblock, VEKTOR, mmax, 0);
  alpha = (REAL *) vmalloc(vmblock, VEKTOR, mmax, 0);
  beta  = (REAL **)vmalloc(vmblock, MATRIX, mmax, mmax);

  if (! vmcomplete(vmblock))         /* lack of memory?               */
  {
    vmfree(vmblock);                 /* free buffers and report error */
    return 3;
  }


  if ((stuetzdat = fopen(file_st, "wb")) == NULL)      /* create node */
  {                                                    /* file for    */
    vmfree(vmblock);                                   /* writing in  */
    return 5;                                          /* binary mode */
  }

  for (m = 1; m <= mmax; m++)  /* generate nodes alpha[j] and weights */
  {                            /* beta[j][k] and A[j] for the orders  */
                               /* from 1 to mmax; store unformatted   */
                               /* in node file                        */

    legendre(m, alpha);          /* compute all roots of the Legendre */
                                 /* polynomial of degree m in the     */
                                 /* interval [-1;1]                   */

    for (i = 0; i < m; i++)                /* transform alpha[i] into */
      alpha[i] = HALF * alpha[i] + HALF;   /* the interval [0;1]      */

    for (j = 0; j < m; j++)             /* compute weights beta[j][k] */
    {                                   /* and A[j]                   */
      zj = ONE;
      for (k = 0; k < j; k++)
        zj *= alpha[j] - alpha[k];
      for (k = j + 1; k < m; k++)
        zj *= alpha[j] - alpha[k];

      for (c[0] = ONE, ng = 0, k = 0;  /* compute coefficients of the */
           k < j;                      /* (j+1)th Lagrange polynomial */
           k++, ng++
          )
      {
        for (alphak = -alpha[k], i = ng; i >= 0; i--)
          c[i + 1] = c[i];
        for (c[0] = alphak * c[1], i = 1; i <= ng; i++)
          c[i] += alphak * c[i + 1];
      }
      for (k = j + 1; k < m; k++, ng++)
      {
        for (alphak = -alpha[k], i = ng; i >= 0; i--)
          c[i + 1] = c[i];
        for (c[0] = alphak * c[1], i = 1; i <= ng; i++)
          c[i] += alphak * c[i + 1];
      }

      for (zj = ONE / zj, Aj = ZERO, k = 0;       /* compute A[j] and */
           k < m;                                 /* all beta[j][.]   */
           k++
          )
      {
        for (alphak = alpha[k], betajk = ZERO, i = m; i >= 1; i--)
          betajk = alphak * (c[i - 1] / i + betajk);
        beta[j][k] =  betajk * zj;
        Aj         += c[k] / (k + 1);
      }
      A[j] = Aj * zj;
    }

    fwrite(alpha, FLTSIZE, m, stuetzdat);       /* store transformed  */
    for (i = 0; i < m; i++)                     /* nodes in node file */
      fwrite(beta[i], FLTSIZE, m, stuetzdat);
    fwrite(A, FLTSIZE, m, stuetzdat);
  }

  fclose(stuetzdat);


  vmfree(vmblock);
  return 0;
}



/* ------------------------------------------------------------------ */

static int vektor_null                /* Is `vektor' the zero vector? */
                      (
                       int  n,                    /* length of vector */
                       REAL vektor[]              /* [0..n-1] vector  */
                      )

{
  while (n-- != 0)
    if (*vektor++ != ZERO)
      return FALSE;

  return TRUE;
}



/* ------------------------------------------------------------------ */

static int matrix_null                /* Is `matrix' the zero matrix? */
                      (
                       int  n,              /* size of matrix         */
                       REAL *matrix[]       /* [0..n-1,0..n-1] matrix */
                      )

{
  int  i,                        /* row counter                       */
       j;                        /* column counter                    */
  REAL *zeile;                   /* address of current matrix element */

  for (i = n; i != 0; i--)
    for (zeile = *matrix++, j = n; j != 0; j--)
      if (*zeile++ != ZERO)
        return FALSE;

  return TRUE;
}



/* ------------------------------------------------------------------ */

static int init_irkv        /* Initialize implruku() .................*/
/*.IX{init\unt irkv}*/
                    (
                     int  n,            /* number of DEs .............*/
                     int  mmax,         /* max. order of IRKMs (>= 5) */
                     char file_st[],    /* name of coefficient file ..*/
                     REAL *epsm,        /* machine constant ..........*/
                     REAL eps_rel,      /* relative error bound or ...*/
                     REAL g[],          /* weights for y .............*/
                     REAL xend,         /* desired final x-value .....*/
                     REAL x0,           /* initial x-value ...........*/
                     REAL y0[],         /* initial y-value at x0 .....*/
                     char file_out[],   /* name of output file .......*/
                     char file_pro[],   /* name of protocol file .....*/
                     FILE **ausdat,     /* output file ...............*/
                     FILE **prodat,     /* protocol file .............*/
                     REAL *summeg       /* sum of weights ............*/
                    )                   /* error code ................*/

/***********************************************************************
* Make preparations for implruku():                                    *
* - check input data of implruku(),                                    *
* - if needed form node file,                                          *
* - if needed compute new value for machine constant epsm,             *
* - initialization of summeg,                                          *
* - if needed make headings in output and protocol files               *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* (see implruku())                                                     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* epsm    possible adjusted value for machine constant                 *
* ausdat  possibly open output file                                    *
* prodat  possibly open protocol file                                  *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 2: wrong input parameter(s)                                        *
* = 4: error opening protocol or output file                           *
* = 5: error opening node file                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* machstuetz, kopf_erzeugen, REAL, FALSE, TRUE, FILE, ZERO,            *
* vektor_null                                                          *
***********************************************************************/

{
  int fehler,                          /* error code of this function */
      i;                               /* loop variable               */


  fehler = 0;                              /* assume valid input data */

  if (n       <     1 ||                                /* n,         */
      mmax    <     5 ||                                /* mmax,      */
      eps_rel <= ZERO ||                                /* eps_rel or */
      vektor_null(n, g)                                 /* g invalid? */
     )
    fehler = 2;
  else                                    /* input parameters valid?  */
    fehler = machstuetz(mmax, file_st);   /* open node file if needed */


  if (! fehler)                           /* no error found so far?   */
  {
    if (*epsm <= ZERO)                    /* improper value for epsm? */
      *epsm = MACH_EPS;                   /* correct it, please       */

    for (*summeg = ZERO, i = 0; i < n; i++)       /* precalculate the */
      *summeg += g[i];                            /* sum of weights   */

    fehler = kopf_erzeugen(n, mmax, xend,            /* initialize    */
                           x0, y0, g, eps_rel,       /* output and    */
                           file_out, file_pro,       /* protocol file */
                           ausdat, prodat            /* if needed     */
                          );
  }


  return fehler;
}


/* ------------------------------------------------------------------ */

static REAL wurzel   /* Compute a square root as needed in implruku() */
/*.IX{wurzel}*/
                  (
                   int  n,        /* number of DEs ...................*/
                   REAL *dfdy[],  /* derivative (df/dy)(x0,y0) .......*/
                   REAL dfdx[],   /* 2nd factor in scalar prod. <g,.> */
                   REAL summeg,   /* sum of weights ..................*/
                   REAL g[],      /* weights for y ...................*/
                   REAL dfdxh[]   /* aux vector ......................*/
                  )               /* square root .....................*/

/***********************************************************************
* Compute the square root                                              *
*         sqrt(summe(g[i] * (df/dx * (df/dy)^k)[i]^2) / summeg)        *
* where the first call is executed with k = 1, the second call with    *
* k = 2, etc. This is achieved by storing the vector                   *
* df/dx * (df/dy)^(k-1)  in dfdx in each call for later use.           *
* These roots are used in implruku() to determine the step size and    *
* mean estimated iteration errors dk and deltak[l], l=0,1,2*m.         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n       number of differential equations                             *
* dfdy    [0..n-1,0..n-1] matrix with the derivatives df/dy of the     *
*         right hand side of the DE system wrt. y at (x0,y0), i.e.,    *
*         dfdy[i][j] = df[i]/dy[j]                                     *
* dfdx    [0..n-1] vector with the derivative df/dx of the right hand  *
*         side of the DE system wrt. x (1st. call) respectively        *
*         df/dx * (df/dy)^(k-1)  at (x0,y0) (kth call)                 *
* g       [0..n-1] weight vector                                       *
* summeg  g[0] + g[1] + g[2] + ... + g[n-1]                            *
* dfdxh   [0..n-1] aux vector                                          *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* dfdx    [0..n-1] vector with the new powers  df/dx * (df/dy)^k       *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* square root                                                          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, sqr, SQRT, ZERO                                                *
***********************************************************************/

{
  REAL summe;                     /* scalar products <dfdy[i,.],dfdx> */
                                  /* and             <dfdxh,g>        */
  REAL tmp;                       /* Function value                   */
  REAL min;                       /* smallest positive machine number */
  int  i, j;                      /* loop variables                   */


  for (i = 0; i < n; i++)
  {
    for (summe = ZERO, j = 0; j < n; j++)
      summe += dfdx[j] * dfdy[i][j];
    dfdxh[i] = summe;
  }

  for (summe = ZERO, j = 0; j < n; j++)
    dfdx[j] =  dfdxh[j],
    summe   += sqr(dfdx[j]) * g[j];


  min = POSMIN;
  tmp = SQRT(summe / summeg);
#ifdef DEBUGg
  printf("min = %30.19"LZP"e\n", min);
  printf("tmp = %30.19"LZP"e\n", tmp);
  if (min == ZERO)
    printf("min is zero!\n");
#endif
  if (FABS(tmp) < min)              /* less than smallest number?     */
    tmp = ZERO;                     /* impossible!                    */
  return tmp;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int implruku  /* Implic. Runge-Kutta meth. for ODE syst. of 1st order */
/*.IX{implruku}*/
            (
             dglsysfnk dgl,          /* right hand side of DE system .*/
             int       n,            /* number of DEs ................*/
             int       mmax,         /* max. order of IRKMs (>= 5) ...*/
             char      file_st[],    /* name of coefficient file .....*/
             char      file_out[],   /* name of output file ..........*/
             char      file_pro[],   /* name of protocol file ........*/
             REAL      epsm,         /* machine constant .............*/
             REAL      *eps_rel,     /* relative error bound or       */
                                     /* maximal relative error .......*/
             long      fmax,         /* max. number of calls of dgl() */
             long      *aufrufe,     /* act. number of calls of dgl() */
             REAL      g[],          /* weights for y ................*/
             REAL      *x0,          /* initial/final x-value ........*/
             REAL      xend,         /* desired final x-value ........*/
             REAL      y0[],         /* initial y-value at x0 ........*/
             REAL      yq[]          /* final y-value at xend ........*/
            )                        /* error code ...................*/

/***********************************************************************
* The subroutine implruku() solves the initial value problem           *
*                 y' = f(x,y),     y(x0) = y0                          *
* for a system of n first order ordinary differential equations.       *
* The solution is found using the implicit Runge-Kutta methods (IRKM). *
* We use step size control and adjust the order of the R-K method.     *
* For this implruku() needs the implicit R-K formulas of orders 1 to   *
* mmax. The subroutines expect to find the coefficients for these      *
* methods in the file named file_st. If this file is not available or  *
* insufficient for mmax, it will be generated internally.              *
.BE*)
*                                                                      *
* To find the optimal order, we check on the anticipated amount of     *
* work in evaluating the right hand side according to the formula      *
*    AW(eps,m)=(n+1+4*m^2)/h(eps,m)   before each Runge-Kutta step.    *
* Here h(eps,m) denotes a theoretical step size for the order m and    *
* the accuracy bound eps. For each step we select the implicit R-K     *
* methods of orders m and m+1 , so that:                               *
*         AW(eps,i) >  AW(eps,i+1)       for  i=1,...,m-1   and        *
*         AW(eps,m) <= AW(eps,m+1);                                    *
* except we set m = mmax - 1, when no such m exists between 1 and      *
* mmax - 1.                                                            *
*                                                                      *
* For the selected m we then choose a step size h, which is            *
* essentially the one from theory.                                     *
* Following this, we execute the Runge-Kutta step twice. The implicit  *
* R-K method of order m has the coefficients  A, alpha and beta, the   *
* one of order m+1 has the coefficients Aq, alphaq and betaq.          *
* The ki are computed iteratively.                                     *
*                                                                      *
* If the mean absolute difference e_rel exceeds eps_rel for the two    *
* approximations y and yq, then we set the step size to be             *
*         h * sf * (0.5 * eps_rel / e) ^ (1 / (2 * m + 1)),            *
* as warranted from theory. Here sf is a fudge factor, e <= e_rel.     *
* Then the last step is repeated for the new h.                        *
*                                                                      *
* If the mean absolute difference of two consecutive approximations    *
* yq(alt) and yq(neu) exceeds the mean theoretical error estimate dk   *
* in one iteration, convergence can not be expected and we reduce h by *
* 6/10 and start anew.                                                 *
*                                                                      *
* We iterate until the relative mean absolute difference of two        *
* approximations becomes less than eps_rel. This ought to happen       *
* theoretically after at most 2*m+1 iterations.                        *
*                                                                      *
* If this does not happen numerically, we increase the order m by 1    *
* and repeat the last step. If this would raise the order beyond       *
* mmax - 1,  we are only left to try to reduce the step size such as   *
* to  0.8 * h.  If the step size then decreases to below 10 * machine  *
* constant, we stop integrating: the problem will not converge on the  *
* given machine.                                                       *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* dgl       pointer to a function which evaluates the right hand side  *
*           of the DE system   y' = f(x,y).                            *
*           The pointer must be of the type `dglsysfnk':               *
*               void (*dglsysfnk)(REAL x, REAL y[], REAL f[]);         *
*           That function must contain statements similar to the       *
*           following:                                                 *
*               ...                                                    *
*               f[0]   := F1(x, y[0], y[1],..., y[n-1]);               *
*               f[1]   := F2(x, y[0], y[1],..., y[n-1]);               *
*               ...       ...                   ...                    *
*               f[n-1] := Fn(x, y[0], y[1],..., y[n-1]);               *
*               ...                                                    *
*           where F1,...,Fn are the n function expressions on the      *
*           right hand side of the explicit DE system.                 *
* n         number of DEs                                              *
* mmax      maximal order allowed. mmax >= 5                           *
*           The maximal order should not be raised too much since the  *
*           quality of the results is linked to the machine constant:  *
*           For the FORTRAN version of this program and the            *
*           Control Data Cyber 175 with a machine constant of ca.      *
*           2.5e-29, the best overall maximal order turns out to be    *
*           mmax = 12. For particular examples a larger value for mmax *
*           may be useful.                                             *
* file_st   file_st is the file with the coefficients for the          *
*           implicit R-K methods up to order mmax. If this file is not *
*           available or does not extend to mmax, it will be generated *
*           internally.                                                *
* file_out  If intermediate output is desired, the file file_out is    *
*           created. If file_out is is an empty string, there is no    *
*           intermediate output.                                       *
* file_pro  protocol file; if file_pro is an empty string, the         *
*           protocol is not kept.                                      *
* epsm      machine constant. If epsm <= 0 we replace this error with  *
*           the correct value. epsm is needed for approximation of     *
*           partial derivatives dfdx and dfdy and for detection of too *
*           small step size.                                           *
* eps_rel   desired relative accuracy                                  *
* fmax      upper bound for the number of allowed function evaluations *
*           of the right hand side via dgl()                           *
* g         [0..n-1] vector with weights for y. The weights g[i] allow *
*           different emphasis of the components y[i] wrt. the         *
*           accuracy bound eps_rel. If all components are to carry the *
*           same weight, choose  g[i] = 1 for all i .                  *
*           W A R N I N G : If g[i] = 0 for one i, we might encounter  *
*           division by zero if the corresponding component of the     *
*           partial derivative of the right hand side is also 0. We    *
*           have not built in any safeguards against this!             *
* x0        initial x-value                                            *
* xend      desired final x-value                                      *
* y0        [0..n-1] vector of initial y-values y(x0) at x0            *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* aufrufe  actual number of calls of dgl() (see fmax)                  *
* eps_rel  estimate for maximal local relative error of successful     *
*          steps                                                       *
* x0       final x-value of integration;                               *
*          normally (fehler = 0) x = xend                              *
* yq       [0..n-1] vector of final y-values at xend for method of     *
*          order m+1                                                   *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 1: no convergence. Possible remedy: increase mmax                  *
* = 2: wrong input parameter(s)                                        *
* = 3: lack of memory                                                  *
* = 4: error opening protocol or output files                          *
* = 5: error opening node file                                         *
* = 6: too many calls of dgl() (aufrufe > fmax)                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* init_irkv, wurzel, REAL, FALSE, TRUE, dglsysfnk, SQRT, POW,          *
* vminit, vmalloc, vmcomplete, vmfree, VEKTOR, MATRIX, fopen, fseek,   *
* SEEK_SET, fread, fprintf, fclose, LZS, ZERO, HALF, ONE, TWO, FIVE,   *
* TEN, LZP, matrix_null, vektor_null                                   *
.BA*)
***********************************************************************/
/*.BE*/

{
  void      *vmblock;     /* List of dynamic allocations              */
  REAL      *fak,         /* [0..mmax-1] vector with factorials       */
                          /* fak[i-1] = (2 * i)!                      */
            *f0,          /* [0..n-1] vector of f-values of the right */
                          /* hand side at (x0,y0) preceeding every    */
                          /* step                                     */
            *f1,          /* [0..n-1] vector of f-values of the right */
                          /* hand side at intermediate points and in  */
                          /* derivative estimates                     */
            *dfdx,        /* [0..n-1] vector with derivatives df/dx   */
                          /* for right hand side                      */
            *deltak,      /* [0..2*mmax] vector for mean iteration    */
                          /* error estimate dk at each step           */
            *heps,        /* [0..mmax-1] vector: heps[m-1] is the     */
                          /* theoretical step size for the implicit   */
                          /* R-K method of order m for eps_rel        */
            *y,           /* [0..n-1] approximate solution vector for */
                          /* order m                                  */
            *yalt,        /* [0..n-1] old approximate solution vector */
                          /* yq kept for estimating accuracy          */
            **dfdy,       /* [0..n-1,0..n-1] derivative matrix for    */
                          /* df/dy:  dfdy[i][j] = df[i]/dy[j]         */
            **ki,         /* [0..n-1,0..mmax-2] matrix for ki of the  */
                          /* implicit R-K method of order m           */
            **kiq,        /* [0..n-1,0..mmax-1] matrix, ditto for     */
                          /* order m+1                                */
            **db,         /* [0..n-1,0..mmax-1] aux matrix for ki and */
                          /* kiq                                      */
            *alpha,       /* [0..mmax-2] vector          \  coeffic.  */
            **beta,       /* [0..mmax-2,0..mmax-2] matrix > of IRM    */
            *A,           /* [0..mmax-2] vector          /  of        */
                          /*                                 order m  */
            *alphaq,      /* [0..mmax-1] vector          \  ditto     */
            **betaq,      /* [0..mmax-1,0..mmax-1] matrix > for       */
            *Aq;          /* [0..mmax-1] vector          /  order m+1 */
  int       fehler,       /* error code of this function              */
            sz,           /* counter for number of intermediate steps */
            m,            /* current R-K method order                 */
            malt,         /* order used in previous step. If m=malt,  */
                          /* the coefficients need not be changed.    */
            mm2,          /* aux variable m*2                         */
            km2,          /* aux variable k*2                         */
            z,            /* used for dynamic adjustments in the      */
                          /* safety factor. z counts the number of    */
                          /* step size reductions due to  e >= eps.   */
            i, j, k,      /* Loop variables                           */
            l,            /* iteration counter                        */
            lh;           /* aux variable for output of l             */
  short int stopp,        /* flag indicating stop of integration      */
            vz,           /* sign of h, determines direction of       */
                          /* integration                              */
            null,         /* Flag:                                    */
                          /*  = 1:  df/dx = 0  or  df/dy = 0          */
                          /*  = 0:  otherwise                         */
            ord_up,       /* flag indicating change in order used.    */
                          /*  = 1:  raise the order for next step     */
                          /*  = 0:  no change in order                */
            cause;        /* indicates reason for change in step size */
                          /* = 0: no change in h                      */
                          /* = 1: e_rel   >= eps                      */
                          /* = 2: delta_g >= dk                       */
                          /* = 3: dg_rel  >= eps                      */
  float     aw1, aw2;     /* variables to estimate amount of work     */
                          /*     AW(eps,m) = (n+1+4*m*m)/h(eps,m)     */
  REAL      epserr,       /* epserr = 10 * epsm                       */
            epslok,       /* estimate for maximal local error         */
            delta,        /* square root of machine constant to       */
                          /* estimate derivatives using forward       */
                          /* differences                              */
            sf,           /* safety factor to reduce theoretical step */
                          /* size a bit.                              */
            x,            /* upper limit for interval of integration  */
                          /* from x0 to x0+h (=x)                     */
            summeg,       /* g[0] + g[1] + ... + g[n-1]               */
            yj,           /* aux variable when estimating the         */
                          /* partials df/dy[j]                        */
            yi,           /* scalar product y[i]  =                   */
                          /*         <ki[i,.],beta[.,k]>    resp.     */
                          /*         <a,ki[i,.]>                      */
            yqi,              /* scalar product yq[i] =               */
                          /*         <kiq[i,.],betaq[.,k]>  resp.     */
                          /*         <aq,kiq[i,.]>                    */
            sum,          /* aux variable for computing df/dx         */
            dk,           /* aux variable for deltak[l], later mean   */
                          /* absolute estimated iteration error.      */
                          /* This is given theoretically in the lth   */
                          /* iteration as                             */
                          /*   h^(l+2)/(l+2)! *                       */
                          /*   sqrt(g*(df/dx*(df/dy)^l)^2/summeg)*5,  */
                          /* where 5 acts as a safety factor          */
            zpot,         /* powers of ten                            */
            h,            /* step size                                */
            faklp1,       /* (l + 1)!                                 */
            x1,           /* nodes x0+h*alpha[k] in the               */
                          /* interval [x0,x]                          */
            yqnorm,       /* L2 norm of yq                            */
            y_diff,       /* aux variable for computing               */
                          /* mean difference                          */
            delta_g,      /* mean difference of two approximations    */
                          /* yq(alt) and yq(neu)                      */
            dg_rel,       /* relative mean difference of two          */
                          /* successive approximations  yq(alt) and   */
                          /* yq(neu)                                  */
            e,            /* mean difference of the two               */
                          /* approximations y and yq resulting from   */
                          /* two methods with different orders        */
            e_rel,        /* relative mean difference for y and yq    */
            genauigkeit,  /* error estimate for approximate solution  */
            hbeg = ZERO,  /* old step size                            */
            hmf0;         /* h * f0[i]                                */
  FILE      *stuetzdat,   /* node file                                */
            *ausdat,      /* output file                              */
            *prodat;      /* protocol file                            */

/* ---- define macro for output of the protocol file ---------------- */

#define PROTOKOLLIEREN(h, l)                                  \
    if (prodat != NULL)                                       \
      {                                                       \
        fprintf(prodat, "\n%4d  %3d  %12.4"LZP"e"             \
                        " %12.4"LZP"e   %12.4"LZP"e"          \
                        "  %10.3"LZP"e %4lu   %4d"            \
                        "   %1d", sz, m, h, x, yq[0],         \
                                  genauigkeit, *aufrufe, l,   \
                                  cause);                     \
        for (i = 1; i < n; i++)                               \
          fprintf(prodat, "\n%51.4"LZP"e", yq[i]);            \
      }


  /*********************************************************************
  *                     P r e p a r a t i o n s                        *
  *********************************************************************/

  if ((fehler = init_irkv(n, mmax, file_st, &epsm, *eps_rel, g,
                          xend, *x0, y0, file_out, file_pro,
                          &ausdat, &prodat, &summeg)) != 0)
    return fehler;

  if ((stuetzdat = fopen(file_st, "rb")) == NULL)   /* open node file */
    return 4;

  vmblock = vminit();                           /* initialize storage */

  /* ---- allocate space for vectors globally for the module  ------- */

  fak    = (REAL *)vmalloc(vmblock, VEKTOR, mmax,     0);
  f0     = (REAL *)vmalloc(vmblock, VEKTOR, n,        0);
  f1     = (REAL *)vmalloc(vmblock, VEKTOR, n,        0);
  dfdx   = (REAL *)vmalloc(vmblock, VEKTOR, n,        0);
  deltak = (REAL *)vmalloc(vmblock, VEKTOR, 2 * mmax, 0);
  heps   = (REAL *)vmalloc(vmblock, VEKTOR, mmax,     0);
  y      = (REAL *)vmalloc(vmblock, VEKTOR, n,        0);
  yalt   = (REAL *)vmalloc(vmblock, VEKTOR, n,        0);
  alpha  = (REAL *)vmalloc(vmblock, VEKTOR, mmax - 1, 0);
  A      = (REAL *)vmalloc(vmblock, VEKTOR, mmax - 1, 0);
  alphaq = (REAL *)vmalloc(vmblock, VEKTOR, mmax,     0);
  Aq     = (REAL *)vmalloc(vmblock, VEKTOR, mmax,     0);

  /* --- allocate space for matrices globally for the module -------- */

  dfdy  = (REAL **)vmalloc(vmblock, MATRIX, n,    n);
  ki    = (REAL **)vmalloc(vmblock, MATRIX, n,    mmax - 1);
  kiq   = (REAL **)vmalloc(vmblock, MATRIX, n,    mmax);
  db    = (REAL **)vmalloc(vmblock, MATRIX, n,    mmax);
  beta  = (REAL **)vmalloc(vmblock, MATRIX, mmax, mmax);
  betaq = (REAL **)vmalloc(vmblock, MATRIX, mmax, mmax);

  if (! vmcomplete(vmblock))                      /* lack of storage? */
  {
    vmfree(vmblock);               /* free buffers and report failure */
    return 3;
  }

  for (fak[0] = TWO, j = 4, i = 1;     /* initialize factorial vector */
       i < mmax;
       i++, j += 2                     /* fak[i] = (2*i)!             */
      )
    fak[i] = fak[i - 1] * j * (j - 1);

  *aufrufe = 0l;
  sz       = 0;

  m        = 0;
  malt     = 1;
  mm2      = 2;

  epserr   = TEN * epsm;
  epslok   = ZERO;
  delta    = SQRT(epsm);
  x        = xend - *x0;
  vz       = (x >= ZERO) ? 1 : -1;
  sf       = (REAL)0.9;


  /*********************************************************************
  *                      I n t e g r a t i o n                         *
  *********************************************************************/

  do                         /* perform one integration step per loop */
  {
    sz++;
    malt  = m;
    stopp = FALSE;

    if (*aufrufe > fmax)      /* too many function evaluations?    */
      stopp  = TRUE,          /* report error                      */
      fehler = 6;

#ifdef OHNE_ROMBERG
    (*dgl)(*x0, y0, f0);      /* find approximate partial derivatives */
    for (j = 0; j < n; j++)   /* of the right hand side wrt. y at     */
    {                         /* (x0,y0) using forward difference     */
      yj    = y0[j];          /* quotients                            */
    (*dgl)(*x0, y0, f0);           /* Naeherungen fuer die partiellen */
    for (j = 0; j < n; j++)        /* Ableitungen der rechten Seite   */
    {                              /* bezueglich y an der Stelle      */
      yj    = y0[j];               /* (x0,y0) mit Hilfe des vorderen  */
      y0[j] = yj + delta;
      (*dgl)(*x0, y0, f1);
      for (i = 0; i < n; i++)
        dfdy[i][j] = (f1[i] - f0[i]) / delta;
      y0[j] = yj;
    }

#else                              /* compute appoximations of the    */
    for (j = 0; j < n; j++)        /* partial derivative of the r.h.  */
    {                              /* side wrt. y at (x0,y0) via      */
      yj    = y0[j];               /* central differences and one step*/
      y0[j] = yj - delta;          /* of a Romberg scheme.            */
      (*dgl)(*x0, y0, f0);
      y0[j] = yj + delta;
      (*dgl)(*x0, y0, f1);
      for (i = 0; i < n; i++)
        dfdy[i][j] = (f1[i] - f0[i]) /
                     (TWO * delta);
      y0[j] = yj - HALF * delta;
      (*dgl)(*x0, y0, f0);
      y0[j] = yj + HALF * delta;
      (*dgl)(*x0, y0, f1);
      for (i = 0; i < n; i++)
        dfdy[i][j] = (FOUR * (f1[i] - f0[i]) /
                      delta - dfdy[i][j]) / THREE;
      y0[j] = yj;
    }

#if 1
    (*dgl)(*x0 - delta, y0, f0);   /* as above: partial derivative    */
    (*dgl)(*x0 + delta, y0, f1);   /* wrt. x, this time.              */
    for (i = 0; i < n; i++)
      dfdx[i] = (f1[i] - f0[i])
                / (TWO * delta);
    (*dgl)(*x0 - HALF * delta, y0, f0);
    (*dgl)(*x0 + HALF * delta, y0, f1);
    for (i = 0; i < n; i++)
      dfdx[i] = (FOUR * (f1[i] - f0[i]) /
                 delta - dfdx[i]) / THREE;
#else
    (*dgl)(*x0 + delta, y0, dfdx); /* ditto wrt x, now using forward  */
    for (i = 0; i < n; i++)        /* difference quotient and one     */
      dfdx[i] = (dfdx[i] - f0[i])  /* step of Romberg.                */
                / delta;
#endif
#endif
#ifdef DEBUG
    zeig(dfdy[0], n);
    if (n > 1)
      zeig(dfdy[1], n);
#endif
    for (i = 0; i < n; i++)
    {
      for (sum = ZERO, j = 0; j < n; j++)
        sum += f0[j] * dfdy[i][j];
      dfdx[i] += sum;
    }
#ifdef DEBUG
    zeig(dfdx,    n);
#endif

#ifdef OHNE_ROMBERG
    *aufrufe += n + 2;            /* count calls of dgl()             */
#else
    *aufrufe += 4 * n + 4;
#endif

    null = matrix_null(n, dfdy) ||  /* check whether all partials wrt.*/
           vektor_null(n, dfdx);    /* y or all partials wrt. x       */
                                    /* vanish                         */

    /*******************************************************************
    *           Find  A M O U N T  of  W O R K                         *
    *******************************************************************/

    /* Here we estimate the mean iteration error, find the step size  */
    /* and determine the amount of work.                              */

    /* The mean estimated iteration error is determined by            */
    /*     h^(l+2)/(l+2)! * sqrt(g*(df/dx * (df/dy)^l)^2 / summeg).   */
    /* The square root is computed via  wurzel() and stored in        */
    /* deltak[l]. dk is given as   dk = h^(l+2) * deltak[l] / (l+2)!. */

    /* The step size for using the method of order m is  heps[m-1] =  */
    /*     *eps_rel * (2*m)! /                                        */
    /*     sqrt(g*(df/dx * (df/dy)^(2*m-2))^2 / summeg)^(1/(2*m)).    */

    /* The amount of work is given by                                 */
    /*     AW(eps,m) = (n+1+4*m*m) / heps[m-1].                       */
    /* We choose the order for which the amount of work relative to   */
    /* the number of function evaluations is minimal. This happens    */
    /* for the first m with  AW(eps,m) < AW(eps,m+1).                 */

    if (! null)                  /* neither partial is zero?          */
    {                            /* determine order and step size and */
      m  = 1;                    /* estimate mean iteration error     */
      dk = ZERO;
      for (i = 0; i < n; i++)
        dk += g[i] * sqr(dfdx[i]);
      deltak[0] = SQRT(dk / summeg);
      heps[0]   = SQRT(*eps_rel * fak[0] / deltak[0]);
      aw2       = (n + 5) / (float)heps[0];
      k         = 2;

      do
      {
        aw1             = aw2;
        km2             = 2 * k;
        deltak[km2 - 3] = wurzel(n, dfdy, dfdx, summeg, g, f1);
        deltak[km2 - 2] = wurzel(n, dfdy, dfdx, summeg, g, f1);
#ifdef DEBUG
        printf("deltak[%2d] = %30.19"LZP"e\n", km2-3, deltak[km2-3]);
        printf("deltak[%2d] = %30.19"LZP"e\n", km2-2, deltak[km2-2]);
#endif

        if (deltak[km2 - 2] != ZERO)
        {
          heps[k - 1] = POW(*eps_rel * fak[k - 1] / deltak[km2 - 2],
                            ONE / (TWO * k));
          aw2         = (n + 1 + 4 * k * k) / (float)heps[k - 1];
          k++;
        }
        else
          null = TRUE;
      }
      while (k <= mmax && aw2 < aw1 && ! null);

      m = k - 2;
    }
#ifdef DEBUGg
    zeig(deltak, 2 * mmax);
#endif

    if (null)         /* We would have to divide by zero in computing */
    {                 /* heps[m-1]. Hence we choose the order to be   */
                      /* three and set                                */
      m    = 3;       /* deltak[l] = eps_rel * 10^(5-l) for l=0,..,5  */
      zpot = ONE;
      for (i = 1; i <= 6; i++)
      {
        deltak[6 - i] =  *eps_rel * zpot;
        zpot          *= TEN;
      }
      h = vz * (REAL)0.1;        /* set the starting step size to 0.1 */
    }
    else
      h = vz * sf * heps[m - 1];/* determine step size to be used now */

    z = 0;                             /* initialize for another step */
    x = *x0 + h;
    if (vz * (xend - x) < ZERO)
      h = xend - *x0,
      x = xend;
    ord_up = FALSE;
    if (m != malt)                                      /* new order? */
    {
      mm2 = m * 2;

      /* direct file pointer to the position where the coefficients   */
      /* of order m are in the node file. To do this, we must skip    */
      /* i(i+2) REAL numbers for i = 1,2,...,m-1. Hence the call of   */
      /* fseek().                                                     */

      fseek(stuetzdat,
            (long)((m * (m - 1) * (2 * m + 5) / 6) * FLTSIZE),
            SEEK_SET);

      fread(alpha, FLTSIZE, m, stuetzdat);
      for (k = 0; k < m; k++)
        fread(beta[k], FLTSIZE, m, stuetzdat);
      fread(A, FLTSIZE, m, stuetzdat);
      fread(alphaq, FLTSIZE, m + 1, stuetzdat);
      for (k = 0; k <= m; k++)
        fread(betaq[k], FLTSIZE, m + 1, stuetzdat);
      fread(Aq, FLTSIZE, m + 1, stuetzdat);
    }

    /*******************************************************************
    *     perform  O N E  S T E P                                      *
    *******************************************************************/

    do               /* perform one step for differing h until the    */
    {                /* desired accuracy has been reached.            */
                     /* EXCEPTION: If the step size h is chosen below */
                     /* epserr, this indicates that the method will   */
                     /* not converge.                                 */
      if (ord_up)                               /* change of order?   */
      {                                         /* raise order by one */
        m++;
        mm2 = m * 2;                           /* set up often used   */
                                               /* constant here       */

        for (i = 0; i < m; i++)       /* copy alphaq, betaq and Aq    */
        {                             /* to  alpha,  beta  and A      */
          alpha[i] = alphaq[i];
          A[i]     = Aq[i];
          for (j = 0; j < m; j++)
            beta[i][j] = betaq[i][j];
        }

        fread(alphaq, FLTSIZE, m + 1, stuetzdat);     /* read new     */
        for (k = 0; k <= m; k++)                      /* coefficients */
          fread(betaq[k], FLTSIZE, m + 1, stuetzdat); /* from         */
        fread(Aq, FLTSIZE, m + 1, stuetzdat);         /* file         */

        if (! null)
        {                      /* determine deltak and heps for new m */

          deltak[mm2 - 1] = wurzel(n, dfdy, dfdx, summeg, g, f1);
          deltak[mm2]     = wurzel(n, dfdy, dfdx, summeg, g, f1);

          if (deltak[mm2 - 2] != ZERO)
          {
            heps[m - 1] = POW(*eps_rel * fak[m - 1] /
                              deltak[mm2 - 2], ONE / mm2);
            h           = vz * sf * heps[m - 1];
          }
          else
            null = TRUE;
        }

        if (null)            /* We would have to divide by zero to    */
        {                    /* find heps[m-1] . Hence we set         */
                             /* deltak[l] = eps_rel * 10^(mm2-1-l)    */
          zpot = ONE;        /* for  l=1,...,2*m-1 .                  */
          for (i = 1; i <= mm2; i++)
          {
            deltak[mm2 - i] =  *eps_rel * zpot;
            zpot            *= TEN;
          }
          h = vz * (REAL)0.1;
        }
      }

      /*****************************************************************
      *                  I t e r a t i o n                             *
      *****************************************************************/

      for (i = 0; i < n; i++)              /* initialize weights ki   */
      {                                    /* and kiq to h*f0         */
        hmf0 = h * f0[i];
        for (j = 0; j < m; j++)
          ki[i][j] = kiq[i][j] = hmf0;
        kiq[i][m] = hmf0;
        yalt[i]   = y0[i] + hmf0;          /* store old approximation */
      }                                    /* for error approximation */

      /* We go through the following loop until the desired accuracy  */
      /* has been reached or the number of iterations exceeds 2*m+1.  */
      /* If the latter happens, we raise the order and repeat the     */
      /* step. If the iteration does not converge or the order can no */
      /* longer be raised, we decrease the step size.                 */

      l      = 0;
      faklp1 = ONE;

      do                                 /* start lth iteration       */
      {
        if (*aufrufe > fmax)             /* excessive function calls? */
          stopp  = TRUE,                 /* report reason of stop     */
          fehler = 6;

        l++;
        faklp1 *= l + 1;
        cause  =  0;

        /***************************************************************
        * lth Iteration for the weights ki                             *
        ***************************************************************/

        for (k = 0; k < m; k++)
        {
          x1 = *x0 + h * alpha[k];      /* compute intermediate nodes */
          for (i = 0; i < n; i++)       /* in the interval [x0, x0+h] */
          {
            y[i] = y0[i];
            for (j = 0; j < m; j++)
              y[i] += ki[i][j] * beta[j][k];
          }

          (*dgl)(x1, y, f1);            /* evaluate r.h. side at      */
                                        /* intermediate x-values      */

          for (i = 0; i < n; i++)       /* perform step               */
            db[i][k] = h * f1[i];
        }

        for (i = 0; i < n; i++)         /* store new weights          */
          for (j = 0; j < m; j++)
            ki[i][j] = db[i][j];

        /***************************************************************
        * lth Iteration for the weights  kiq                           *
        ***************************************************************/

        for (k = 0; k <= m; k++)
        {
          x1 = *x0 + h * alphaq[k];     /* intermediate nodes         */
          for (i = 0; i < n; i++)
          {
            yq[i] = y0[i];
            for (j = 0; j <= m; j++)
              yq[i] += kiq[i][j] * betaq[j][k];
          }
          (*dgl)(x1, yq, f1);           /* evaluate r.h. side at      */
                                        /* intermediate nodes         */
          for (i = 0; i < n; i++)       /* perform step               */
            db[i][k] = h * f1[i];
        }
        for (i = 0; i < n; i++)         /* store new weights kiq      */
          for (j = 0; j <= m; j++)
            kiq[i][j] = db[i][j];

        *aufrufe += mm2 + 1;        /* count number of calls of dgl() */

        /***************************************************************
        *                  Approximations of the  lth iteration        *
        ***************************************************************/

        for (i = 0; i < n; i++)      /* compute new approximations y  */
        {                            /* (of order m) and yq (of order */
          yi = yqi = y0[i];          /* m+1) in the lth iteration     */
          for (j = 0; j < m; j++)    /* according to Runge-Kutta      */
            yi  += A[j]  * ki[i][j],
            yqi += Aq[j] * kiq[i][j];
          y[i]  = yi;
          yq[i] = yqi + Aq[m] * kiq[i][m];
        }

        /* Compute the mean absolute and relative difference          */
        /* delta_g and dg_rel of two successive iterations as well as */
        /* the absolute and relative mean difference e and e_rel      */
        /* of the two approximations with different orders.           */

        for (delta_g = e = yqnorm = ZERO, i = 0; i < n; i++)
        {
          yqnorm  += sqr(yq[i]);
          y_diff  =  yq[i] - yalt[i];
          delta_g += g[i] * y_diff * y_diff;
          y_diff  =  y[i] - yq[i];
          e       += g[i] * y_diff * y_diff;
          yalt[i] =  yq[i];
        }
        delta_g = SQRT(delta_g / summeg);
        e       = SQRT(e       / summeg);
        if (yqnorm > ZERO)
          yqnorm = SQRT(yqnorm),
          dg_rel = delta_g / yqnorm,
          e_rel  = e       / yqnorm;
        else
          dg_rel = delta_g,
          e_rel  = e;
        genauigkeit = max(e_rel, dg_rel);

        /***************************************************************
        *         check break-off criteria                             *
        ***************************************************************/

        lh = l;
        if (e_rel >= *eps_rel)       /* The two different order       */
        {                            /* approximations differ by the  */
          cause = 1;                 /* same order as the approximate */
          if (z != 0)                /* solution. Otherwise the step  */
            sf *= (REAL)0.9;         /* size is reduced               */
          z++;
          hbeg = h;

          h *= sf * POW((HALF * *eps_rel / e), ONE / (mm2 + 1));
          l =  0;
        }
        else
        {
          if (null)                         /* compute mean estimated */
            dk = deltak[l - 1];             /* iteration error        */
          else
            dk = FIVE *
                 POW(h, (REAL)(l + 1)) *
                 deltak[l - 1] / faklp1;

          if (delta_g >= dk)        /* As only the local truncation   */
          {                         /* error is used for computing h, */
            cause =  2;             /* we must now check the          */
            hbeg  =  h;             /* convergence of the iterations. */
            h     *= (REAL)0.6;
            sf    *= (REAL)0.8;
            l     =  0;
          }
          else if (dg_rel >= *eps_rel &&  /* If the maximal number of */
                   m      >= mmax - 1 &&  /* iterations has been      */
                   l      >= mm2 - 1      /* reached without the      */
                  )                       /* desired accuracy, we     */
          {                               /* ought to increase the    */
            cause =  3;                   /* order. Since this is     */
            hbeg  =  h;                   /* impossible here, we can  */
            h     *= (REAL)0.8;           /* only decrease the step   */
            l     =  0;                   /* size.                    */
          }
        }

        if (l == 0)  /* l = 0 indicates that the step is unsuccessful */
        {

          PROTOKOLLIEREN(hbeg, lh);

          if (h < epserr)               /* Method does not converge?  */
            stopp  = TRUE,
            fehler = 1;
          else                          /* Method does converge?      */
          {                             /* erase last step and repeat */
            faklp1 = ONE;
            x      = *x0 + h;
            if (vz * (xend - x) < ZERO)
              h = xend - *x0,
              x = xend;
            for (i = 0; i < n; i++)
            {
              hmf0 = h * f0[i];
              for (j = 0; j < m; j++)
                ki[i][j] = kiq[i][j] = hmf0;
              kiq[i][m] = hmf0;
              yalt[i]   = y0[i] + hmf0;
            }
          }
        }
      }                                   /* end of the lth iteration */
      while (((l < mm2 - 1 && dg_rel >= *eps_rel) ||
              l == 0) &&
             ! stopp
            );

      if (dg_rel >= *eps_rel &&    /* acuracy not reached after the   */
          m + 1  <  mmax     &&    /* theoretical number of steps:    */
          ! stopp                  /* Increase order and repeat ...   */
         )
      {
        ord_up = TRUE;
        cause  = 4;
        PROTOKOLLIEREN(h, l);
      }
      else
        ord_up = FALSE;
    }                               /* until desired accuracy reached */
    while (ord_up);

    if (! stopp)                    /* Step successful?               */
    {
      PROTOKOLLIEREN(h, l);

      if (ausdat != NULL)
      {
        fprintf(ausdat, "\n%5d %22.15"LZP"e    1"
                        "    %23.15"LZP"e  %16.9"LZP"e",
                        sz, x, yq[0], genauigkeit);
        for (i = 1; i < n; i++)
          fprintf(ausdat, "\n%33d    %23.15"LZP"e", i, yq[i]);
      }
    }

    epslok = max(epslok, genauigkeit);
    if (x != xend && ! stopp)                  /* did not reach xend? */
    {
      if (z > 1)
        sf /= (REAL)0.97;
      *x0 = x;
      for (i = 0; i < n; i++)
        y0[i] = yq[i];
    }

    else                  /* stop subroutine since the IVP has been   */
    {                     /* solved (fehler = 0), or since the method */
      if (! stopp)        /* will not converge (fehler = 1), or since */
        *x0 = x;          /* the maximal number of calls of dgl()     */
                          /* (fmax) has been exceeded (fehler = 6)    */
      else                        /* xend not reached? => copy result */
        for (i = 0; i < n; i++)   /* of last successful step to yq    */
          yq[i] = y0[i];

      stopp    = TRUE;
      *eps_rel = epslok;
    }
  }
  while (! stopp);                  /* finished with the integration? */

  fclose(stuetzdat);

  if (*file_out != '\0')           /* Did we form an output file?     */
    fclose(ausdat);
  if (*file_pro != '\0')           /* Did we form a protocol file?    */
    fclose(prodat);

  vmfree(vmblock);
  return fehler;
}

/* -------------------------- END implruku.c ------------------------ */
