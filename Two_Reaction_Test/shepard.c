#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE shepard.c ------------------------ */

/***********************************************************************
*                                                                      *
* Shepard interpolation of surfaces in  R3                             *
* ----------------------------------------                             *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Bjoern Terwege, 6.12.1995 (FORTRAN 77)         *
* Adapted by :          Juergen Dietel, Computing Center, RWTH Aachen  *
* Source:               FORTRAN 77 source code                         *
* Date:                 2.21.1996                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>     /*  for  REAL, POSMIN, ZERO, TWO, SQRT, sqr,   */
                       /*       POW, ONE                              */
#include <vmblock.h>   /* vminit, vmalloc, VEKTOR, vmcomplete, vmfree */
#include <shepard.h>   /*  for  shepard                               */


/* put out a REAL vector of length `n' (and its name)                 */
#define zeig(v, n)                   \
  {                                  \
    int i;                           \
    printf("%-9s", #v": ");          \
    for (i = 0; i < n; i++)          \
      printf("%11.4"LZP"e", v[i]);   \
    printf("\n");                    \
  }



/*--------------------------------------------------------------------*/
/*.BA*/

int shepard       /* Shepard interpolation (global,local, F-L weights)*/
/*.IX{shepard}*/
    (
     REAL x0,           /* (x0,y0) = Interpolation point  ............*/
     REAL y0,
     REAL x[],          /* (x[i],y[i]) = nodes         ...............*/
     REAL y[],
     REAL f[],          /* f[i] = function values at nodes ...........*/
     int  n,            /* number of nodes  - 1         ..............*/
     REAL mue,          /* adjustable Shepard parameter (>0) .........*/
     int  methode,      /* global, local, F.-Little weights: (0,1,2)? */
     REAL R,            /* Radius for circle around (x0,y0) ..........*/
     REAL *PHI          /* Interpolation value at (x0,y0) ............*/
    )                   /* error code  ...............................*/

/***********************************************************************
*   This function computes one functional value at (X,Y) for given     *
*   nodes using the  Shepard method. The user can select betwen the    *
*   global, the local and the local method with Franke-Little weights. *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0  \    (x0,y0) point for desired Shepard interpolation             *
* y0  /            evaluation                                          *
* x        [0..n] vector with x-coordinates of nodes                   *
* y        [0..n] vector with y-coordinates of nodes                   *
* f        [0..n] vector with corresponding f values                   *
* n        number of nodes - 1   (n = index for last node)             *
* mue      freely choosable parameter for the Shepard method, used for *
*          setting the exponents for the weights  (mue > 0);           *
*          for goood results we recommend to set  2 < mue < 6.         *
*          If mue <= 0, we set mue = 2 internally.                     *
* methode  Number for the variant of Shepard method to be used         *
*          = 0: global method                                          *
*          = 1: local method                                           *
*          = 2: local method with Franke-Little weights                *
* R        Radius for the local method; determines the circle around   *
*          (x0,y0), in which nodes are used for the interpolation;     *
*          nodes outside this circle are ignored.                      *
*          R should be chosen large enough so that the circle contains *
*          sufficiently many nodes.                                    *
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* PHI      value of the interpolation at (x0,y0)                       *
*                                                                      *
* Funktion value:                                                      *
* ==============                                                       *
* Error code:                                                          *
* = 0: all is ok                                                       *
* = 1: invalid input parameter:                                        *
*      n < 0  or  methode != 0,1,2  or  R <= 0                         *
* = 2: All weights w[i] are zero.                                      *
* = 3: lack of disk space                                              *
*                                                                      *
* global variables:                                                    *
* =================                                                    *
* REAL, POSMIN, ZERO, vminit, vmalloc, VEKTOR, vmcomplete, vmfree,     *
* TWO, SQRT, sqr, POW, ONE                                             *
.BA*)
***********************************************************************/
/*.BE*/

{
  REAL *r;          /* [0..n] vector with euclidean distances of the  */
                    /* nodes to (x0,y0)                               */
  REAL *w;          /* [0..n] vector of weights, depends on r and mue */
  REAL *psi = NULL; /* [0..n] aux vector for weights in local method  */
  REAL *xi = NULL;  /* [0..n] aux vector for Franke-Little weights    */
  void *vmblock;    /* list of dynamically allocated vectors          */
  REAL norm;        /* 1-norm of weight vector before normalisation   */
  int  j;           /* Loop variable                                  */


  *PHI = POSMIN;                   /* Output parameter for faulty run */
                                   /* preset to a nonsense value      */

  if (n < 0 ||                      /* invalid value for n?           */
      methode < 0 || methode > 2)   /* invalid method code?           */
    return 1;                       /* return error                   */

  if (methode != 0)                 /* not the global method?         */
    if (R <= ZERO)                  /* improper value for R?          */
      if (methode == 2)             /* local method with Franke-Little*/
        R = (REAL)0.1;              /* correct R                      */
      else                          /* different method?              */
        return 1;                   /* return error                   */

  vmblock = vminit();               /* initialize storage             */
  r   = (REAL *)vmalloc(vmblock,    /* storage for three vectors      */
        VEKTOR, n + 1, 0);
  w   = (REAL *)vmalloc(vmblock,
        VEKTOR, n + 1, 0);
  if (methode != 0)
    psi = (REAL *)vmalloc(vmblock,
          VEKTOR, n + 1, 0),
    xi  = psi;                      /* another name for storage       */
  if (! vmcomplete(vmblock))        /* lack of memory?                */
  {
    vmfree(vmblock);           /* free dynamically allocated storage  */
    return 3;                       /* return error                   */
  }

  if (mue <= ZERO)                  /* improper value for mue?        */
    mue = TWO;                      /* use default value = 2          */

  for (j = 0; j <= n; j++)          /* compute distances  r[j]        */
  {
    r[j] = SQRT(sqr(x0 - x[j]) +
                sqr(y0 - y[j]));
    if (r[j] == ZERO)               /* one distance is zero, i.e.,    */
    {                               /* (x0,y0) is a node              */
      *PHI = f[j];                  /* return this node and its       */
      vmfree(vmblock);              /* function value                 */
      return 0;
    }
  }


  switch (methode)         /* set up weight vector, find 1-norm,      */
  {                        /* global method ?                         */
    case 0:
      for (j = 0, norm = ZERO; j <= n; j++)
        w[j] = ONE / POW(r[j], mue),
        norm += w[j];
      break;

    case 1:                                    /* local method?       */
      for (j = 0; j <= n; j++)                 /* compute  psi[j]     */
        if (r[j] >= R)
          psi[j] = ZERO;
        else
          psi[j] = (R / r[j]) - ONE;

      for (j = 0, norm = ZERO; j <= n; j++)    /* compute unnormalise */
        if (psi[j] != ZERO)                    /* weights from psi,   */
          w[j] = ONE / POW(psi[j], mue),
          norm += w[j];                        /* and their sum       */
        else
          w[j] = ZERO;
      break;

    case 2:                                    /* local methodD with  */
                                            /* Franke-Little weights? */
      for (j = 0; j <= n; j++)                 /* compute xi[j]       */
        if (r[j] >= R)
          xi[j] = ZERO;
        else
          xi[j] = ONE - r[j] / R;

      for (j = 0, norm = ZERO; j <= n; j++)    /* compute unnormalised*/
        w[j] =  POW(xi[j], mue),               /* weights from xi     */
        norm  += w[j];                         /* and their sum       */
      break;

    default:                                   /* impossible !        */
      vmfree(vmblock);
      return 1;
  }

  if (norm == ZERO)                            /* All weights  w[j]   */
  {                                            /* are zero?           */
    vmfree(vmblock);
    return 2;                                  /* return error        */
  }
  for (j = 0; j <= n; j++)                     /* normalize weights   */
    if (w[j] != ZERO)
      w[j] /= norm;
#ifdef DEBUG
  zeig(w, n + 1);
  printf("  norm=%11.4"LZP"e\n", norm);
#endif


  *PHI = ZERO;                     /* compute value of function at    */
  for (j = 0; j <= n; j++)         /* (x0,y0)                         */
    *PHI += w[j] * f[j];


  vmfree(vmblock);
  return 0;
}

/* ------------------------- END  shepard.c ------------------------- */
