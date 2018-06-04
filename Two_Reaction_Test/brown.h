/*.BA*/



/*.FE{C 6.2.4}{Brown's Method for Nonlinear Systems}
              {Brown's Method for Nonlinear Systems}*/

/*.BE*/
/* ---------------------- DECLARATIONS brown.h ---------------------- */

#ifndef BROWN_H_INCLUDED
#define BROWN_H_INCLUDED

/* Type of function, which evaluates the right hand side f of the kth */
/* equation of a (normally nonlinear) system of equations at x        */
typedef int (*nlgls)(int k, REAL x[], REAL *f);
/*.IX{nlgls}*/

int brown   /* Brown's method for nonlinear systems of equations .....*/
         (
          nlgls   fkt,      /* Function ..............................*/
          int     n,        /* number of equations ...................*/
          REAL    x0[],     /* Starting value for iteration ..........*/
          REAL    eps,      /* error bound ...........................*/
          int     prot,     /* Protokol switch .......................*/
          int     maxit,    /* maximal number of steps ...............*/
          REAL    x1[],     /* solution ..............................*/
          int     *itanz    /* actual steps performed ................*/
         );                 /* error code ............................*/

#endif

/* -------------------------- END brown.h --------------------------- */
