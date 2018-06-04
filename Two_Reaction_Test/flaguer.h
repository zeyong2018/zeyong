/*.BA*/



/*.FE{C 3.3.5}{The Laguerre Method}{The Laguerre Method}*/

/*.BE*/
/* --------------------- DECLARATIONS flaguer.h --------------------- */

int laguerre             /* real polynomial roots, Method of Laguerre */
            (
             int  n,         /* Polynomial degree (>= 3) .............*/
             REAL a[],       /* Polynomial coefficients (ascending) ..*/
             REAL abserr,    /* absolute error bound .................*/
             REAL relerr,    /* relative error bound .................*/
             int  maxit,     /* maximal number of iterations .........*/
             REAL x[],       /* real roots ...........................*/
             int  iter[],    /* Iterations per root ..................*/
             int  *nulanz    /* Number of found roots ................*/
            );               /* Error code ...........................*/

/* ------------------------- END flaguer.h -------------------------- */
