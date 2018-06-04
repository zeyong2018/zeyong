/*.BA*/



/*.FE{C 2.8.5}{Zeroin Method}{Zeroin Method}*/

/*.BE*/
/* --------------------- DECLARATIONS fzeroin.h --------------------- */

int zeroin                       /* Find roots with the Zeroin method */
          (
           REALFCT fkt,         /* Function ..........................*/
           REAL    *abserr,     /* absolute error bound ..............*/
           REAL    *relerr,     /* relative error bound ..............*/
           int     fmax,        /* maximal number of calls for fkt() .*/
           char    *protnam,    /* Name of the log file ..............*/
           REAL    a,           /* [a,b] = inclusion interval ........*/
           REAL    *b,          /* right endpoint or zero ............*/
           REAL    *fb,         /* Function value at the root b ......*/
           int     *fanz        /* number of actual function calls ...*/
          );                    /* error code ........................*/

/* ------------------------- END fzeroin.h -------------------------- */
