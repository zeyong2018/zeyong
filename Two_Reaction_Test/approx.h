/*.BA*/
/*.KA{C 8}{Linear and Nonlinear Approximation}
          {Linear and Nonlinear Approximation}*/
/*.BE*/
/* --------------------- DECLARATIONS approx.h ---------------------- */

int gfq            /* Normal equations for least square approximation */
    (
     int  n,                 /* degree of approximation polynomial ...*/
     int  m,                 /* number of nodes - 1 ..................*/
     REAL x[],               /* nodes: x-values ......................*/
     REAL y[],               /*        y-values ......................*/
     REAL w[],               /* weights ..............................*/
     REAL c[]                /* coeffic. of approximating polynomial  */
    );                       /* error code ...........................*/

int pol_appr       /* discr. lin. least squares via orthog. polynom.  */
    (
     int  n,                 /* degree of least square appr. polynom. */
     int  m,                 /* number of nodes - 1 ..................*/
     REAL x[],               /* nodes: x-values ......................*/
     REAL y[],               /*        y-values ......................*/
     REAL w[],               /* weights ..............................*/
     REAL c[],               /* coefficients of optimal polynomial ...*/
     REAL b[],               /* aux variables for lower degree .......*/
     REAL d[]                /* orthogonal polynomials ...............*/
    );                       /* error code ...........................*/

REAL opolwert      /* Evaluate the polynomial from  pol_appr() .......*/
    (
     int  n,                 /* degree of polynomial .................*/
     REAL x,                 /* x-value ..............................*/
     REAL b[],               /* coeffici. for orthogonal polynomials  */
     REAL d[],
     REAL c[]                /* coefficients of optimal polynomial ...*/
    );                       /* value of polynomial at x .............*/

int opolkoeff      /* orthogonal coeffic.  -->  standard coeffic. ....*/
    (
     int  n,                 /* degree of polynomial .................*/
     REAL b[],               /* coefficients for computation .........*/
     REAL d[],               /* of orthogonal polynomials ............*/
     REAL c[],               /* coeff. for c[0]*Q_0 + ... + c[n]*Q_n  */
     REAL a[]                /* coeff. for a[0]*x^0 + ... + a[n]*x^n  */
    );                       /* error code ...........................*/

int lin_happr      /* lin. least squares via Householder transform. ..*/
    (
     int       m,            /* number of nodes ......................*/
     int       n,            /* number of functions -1 ...............*/
     REAL      x[],          /* nodes : x-values .....................*/
     REAL      y[],          /*         y-values .....................*/
     REAL      w[],          /* positive weights .....................*/
     ansatzfnk phi,          /* model functions ......................*/
     REAL      c[],          /* optimal coefficients .................*/
     REAL      *r            /* error of the opt. solution ...........*/
    );                       /* error code ...........................*/

REAL lin_hwert     /* Evaluate function from lin_happr() .............*/
    (
     REAL      x0,           /* x-value ..............................*/
     int       n,            /* number of functions - 1 ..............*/
     ansatzfnk phi,          /* model functions ......................*/
     REAL      c[]           /* optimal coefficients .................*/
    );                       /* return value .........................*/

/* -------------------------- END approx.h -------------------------- */
