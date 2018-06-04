/*.BA*/
/*.KA{C 13}{Akima and Renner Subsplines}
           {Akima and Renner Subsplines}*/
/*.BE*/
/* -------------------- DECLARATIONS subsplin.h --------------------- */

int akima                 /* compute coeffic. of an Akima subspline ..*/
        (
         int  *n,                  /* Number of final node ...........*/
         int  nmax,                /* upper index limit for nodes x, y*/
         REAL x[],                 /* nodes: x-values ................*/
         REAL y[],                 /*        y-values ............... */
         int  perio,               /* periodic interpolation? ........*/
         REAL beta,                /* rounding parameter .............*/
         REAL b[],                 /* Spline coefficients ............*/
         REAL c[],
         REAL d[]
        );                         /* error code                      */

int renner                /* compute renner subspline coefficients ...*/
        (
         int  *n,                  /* number of last node ............*/
         int  nmax,                /* upper index limit for P ........*/
         int  dim,                 /* space dimension (2 or 3) .......*/
         REAL *P[],                /* nodes ..........................*/
         REAL beta,                /* rounding parameter .............*/
         REAL T[],                 /* lengths of parameter intervals  */
         REAL *b[],                /* coeffizients for t^1 ...........*/
         REAL *c[],                /* coeffizients for t^2 ...........*/
         REAL *d[]                 /* Koeffizients for t^3 ...........*/
        );                         /* error code .....................*/

void rennwert             /* evaluation of a Renner subspline ........*/
        (
         int  n,                   /* number of spline pieces ........*/
         int  dim,                 /* space dimension (2 or 3) .......*/
         REAL twert,               /* point of evaluation ............*/
         REAL t[],                 /* nodes ..........................*/
         REAL *a[],                /* spline coeff. of (t-t[i])^0 ....*/
         REAL *b[],                /* spline coeff. of (t-t[i])^1 ....*/
         REAL *c[],                /* spline coeff. of (t-t[i])^2 ....*/
         REAL *d[],                /* spline coeff. of (t-t[i])^3 ....*/
         REAL s[],                 /* point on the spline curve ......*/
         REAL ds[4][3]             /* 0th - 3rd derivative of spline  */
        );

int renntab               /* table of values of a Renner subspline ...*/
        (
         int  n,                   /* number of spline pieces ........*/
         REAL tanf,                /* left end point of interval .....*/
         REAL tend,                /* right end point of interval ....*/
         REAL delt,                /* step size ......................*/
         int  anzahl,              /* maximal size of table ..........*/
         REAL t[],                 /* nodes ..........................*/
         int  dim,                 /* space dimension (2 or 3) .......*/
         REAL *a[],                /* spline coeffic. of (t-t[i])^0 ..*/
         REAL *b[],                /* spline coeffic. of (t-t[i])^1 ..*/
         REAL *c[],                /* spline coeffic. of (t-t[i])^2 ..*/
         REAL *d[],                /* spline coeffic. of (t-t[i])^3 ..*/
         REAL *sptab[],            /* table of spline values .........*/
         int  *lentab              /* actual table length ............*/
        );                         /* error code .....................*/

/* ------------------------- END subsplin.h ------------------------- */
