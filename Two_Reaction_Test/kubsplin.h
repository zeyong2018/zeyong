/*.BA*/
/*.KA{C 10}{Interpolating Polynomial Splines}
           {Interpolating Polynomial Splines for Constructing Smooth
            Curves}*/
/*.BE*/
/*.FE{C 10.1}
     {Cubic Polynomial Splines}
     {Cubic Polynomial Splines}*/

/* -------------------- DECLARATIONS kubsplin.h --------------------- */

int spline      /* non-parametric cubic splines ......................*/
          (
           int  m,            /* number of nodes .....................*/
           REAL x[],          /* nodes : x-values ....................*/
           REAL y[],          /*         y-values ....................*/
           int  marg_cond,    /* type of end point condition .........*/
           REAL marg_0,       /* left end point condition ............*/
           REAL marg_n,       /* right end point condition ...........*/
           int  save,         /* save aux vectors ? ..................*/
           REAL b[],          /* Spline coefficients of (x-x[i]) .....*/
           REAL c[],          /* Spline coefficients of (x-x[i])^2 ...*/
           REAL d[]           /* Spline coefficients of (x-x[i])^3 ...*/
          );                  /* error code ..........................*/

int parspl           /* parametric cubic splines .....................*/
          (
           int  m,            /* number of nodes .....................*/
           REAL x[],          /* nodes : x-values ....................*/
           REAL y[],          /*         y-values ....................*/
           int  marg_cond,    /* type of end point condition .........*/
           REAL marg_0[],     /* left end point condition ............*/
           REAL marg_n[],     /* right end point condition ...........*/
           int  cond_t,       /* Parameter nodes given ? .............*/
           REAL t[],          /* Parameter nodes .....................*/
           REAL bx[],         /* x spline coeffic. for (t-t[i]) ......*/
           REAL cx[],         /* x spline coeffic. for (t-t[i])^2 ....*/
           REAL dx[],         /* x spline coeffic. for (t-t[i])^3 ....*/
           REAL by[],         /* y spline coeffic. for (t-t[i]) ......*/
           REAL cy[],         /* y spline coeffic. for (t-t[i])^2 ....*/
           REAL dy[]          /* y spline coeffic. for (t-t[i])^3 ....*/
          );                  /* error code ..........................*/

int spltrans     /* transformed parametric cubic splines .............*/
            (
             int  m,         /* number of nodes ......................*/
             REAL x[],       /* nodes : x-values .....................*/
             REAL y[],       /*         y-values .....................*/
             int  mv,        /* type of transformation of origin .....*/
             REAL px[],      /* coordinats of the transformation      */
             REAL py[],      /* vector  P ............................*/
             REAL a[],       /* Spline coeff. of (phi-phin[i])^0 .....*/
             REAL b[],       /* Spline coeff. of (phi-phin[i]) .......*/
             REAL c[],       /* Spline coeff. of (phi-phin[i])^2 .....*/
             REAL d[],       /* Spline coeff. of (phi-phin[i])^3 .....*/
             REAL phin[],    /* angular coordinates of nodes .........*/
             REAL *phid      /* angle of rotation of coordinate system*/
            );               /* error code ...........................*/

/* ------------------------- END kubsplin.h ------------------------- */
