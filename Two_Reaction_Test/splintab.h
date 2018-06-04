/*.BA*/



/*.FE{}{Tabulating of Polynomial Splines}
       {Tabulating of Polynomial Splines}*/

/*.BE*/
/* -------------------- DECLARATIONS splintab.h --------------------- */

int sptab       /* Table of values of a cubic spline .................*/
         (
          int  n,         /* number of spline pieces ( = # nodes - 1) */
          REAL xanf,      /* left end point of interval ..............*/
          REAL xend,      /* right end point of interval .............*/
          REAL deltx,     /* step size ...............................*/
          int  anzahl,    /* maximal size of table ...................*/
          REAL x[],       /* nodes ...................................*/
          REAL a[],       /* Spline coefficients for (x-x[i])^0 ......*/
          REAL b[],       /* Spline coefficients for (x-x[i])^1 ......*/
          REAL c[],       /* Spline coefficients for (x-x[i])^2 ......*/
          REAL d[],       /* Spline coefficients for (x-x[i])^3 ......*/
          REAL xtab[],    /* x-coordinates of the table ..............*/
          REAL ytab[],    /* y-coordinates of the table ..............*/
          int  *lentab    /* actual table length .....................*/
         );               /* error code ..............................*/

int partab    /* Table of values for a parametric cubic spline .......*/
          (
           int  n,         /* number of spline pieces ................*/
           REAL tanf,      /* left end point of interval .............*/
           REAL tend,      /* right end point of interval ............*/
           REAL delt,      /* step size ..............................*/
           int  anzahl,    /* maximal size of table ..................*/
           REAL t[],       /* parameter nodes ........................*/
           REAL ax[],      /* x spline coefficients for (t-t[i])^0 ...*/
           REAL bx[],      /* x spline coefficients for (t-t[i])^1 ...*/
           REAL cx[],      /* x spline coefficients for (t-t[i])^2 ...*/
           REAL dx[],      /* x spline coefficients for (t-t[i])^3 ...*/
           REAL ay[],      /* y spline coefficients for (t-t[i])^0 ...*/
           REAL by[],      /* y spline coefficients for (t-t[i])^1 ...*/
           REAL cy[],      /* y spline coefficients for (t-t[i])^2 ...*/
           REAL dy[],      /* y spline coefficients for (t-t[i])^3 ...*/
           REAL xtab[],    /* x-coordinates of table .................*/
           REAL ytab[],    /* y-coordinates of table .................*/
           int  *lentab    /* actual size of table ...................*/
          );               /* error code .............................*/

int hmtab         /* Table of values for a Hermite spline ............*/
         (
          int  n,         /* number of spline pieces .................*/
          REAL xanf,      /* left end point for tabulating interval ..*/
          REAL xend,      /* right end point .........................*/
          REAL deltx,     /* step size ...............................*/
          int  anzahl,    /* maximal length of table .................*/
          REAL x[],       /* nodes ...................................*/
          REAL a[],       /* Spline coefficients for (x-x[i])^0 ......*/
          REAL b[],       /* Spline coefficients for (x-x[i])^1 ......*/
          REAL c[],       /* Spline coefficients for (x-x[i])^2 ......*/
          REAL d[],       /* Spline coefficients for (x-x[i])^3 ......*/
          REAL e[],       /* Spline coefficients for (x-x[i])^4 ......*/
          REAL f[],       /* Spline coefficients for (x-x[i])^5 ......*/
          REAL xtab[],    /* x-coordinates of table ..................*/
          REAL ytab[],    /* y-coordinates of table ..................*/
          int  *lentab    /* actual size of table ....................*/
         );               /* error code ..............................*/

int pmtab   /* Table of values for a parametric Hermite spline .......*/
         (
          int  n,         /* number of spline pieces .................*/
          REAL tanf,      /* left end point of interval ..............*/
          REAL tend,      /* right end point .........................*/
          REAL delt,      /* step size ...............................*/
          int  anzahl,    /* maximal size of table ...................*/
          REAL t[],       /* nodes ...................................*/
          REAL ax[],      /* x spline coefficients for (t-t[i])^0 ....*/
          REAL bx[],      /* x spline coefficients for (t-t[i])^1 ....*/
          REAL cx[],      /* x spline coefficients for (t-t[i])^2 ....*/
          REAL dx[],      /* x spline coefficients for (t-t[i])^3 ....*/
          REAL ex[],      /* x spline coefficients for (t-t[i])^4 ....*/
          REAL fx[],      /* x spline coefficients for (t-t[i])^5 ....*/
          REAL ay[],      /* y spline coefficients for (t-t[i])^0 ....*/
          REAL by[],      /* y spline coefficients for (t-t[i])^1 ....*/
          REAL cy[],      /* y spline coefficients for (t-t[i])^2 ....*/
          REAL dy[],      /* y spline coefficients for (t-t[i])^3 ....*/
          REAL ey[],      /* y spline coefficients for (t-t[i])^4 ....*/
          REAL fy[],      /* y spline coefficients for (t-t[i])^5 ....*/
          REAL xtab[],    /* x-coordinates of spline in table ........*/
          REAL ytab[],    /* y-coordinates of spline .................*/
          int  *lentab    /* actual size of table ....................*/
         );               /* error code ..............................*/

int strtab /* Table of values for transformed parametric cubic spline */
          (
           int  n,         /* number of spline pieces ................*/
           REAL panf,      /* starting angle for table ...............*/
           REAL pend,      /* final angle of table ...................*/
           REAL phin[],    /* angular nodes ..........................*/
           REAL a[],       /* Spline coeff. for (phi-phin[i])^0 ......*/
           REAL b[],       /* Spline coeff. for (phi-phin[i])^1 ......*/
           REAL c[],       /* Spline coeff. for (phi-phin[i])^2 ......*/
           REAL d[],       /* Spline coeff. for (phi-phin[i])^3 ......*/
           REAL phid,      /* angle of rotation of coordinates .......*/
           REAL px,        /* x-coordinate,                           */
           REAL py,        /* y-coordinate of translation vector .....*/
           REAL x[],       /* nodes: x-values ........................*/
           REAL y[],       /*        y-values ........................*/
           int  nl,        /* maximal length of table ................*/
           int  *nt,       /* actual length of table .................*/
           REAL xtab[],    /* x-coordinates in table .................*/
           REAL ytab[]     /* y-coordinates in table .................*/
          );               /* error code .............................*/

/* ------------------------- END splintab.h ------------------------- */
