/*.BA*/



/*.FE{}{Evaluation of Polynomial Splines}
       {Evaluation of Polynomial Splines}*/

/*.BE*/
/* -------------------- DECLARATIONS spliwert.h --------------------- */

REAL spwert       /* evaluate a cubic spline .........................*/
           (
            int  n,        /* number of spline pieces ................*/
            REAL xwert,    /* x-value ................................*/
            REAL a[],      /* Spline coefficients of (x-x[i])^0 ......*/
            REAL b[],      /* Spline coefficients of (x-x[i])^1 ......*/
            REAL c[],      /* Spline coefficients of (x-x[i])^2 ......*/
            REAL d[],      /* Spline coefficients of (x-x[i])^3 ......*/
            REAL x[],      /* nodes (x-values) .......................*/
            REAL ausg[]    /* 1st, 2nd and 3rd derivatives of spline .*/
           );              /* Functional value for spline ............*/

void pspwert  /* Evaluate a parametric cubic spline ..................*/
            (
             int      n,        /* number of nodes ...................*/
             REAL     twert,    /* place for evaluation ..............*/
             REAL     t[],      /* nodes .............................*/
             REAL     ax[],     /* x spline coeff. of (t-t[i])^0 .....*/
             REAL     bx[],     /* x spline coeff. of (t-t[i])^1 .....*/
             REAL     cx[],     /* x spline coeff. of (t-t[i])^2 .....*/
             REAL     dx[],     /* x spline coeff. of (t-t[i])^3 .....*/
             REAL     ay[],     /* y spline coeff. of (t-t[i])^0 .....*/
             REAL     by[],     /* y spline coeff. of (t-t[i])^1 .....*/
             REAL     cy[],     /* y spline coeff. of (t-t[i])^2 .....*/
             REAL     dy[],     /* y spline coeff. of (t-t[i])^3 .....*/
             REAL     *sx,      /* x-coordinate, .....................*/
             REAL     *sy,      /* y-coordinate of spline value ......*/
             abl_mat1 ausp      /* 0 to third derivatives of spline ..*/
            );

REAL hmtwert        /* Evaluate a Hermite spline .....................*/
            (
             int  n,        /* number of nodes - 1 ...................*/
             REAL x0,       /* place of evaluation ...................*/
             REAL a[],      /* Spline coefficient of (x-x[i])^0 ......*/
             REAL b[],      /* Spline coefficient of (x-x[i])^1 ......*/
             REAL c[],      /* Spline coefficient of (x-x[i])^2 ......*/
             REAL d[],      /* Spline coefficient of (x-x[i])^3 ......*/
             REAL e[],      /* Spline coefficient of (x-x[i])^4 ......*/
             REAL f[],      /* Spline coefficient of (x-x[i])^5 ......*/
             REAL x[],      /* n+1 nodes .............................*/
             REAL ausg[]    /* 1st to 5th derivatives of spline ......*/
            );              /* Function value of spline ..............*/

void pmtwert  /* Evaluate a  parametric Hermite spline ...............*/
            (
             int      n,        /* number of nodes - 1 ...............*/
             REAL     twert,    /* place of evaluation ...............*/
             REAL     t[],      /* nodes (x-values) ..................*/
             REAL     ax[],     /* x spline coeff. of (t-t[i])^0 .....*/
             REAL     bx[],     /* x spline coeff. of (t-t[i])^1 .....*/
             REAL     cx[],     /* x spline coeff. of (t-t[i])^2 .....*/
             REAL     dx[],     /* x spline coeff. of (t-t[i])^3 .....*/
             REAL     ex[],     /* x spline coeff. of (t-t[i])^4 .....*/
             REAL     fx[],     /* x spline coeff. of (t-t[i])^5 .....*/
             REAL     ay[],     /* y spline coeff. of (t-t[i])^0 .....*/
             REAL     by[],     /* y spline coeff. of (t-t[i])^1 .....*/
             REAL     cy[],     /* y spline coeff. of (t-t[i])^2 .....*/
             REAL     dy[],     /* y spline coeff. of (t-t[i])^3 .....*/
             REAL     ey[],     /* y spline coeff. of (t-t[i])^4 .....*/
             REAL     fy[],     /* y spline coeff. of (t-t[i])^5 .....*/
             REAL     *sx,      /* x-coordinate, .....................*/
             REAL     *sy,      /* y-coordinate of spline value ......*/
             abl_mat2 ausp      /* 0th to fifth derivatives of spline */
            );

int strwert  /* Evaluate a transformed parametric cubic spline .......*/
           (
            REAL phi,        /* place of evaluation  .................*/
            int  n,          /* number of nodes - 1 ..................*/
            REAL phin[],     /* angular nodes ........................*/
            REAL a[],        /* Spline coeff. of (phi-phin[i])^0 .....*/
            REAL b[],        /* Spline coeff. of (phi-phin[i])^1 .....*/
            REAL c[],        /* Spline coeff. of (phi-phin[i])^2 .....*/
            REAL d[],        /* Spline coeff. of (phi-phin[i])^3 .....*/
            REAL phid,       /* angle of plane rotation ..............*/
            REAL px,         /* coordinates of translation vector P ..*/
            REAL py,
            REAL ablei[],    /* 0th to third derivatives wrt. x ......*/
            REAL *xk,        /* x-coordinate,                         */
            REAL *yk,        /* y-coordinate of spline at phi ........*/
            REAL *c1,        /* 1st derivative (dr/dphi) at phi ......*/
            REAL *ckr        /* curvature of spline at phi ...........*/
           );                /* error code ...........................*/

/* ------------------------- END spliwert.h ------------------------- */
