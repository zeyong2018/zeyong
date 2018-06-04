


/*.FE{C 10.2}
     {Hermite Splines of Fifth Degree}
     {Hermite Splines of Fifth Degree}*/

/* -------------------- DECLARATIONS hrmsplin.h --------------------- */

int hermit         /* non-parametric Hermite spline ..................*/
          (
           int  m,            /* number of  nodes ....................*/
           REAL x[],          /* nodes: x-values .....................*/
           REAL y[],          /*        y-values .....................*/
           REAL y1[],         /* first derivative vector at nodes ....*/
           int  marg_cond,    /* type of boundary condition ..........*/
           REAL marg_0,       /* left boundary condition .............*/
           REAL marg_n,       /* right boundary condition ............*/
           int  save,         /* save dynamic aux arrays ? ...........*/
           REAL c[],          /* Spline coefficients of (x-x[i])^2 ...*/
           REAL d[],          /* Spline coefficients of (x-x[i])^3 ...*/
           REAL e[],          /* Spline coefficients of (x-x[i])^4 ...*/
           REAL f[]           /* Spline coefficients of (x-x[i])^5 ...*/
          );                  /* error code ..........................*/

int parmit              /* parametric Hermite splines ................*/
          (
           int  m,           /* number of nodes ......................*/
           REAL x[],         /* nodes : x-values .....................*/
           REAL y[],         /*         y-values .....................*/
           int  richt,       /* type of derivative ...................*/
           REAL xricht[],    /* Tangent or normal direction or only   */
           REAL yricht[],    /* dy/dx in yricht ......................*/
           int  marg,        /* type of end point condition ..........*/
           REAL corn_1[],    /* left hand end point condition ........*/
           REAL corn_n[],    /* right hand end point condition .......*/
           REAL cx[],        /* x spline coeffic. for (t-t[i])^2 .....*/
           REAL dx[],        /* x spline coeffic. for (t-t[i])^3 .....*/
           REAL ex[],        /* x spline coeffic. for (t-t[i])^4 .....*/
           REAL fx[],        /* x spline coeffic. for (t-t[i])^5 .....*/
           REAL cy[],        /* y spline coeffic. for (t-t[i])^2 .....*/
           REAL dy[],        /* y spline coeffic. for (t-t[i])^3 .....*/
           REAL ey[],        /* y spline coeffic. for (t-t[i])^4 .....*/
           REAL fy[],        /* y spline coeffic. for (t-t[i])^4 .....*/
           REAL t[],         /* Parameters nodes .....................*/
           REAL xt[],        /* normalized tangent vectors (x comp.) .*/
           REAL yt[]         /* normalized tangent vectors (y comp.) .*/
          );                 /* error code ...........................*/

/* ------------------------- END hrmsplin.h ------------------------- */
