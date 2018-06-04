/*.BA*/



/*.FE{C 12.2}{Two-Dimensional Interpolating Surface Splines}
             {Two-Dimensional Interpolating Surface Splines}*/

/*.BE*/
/* -------------------- DECLARATIONS thinplat.h --------------------- */

int prob2    /* compute two-dimensional surface splines ..............*/
         (
          int  NX,         /* number of spline nodes .................*/
          REAL x[],        /* nodes: x-coordinates ...................*/
          REAL y[],        /*        y-coordinates ...................*/
          REAL z[],        /* values to be smoothed at (x[i],y[i]) ...*/
          int  M,          /* derivative order .......................*/
          REAL rho,        /* smoothing parameter (>= 0) .............*/
          REAL w[],        /* weights ................................*/
          REAL c[]         /* Spline coefficients ....................*/
         );                /* error code .............................*/

REAL apprx2    /* Compute functional value of a surface spline .......*/
           (
            REAL x0,           /* (x0,y0) = place for evaluation .....*/
            REAL y0,
            int  NX,           /* number of spline nodes .............*/
            int  M,            /* derivative order ...................*/
            REAL x[],          /* nodes: x-coordinates ...............*/
            REAL y[],          /*        y-coordinates ...............*/
            REAL c[]           /* Spline coefficients ................*/
           );                  /* error code .........................*/

void ekreistrafo       /* Transformation to unit circle ..............*/
                (
                 REAL x[],    /* original or transformed coordinates .*/
                 REAL y[],
                 int  n       /* number of transformed points ........*/
                );

/* ------------------------- END thinplat.h ------------------------- */
