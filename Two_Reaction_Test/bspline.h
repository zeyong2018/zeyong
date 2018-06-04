


/*.FE{C 12.4}{B Splines}{B Splines}*/

/* --------------------- DECLARATIONS bspline.h --------------------- */

int bspline     /* Compute points on a B spline curve ................*/
           (
            REAL *d[],     /* given de Boor points ...................*/
            int  n,        /* number of de Boor points ...............*/
            int  k,        /* Order of the B spline (3<=k<=n) ........*/
            int  m,        /* dimension of space .....................*/
            int  offen,    /* open curve? ............................*/
            REAL *c[],     /* computed points of the curve ...........*/
            int  *nc       /* maximal or actual number of points .....*/
           );              /* error code .............................*/

int bspline2    /* compute a series of w curves on a B spline surface */
            (
             REAL **d[],   /* given de Boor points ...................*/
             int  m,       /* number of de Boor points in v direction */
             int  n,       /* number of de Boor points in w direction */
             int  k,       /* order of the B spline ..................*/
             int  voffen,  /* open v curves? .........................*/
             int  woffen,  /* open w curves? .........................*/
             REAL **c[],   /* computed points on the surface .........*/
             int  *nv,     /* max./act. number of points in v direct. */
             int  *nw      /* max./act. number of points in w direct. */
            );             /* error code .............................*/

/* ------------------------- END bspline.h -------------------------- */
