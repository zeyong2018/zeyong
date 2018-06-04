/*.BA*/



/*.FE{C 9.8.2}{Shepard Interpolation}{Shepard Interpolation}*/

/*.BE*/
/* --------------------- DECLARATIONS shepard.h --------------------- */

int shepard       /* Shepard interpolation (global,local, F-L weights)*/
    (
     REAL x0,           /* (x0,y0) = Interpolation point  ............*/
     REAL y0,
     REAL x[],          /* (x[i],y[i]) = nodes         ...............*/
     REAL y[],
     REAL f[],          /* f[i] = function values at nodes ...........*/
     int  n,            /* number of nodes  - 1         ..............*/
     REAL mue,          /* adjustable Shepard parameter (>0) .........*/
     int  methode,      /* global, local, F.-Little weights: (0,1,2)? */
     REAL R,            /* Radius for circle around (x0,y0) ..........*/
     REAL *PHI          /* Interpolation value at (x0,y0) ............*/
    );                  /* error code  ...............................*/

/* ------------------------- END  shepard.h ------------------------- */
