/*.BA*/



/*.FE{C 17.7.3}
     {Gear's method for integrating stiff systems of DEs}
     {Gear's method for integrating stiff systems of DEs}*/

/*.BE*/
/* ---------------------- DECLARATIONS gear.h ----------------------- */

int gear4          /* Gear  method of 4th order for DESs of 1st order */
    (
     REAL      *x,             /* starting or end point      .........*/
     REAL      xend,           /* desired end point (> x)     ........*/
     int       n,              /* number of DEs    ...................*/
     dglsysfnk dgl,            /* right hand side of system of DEs ...*/
     REAL      y[],            /* initial value or solution ..........*/
     REAL      epsabs,         /* absolute error bound    ............*/
     REAL      epsrel,         /* relative error bound    ............*/
     REAL      *h,             /* starting or final step size   ......*/
     long      fmax,           /* maximal number of calls of dgl() ...*/
     long      *aufrufe        /* actual number of calls of dgl() ....*/
    );                         /* error code .........................*/

char *gear_fehlertext          /* find error message and error class  */
    (
     int       fehlercode,     /* error number from gear()       .....*/
     fehler_t  *fehlerart      /* severity of error from gear()  .....*/
    );                         /* error text .........................*/

/* --------------------------- END  gear.h -------------------------- */
