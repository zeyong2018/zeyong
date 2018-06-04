/*.BA*/



/*.FE{C 17.5}
     {Bulirsch-Stoer-Gragg Extrapolation}
     {Bulirsch-Stoer-Gragg Extrapolation}*/

/*.BE*/
/* --------------------- DECLARATIONS bulirsch.h -------------------- */

int bul_stoe       /* Extrapolation method for 1st order DE systems ..*/
    (
     REAL      *x,             /* initial x-value/ final x-value .....*/
     REAL      xend,           /* desired end point ..................*/
     int       n,              /* number of DEs ......................*/
     dglsysfnk dgl,            /* right hand side of DE system .......*/
     REAL      y[],            /* initial y-value/ final y-value .....*/
     REAL      epsabs,         /* absolute error bound ...............*/
     REAL      epsrel,         /* relative error bound ...............*/
     REAL      *h,             /* initial/final step size ............*/
     REAL      hmax,           /* maximal step size ..................*/
     int       neu,            /* use an outside given x ? ...........*/
     long      fmax,           /* maximal # of calls of  dgl() .......*/
     long      *aufrufe        /* actual # of calls of dgl() .........*/
    );                         /* error code .........................*/

/* -------------------------- END bulirsch.h ------------------------ */
