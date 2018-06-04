/*.BA*/



/*.FE{C 17.3.5}
     {Implicit Runge-Kutta Methods of Gaussian Type}
     {Implicit Runge-Kutta Methods of Gaussian Type}*/

/*.BE*/
/* --------------------- DECLARATIONS implruku.h -------------------- */

int implruku  /* Implic. Runge-Kutta meth. for ODE syst. of 1st order */
            (
             dglsysfnk dgl,          /* right hand side of DE system .*/
             int       n,            /* number of DEs ................*/
             int       mmax,         /* max. order of IRKMs (>= 5) ...*/
             char      file_st[],    /* Name of coefficient file .....*/
             char      file_out[],   /* Name of output file ..........*/
             char      file_pro[],   /* Name of protocol file ........*/
             REAL      epsm,         /* Machine constant .............*/
             REAL      *eps_rel,     /* relative error bound or       */
                                     /* maximal relative error .......*/
             long      fmax,         /* max. number of calls of dgl() */
             long      *aufrufe,     /* act. number of calls of dgl() */
             REAL      g[],          /* weights for  y ...............*/
             REAL      *x0,          /* initial/final x-value ........*/
             REAL      xend,         /* desired final x-value ........*/
             REAL      y0[],         /* initial y-value at x0 ........*/
             REAL      yq[]          /* final y-value at xend ........*/
            );                       /* error code ...................*/

/* -------------------------- END implruku.h ------------------------ */
