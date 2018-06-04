/* ---------------------- DECLARATIONS rk_fehl.h -------------------- */

int rk_fehl   /* Runge-Kutta-Fehlberg method for 1st order DE systems */
           (
            REAL      *x,        /* initial / final x-value ..........*/
            REAL      xend,      /* desired x-value ..................*/
            int       n,         /* number of DEs ....................*/
            REAL      y[],       /* initial/final y-value ............*/
            dglsysfnk dgl,       /* righ thand side for the system ...*/
            REAL      *h,        /* initial/final step size ..........*/
            REAL      hmax,      /* maximal step size ................*/
            REAL      epsabs,    /* absolute error bound .............*/
            REAL      epsrel     /* relative error bound .............*/
           );                    /* error code .......................*/

/* -------------------------- END rk_fehl.h ------------------------- */
