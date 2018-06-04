/*.BA*/



/*.FE{C 17.3.4.4}{Embedding Formulas}{Embedding Formulas}*/

/*.BE*/
/* ---------------------- DECLARATIONS einb_rk.h -------------------- */

int einb_rk      /* solve DE system with one of 15 embedding formulas */
           (
            REAL      *x,        /* initial/final x-value ............*/
            REAL      beta,      /* desired final x-value ............*/
            int       n,         /* number of DEs ....................*/
            dglsysfnk dgl,       /* right hand side of system ........*/
            REAL      y[],       /* initial y-value at x .............*/
            REAL      abserr,    /* absolute error bound .............*/
            REAL      relerr,    /* relative error bound .............*/
            int       neinb,     /* Number of embedding formula ......*/
            int       hullstp,   /* step size control according to    */
                                 /* Hull ? ...........................*/
            int       neu,       /* do not use old data ? ............*/
            int       save,      /* save data for a future call ? ....*/
            long      fmax,      /* maximal # of calls of  dgl() .....*/
            long      *aufrufe   /* actual # of calls of dgl() .......*/
           );                    /* error code .......................*/

/* -------------------------- END einb_rk.h ------------------------- */
