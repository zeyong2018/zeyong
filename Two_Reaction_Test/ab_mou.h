/*.BA*/



/*.FE{C 17.4.3}
     {The Predictor-Corrector Method of Adams-Moulton}
     {The Predictor-Corrector Method of Adams-Moulton}*/

/*.BE*/
/* ---------------------- DECLARATIONS ab_mou.h --------------------- */

int prae_korr   /* Predictor-corrector meth. for 1st order DE systems */
             (
              REAL      *x,          /* initial/final x-value ........*/
              REAL      y[],         /* initial value/ solution ......*/
              int       n,           /* number of DEs ................*/
              dglsysfnk dgl,         /* right hand side for DEs ......*/
              REAL      xend,        /* desired final x-value ........*/
              REAL      *h,          /* starting/final step size .....*/
              REAL      epsabs,      /* absolute error bound .........*/
              REAL      epsrel,      /* relative error bound .........*/
              long      fmax,        /* maximal # of calls for dgl() .*/
              long      *aufrufe,    /* actual # of calls of dgl() ...*/
              REAL      hmax,        /* maximal step size ............*/
              boolean   neu          /* delete old data ? ............*/
             );                      /* error code ...................*/

/* --------------------------- END ab_mou.h ------------------------- */
