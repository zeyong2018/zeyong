/*.BA*/
/*.KA{C 17}{Initial Value Problems}
           {Initial Value Problems for Ordinary Differential
            Equations}*/
/*.FE{C 17.3}{One-Step Methods}{One-Step Methods}*/

/*.BE*/
/* ---------------------- DECLARATIONS einschr.h -------------------- */

int dglesv         /* One-step methods for 1st order DEs .............*/
          (
           REAL   *x,        /* initial/final x-value ................*/
           REAL   *y,        /* ditto for y ..........................*/
           dglfnk dgl,       /* righ thand side of DE ................*/
           REAL   xend,      /* desired final x-value ................*/
           REAL   *h,        /* initial/final step size ..............*/
           REAL   epsabs,    /* absolute error bound .................*/
           REAL   epsrel,    /* relative error bound .................*/
           int    intpol,    /* use interpolation ? ..................*/
           int    methode,   /* select method  (0, 1, 2) .............*/
           int    rand,      /* xend not reached ? ...................*/
           int    neu        /* pass on old data ? ...................*/
          );                 /* error code ...........................*/

/* -------------------------- END einschr.h ------------------------- */
