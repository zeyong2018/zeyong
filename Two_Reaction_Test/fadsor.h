/*.BA*/



/*.FE{C 5.6.2}
     {Adaptive SOR Method}
     {Estimate for the Optimal Relaxation Coefficient, an Adaptive
      SOR Method}*/

/*.BE*/
/* --------------------- DECLARATIONS fadsor.h ---------------------- */

int adsor                         /* adaptive SOR method .............*/
         (
          int  crit,              /* Convergence criterion (0,1,2,3) .*/
          int  n,                 /* size of matrix ..................*/
          REAL *mat[],            /* matrix ..........................*/
          REAL b[],               /* right hand side .................*/
          REAL *omega,            /* Relaxation coefficient ..........*/
          REAL x[],               /* solution ........................*/
          REAL residu[],          /* Residuum ........................*/
          int  *iter,             /* number of iterations ............*/
          int  l,                 /* number of steps before adapting  */
                                  /* new coefficient .................*/
          REAL eps,               /* error bound  ....................*/
          int  maxit,             /* Maximal number of iterations ....*/
          int  methode            /* method :  (0,1,2) ...............*/
         );                       /* error code ......................*/

/* -------------------------- END fadsor.h -------------------------- */
