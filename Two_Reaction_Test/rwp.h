/*.BA*/
/*.KA{C 18}{Boundary Value Problems}
           {Boundary Value Problems for Ordinary Differential
            Equations}*/
/*.FE{C 18.2}
     {Reduct.\ of Boundary Value Problems to Initial Value Pr.}
     {Reduction of Boundary Value Problems to Initial Value Problems}*/

/*.BE*/
/* ----------------------- DECLARATIONS rwp.h ---------------------- -*/

int rwp    /* Shooting method for boundary value problem of 1st order */
       (
        REAL      a,           /* left end point .....................*/
        REAL      b,           /* right end point ....................*/
        REAL      h,           /* starting step size .................*/
        REAL      y_start[],   /* initial approximation or solution ..*/
                               /* of initial value problem:  y(a) ....*/
        int       n,           /* number of differntial equations ....*/
        dglsysfnk dgl,         /* right hand side for the system .....*/
        rndbedfnk rand,        /* Function for the boundary conditions*/
        int       awpnumm,     /* desired number for IVP solver ......*/
        REAL      epsawp,      /* error bound for initial value       */
                               /* problem ............................*/
        REAL      epsrb,       /* error bound for boundary value      */
                               /* problem ............................*/
        long      fmax,        /* maximal number of calls of dgl() ...*/
        int       itmax,       /* maximal number of Newton iterations */
        int       *act_iter    /* actual number of Newton steps ......*/
       );                      /* error code .........................*/

/* --------------------------- END rwp.h ---------------------------- */
