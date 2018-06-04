/*.BA*/



/*.FE{C 4.7.2}
     {The Conjugate Gradient Method}{The Conjugate Gradient Method}*/

/*.BE*/
/* ----------------------- DECLARATIONS cg.h ------------------------ */

int cg_verfahren             /* Conjugate Gradient Method ............*/
                (
                 int  n,     /* Size of the linear system ............*/
                 REAL *a[],  /* System matrix ........................*/
                 REAL y[],   /* right hand side ......................*/
                 REAL x[]    /* solution vector ......................*/
                );           /* Error code ...........................*/

/* ---------------------------- END cg.h ---------------------------- */
