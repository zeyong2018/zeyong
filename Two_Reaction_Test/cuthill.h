/*.BA*/



/*.FE{C 4.17}
     {The Cuthill-McKee Algorithm}
     {The Algorithm of Cuthill-McKee for Sparse Symmetric Matrices}*/

/*.BE*/
/* --------------------- DECLARATIONS cuthill.h --------------------- */

int cutgaucho       /* sparse matrices via  Cuthill-McKee ............*/
    (
     boolean gauss,       /* Flag: Gauss or Cholesky .................*/
     int     n,           /* size of the sparse matrix ...............*/
     int     nv,          /* number of nonzero matrix elements + n ...*/
     int     ic[],        /* vector with column indices of v-elements */
     REAL    v[],         /* vector with nonzero matrix elements .....*/
     int     nrs,         /* number of right hand sides ..............*/
     REAL    rs[],        /* vector with all right hand sides ........*/
     REAL    x[],         /* vector with all solutions ...............*/
     int     *m           /* half band width of matrix ...............*/
    );                    /* error code ..............................*/

/* ------------------------- END cuthill.h -------------------------- */
