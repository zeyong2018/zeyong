/*.BA*/



/*.FE{C 4.10.2}
     {Systems with Trid.\ Symm.\ Strongly Nonsing.\ Matrices}
     {Systems with Tridiagonal Symmetric
      Strongly Nonsingular Matrices}*/

/*.BE*/
/* -------------------- DECLARATIONS ftrdiasy.h --------------------- */

int trdiasy      /* tridiagonal, symmetric linear system .............*/
           (
            int  modus,   /* kind of call: 0, 1, 2 ...................*/
            int  n,       /* size of the system ......................*/
            REAL diag[],  /* main diagonal of A, or D in A = R'DR ....*/
            REAL oben[],  /* co-diagonal of A, or co-diagonal of R ...*/
            REAL rs[]     /* right hand side, or solution ............*/
           );             /* error code ..............................*/

/* ------------------------- END ftrdiasy.h ------------------------- */
