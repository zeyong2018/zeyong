/*.BA*/



/*.FE{}{The Cholesky Decomposition for Band Matrices}
       {The Cholesky Decomposition for Band Matrices}*/

/*.BE*/
/* --------------------- DECLARATIONS choband.h --------------------- */

int chobnd         /* Cholesky method for condensed band matrices ....*/
          (
           int  modus,    /* type of call: 0, 1, 2 ...................*/
           int  n,        /* size of the matrix ......................*/
           int  m,        /* half band width .........................*/
           REAL *ap[],    /* condensed matrix: Input or factorization */
           REAL rs[]      /* right hand side or solution .............*/
          );              /* Error code ..............................*/

/* ------------------------- END choband.h -------------------------- */
