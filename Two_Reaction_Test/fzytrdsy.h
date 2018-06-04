/*.BA*/



/*.FE{C 4.11.2}
     {Systems with Symm.\ Cycl.\ Trid.\ Str.\ Nonsing.\ Matrices}
     {Systems with Symmetric Cyclically
      Tridiagonal Strongly Nonsingular Matrices}*/

/*.BE*/
/* -------------------- DECLARATIONS fzytrdsy.h --------------------- */

int zytrdsy    /* cyclic tridiagonal symmetric linear system .........*/
           (
            int  modus,     /* Modus of call: 0, 1, 2 ................*/
            int  n,         /* size of matrix ........................*/
            REAL diag[],    /* main diagonal or D in R'*D*R ..........*/
            REAL oben[],    /* first codiagonal in R in R'*D*R .......*/
            REAL rechts[],  /* right most column of R in R'*D*R ......*/
            REAL rs[]       /* right hand side/solution ..............*/
           );               /* error code ............................*/

/* ------------------------- END fzytrdsy.h ------------------------- */
