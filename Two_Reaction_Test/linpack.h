/* --------------------- DECLARATIONS linpack.h --------------------- */

int sspfa     /* decompose a condensed symmetric matrix ..............*/
          (
           REAL ap[],     /* upper triangle of matrix, condensed .....*/
           int  n,        /* size of matrix ..........................*/
           int  pvt[]     /* Pivot vector ............................*/
          );              /* singular pivot blocks ? .................*/

void sspsl    /* Solve linear system for a symmetric condensed matrix */
          (
           REAL ap[],      /* Vector with condensed factorization ....*/
           int  n,         /* size of matrix .........................*/
           int  pvt[],     /* Pivot indices ..........................*/
           REAL b[]        /* right hand side/solution vector ........*/
          );

int sspco    /* factor condensed symmetric matrix, estimate condition */
         (
          REAL ap[],     /* upper triangle of matrix, condensed ......*/
          int  n,        /* size of matrix ...........................*/
          int  pvt[],    /* Pivot indices ............................*/
          REAL *rcond    /* estimate for reciprocal of condition # ...*/
         );              /* error code ...............................*/

/* ------------------------- END linpack.h -------------------------- */
