/* -------------------- DECLARATIONS zeigkrv2.h --------------------- */

int zeigkrv2           /* plot a function table in R2 ................*/
            (
             int  nk,              /* size of function table .........*/
             REAL *kurve[],        /* function table .................*/
             int  ns,              /* number of nodes ................*/
             REAL *stuetz[]        /* nodes ..........................*/
            );                     /* error code .....................*/

int zeigflaeche      /* plot a function table in R3 ..................*/
        (
         REAL **c[],               /* function table .................*/
         int  nv,                  /* The table is a                  */
         int  nw,                  /* [0..nv-1,0..nw-1] matrix. ......*/
         REAL **d[],               /* nodes ..........................*/
         int  m,                   /* The nodes form a                */
         int  n,                   /* [0..m-1,0..n-1] matrix. ........*/
         int  voffen,              /* open v curves? .................*/
         int  woffen,              /* open w curves? .................*/
         int  st_gitter            /* grid through nodes? ............*/
        );                         /* error code .....................*/

/* ------------------------- END zeigkrv2.h ------------------------- */
