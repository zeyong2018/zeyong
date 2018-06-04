/* ----------------------- DECLARATIONS fft.h ----------------------- */

#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED

typedef struct { REAL x; REAL y; } complex;
/*.IX{complex}*/


int rfft      /* fast real Fourier transform .........................*/
    (
     int  tau,            /* 2^tau = number of nodes .................*/
     REAL y[],            /* nodes / Fourier coefficients ............*/
     int  synthese        /* direction of transform ..................*/
    );                    /* error code ..............................*/

int fft       /* fast complex Fourier transform ......................*/
    (
     int     tau,         /* 2^tau = number of nodes .................*/
     complex y[],         /* node vector or Fourier coefficients .....*/
     int     synthese     /* direction of transform ..................*/
    );                    /* error code ..............................*/

int fftb      /* complex FFT for arbitary many nodes .................*/
    (
     int     N,           /* number of nodes .........................*/
     complex y[],         /* node vector or Fourier coefficients .....*/
     int     synthese     /* direction for transformation ............*/
    );                    /* error code ..............................*/

int fdicht    /* trigonomertric interpolating polynomial at shifted   */
              /* nodes ...............................................*/
    (
     int     M,           /* length of the table of values ...........*/
     complex F[],         /* values for f ............................*/
     REAL    p,           /* Period of f .............................*/
     REAL    theta        /* size of shift ...........................*/
    );                    /* error code ..............................*/

int fourn     /* Fourier transform for non-periodic function .........*/
    (
     int     M,           /* length of the table of values ...........*/
     complex F[],         /* values for f ............................*/
     REAL    a,           /* left edge of node interval ..............*/
     REAL    deltax       /* node distance ...........................*/
    );                    /* errorcode ...............................*/

int ffako     /* discrete cyclic convolution or correlation, periodic */
    (
     int     M,           /* length of the table of values ...........*/
     complex F[],         /* values for f ............................*/
     complex H[]          /* function values for h ...................*/
    );                    /* errorcode ...............................*/

int ffakon    /* discrete cyclic convolution or correlation,          */
              /* non-periodic ........................................*/
    (
     int     M,           /* number - 1 of f values ..................*/
     complex F[],         /* values for f ............................*/
     int     N,           /* number - 1 of h values ..................*/
     complex H[],         /* function values for h ...................*/
     int     tau,         /* L = 2^tau = length of F and H ...........*/
     REAL    deltax,      /* node distance ...........................*/
     int     faltung      /* Compute a convolution or correlation ? ..*/
    );                    /* errorcode ...............................*/

int nli_happr    /* non-linear least squares .........................*/
             (
              int       m,         /* number of nodes ................*/
              int       n,         /* number of coefficients - 1 .....*/
              REAL      x[],       /* x-values .......................*/
              REAL      y[],       /* y-values .......................*/
              REAL      w[],       /* positive weights ...............*/
              approxfnk PHI,       /* approximating function .........*/
              int       ablOK,     /* flag for derivatives via ABL ...*/
              ableitfnk ABL,       /* partial derivative of PHI ......*/
              int       *maxIt,    /* max/current number of iterations*/
              REAL      RelEps,    /* relative error bound ...........*/
              REAL      c[],       /* Starting/optimal coefficient ...*/
              REAL      *MiQuFe    /* mean square error ..............*/
             );                    /* error code .....................*/

#endif

/* --------------------------- END fft.h ---------------------------- */
