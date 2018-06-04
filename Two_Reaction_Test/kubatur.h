/*.BA*/
/*.KA{C 16}{Numerical Cubature}{Numerical Cubature}*/
/*.BE*/
/* --------------------- DECLARATIONS kubatur.h --------------------- */

#ifndef KUBATUR_H_INCLUDED
#define KUBATUR_H_INCLUDED

int  Kub4NeCn  (REAL a, REAL b, int Nx, REAL c, REAL d, int Ny,
                int Verfahren, REAL f(REAL,REAL),
                REAL* Wert, int Schaetzen, REAL* FehlerSch);
int  Kub3NeC3  (REAL Px, REAL Py, REAL Qx, REAL Qy, REAL Rx, REAL Ry,
                int n, REAL f(REAL,REAL), REAL* Wert);
int  Kub4RoRi  (REAL a, REAL b, int Nx, REAL c, REAL d, int Ny,
                int nST, REAL f(REAL,REAL),
                REAL* Wert, REAL* FehlerSch);
int  Kub4BuRi  (REAL a, REAL b, int Nx, REAL c, REAL d,
                int Ny, int nST, REAL f(REAL,REAL),
                REAL* Wert, REAL* FehlerSch);
int  Kub3RoRi  (REAL Px, REAL Py, REAL Qx, REAL Qy, REAL Rx, REAL Ry,
                int n, REAL f(REAL,REAL), REAL* Wert, REAL* FehlerSch);
int  Kub4GauE  (REAL a, REAL b, int Nx, REAL c, REAL d,
                int Ny, int Verf, REAL f(REAL,REAL),
                REAL* Wert, int Schaetzen, REAL* FehlerSch);
int  Kub4GauV  (REAL* x, int Nx, REAL* y, int Ny,
                int Verf, REAL f(REAL,REAL),
                REAL* Wert, int Schaetzen, REAL* FehlerSch);
int  Kub3GauN  (REAL Px, REAL Py, REAL Qx, REAL Qy,
                REAL Rx, REAL Ry, int n, int m,
                REAL f(REAL,REAL), REAL* Wert);

REAL fibiku (int n, int m, mat4x4** a, REAL* x, REAL* y);
int  fibik2 (int n, int m, mat4x4** a, REAL* x, REAL* y,
             REAL xlow, REAL ylow, REAL xup, REAL yup,
             REAL* value);
#endif

/* ------------------------- END kubatur.h -------------------------- */
