/*.BA*/
/*.KA{C 12}{Two-Dim., B\'ezier, Surface, B Splines}
           {Two-Dimensional Splines, Surface Splines, B\'ezier Splines,
            B Splines}*/
/*.BE*/
/* ---------------------- DECLARATIONS bikub.h ---------------------- */

/***********************************************************************
* include file for files bikub.c, bezier.c                             *
***********************************************************************/

#ifndef BIKUB_H_INCLUDED
#define BIKUB_H_INCLUDED

int  bikub1   (int n, int m, mat4x4** mat, REAL* x, REAL* y);
int  bikub2   (int n, int m, mat4x4** mat, REAL* x, REAL* y);
int  bikub3   (int n, int m, mat4x4** mat,
               REAL* x, REAL* y, REAL*** fn);
int  bsval  (int n, int m, mat4x4** mat,
             REAL* x, REAL* y, REAL xcoord,
             REAL ycoord, REAL* value);
int  xyintv (int n,  int m,  REAL* x,  REAL* y,
             int* i, int* j, REAL* xi, REAL* yj,
             REAL xcoord, REAL ycoord);

int kubbez     /* compute Bezier points of a bezier spline curve .....*/
        (
         REAL   *b[],                /* weight points ................*/
         REAL   *d[],                /* Bezier points ................*/
         int    m,                   /* number of curve pieces .......*/
         int    dim                  /* 2,3 for planar, spatial curve */
        );                           /* error code ...................*/

int valbez     /* evaluation of a Bezier spline curve ................*/
        (
         REAL   t,                   /* parameter value t from [0,1]  */
         int    m,                   /* number of curve pieces        */
         int    dim,                 /* 2,3 for planar, spatial curve */
         REAL   *b[],                /* Bezier points                 */
         REAL   *x,                  /* coordinates of curve point    */
         REAL   *y,
         REAL   *z
        );                           /* error code                    */

int  bezier(REAL*** b, REAL*** d, int typ, int m, int n, REAL eps);
int  rechvp(REAL*** b, int m, int n, REAL vp,
            int num, REAL *points[]);
int  rechwp(REAL*** b, int m, int n, REAL wp,
            int num, REAL *points[]);

int mokube     /* compute Bezier points on a interpolating curve .....*/
        (
         REAL   *b[],                /* weight points ................*/
         REAL   *d[],                /* Bezier points ................*/
         int    m,                   /* number of spline segments ....*/
         int    dim,                 /* 2,3 for planar, spatial curve */
         REAL   eps                  /* accuracy of interpolation ....*/
        );                           /* error code ...................*/
#endif

/* -------------------------- END bikub.h --------------------------- */
