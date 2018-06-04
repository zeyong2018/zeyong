/*.BA*/
/*.KA{C 11}
     {Cubic Fitting Splines}
     {Cubic Fitting Splines for Constructing Smooth Curves}*/
/*.BE*/
/* ---------------------- DECLARATIONS glsp.h ----------------------- */

/***********************************************************************
* header file for routines for spline computation       Egg, 2.24.1991 *
***********************************************************************/
#ifndef GLSP_H_INCLUDED
#define GLSP_H_INCLUDED

int glspnp  (int   n,      REAL* x,         REAL* f,
             REAL* w,      int   marg_cond, REAL marg_0,
             REAL  marg_n, REAL* a,         REAL* b,
             REAL* c,      REAL* d);

int glsptr  (int   n,  REAL* x,         REAL* f,
             REAL* w,  int   Verschieb, REAL* px,
             REAL* py, REAL* A,         REAL* B,
             REAL* C,  REAL* D,         REAL* phi,
             REAL* r,  REAL* Dphi,      REAL* help);

int glsppa  (int   n,       REAL* x,       REAL* f,
             REAL* wx,      REAL* wf,      REAL* t,
             int   marke_t, int   rand,    REAL* alpha,
             REAL* beta,    int   marke_w, REAL* ax,
             REAL* bx,      REAL* cx,      REAL* dx,
             REAL* ay,      REAL* by,      REAL* cy,
             REAL* dy,      REAL* help);

int glsppe  (int   n,     REAL* x,     REAL* f,
             REAL* w,     int   rep,   REAL* a,
             REAL* b,     REAL* c,     REAL* d,
             REAL* h,     REAL* h1,    REAL* h2,
             REAL* h3,    REAL* rs,    REAL* hup);

int fzyfsy (int   n,
            REAL* md,     REAL* ud1,    REAL* ud2,
            REAL* rs,
            REAL* x,
            REAL* cmd,
            REAL* cld_1,  REAL* cld_2,
            REAL* cld_l2, REAL* cld_l1,
            REAL* bud_1,  REAL* bud_2,
            REAL* brs_2,  REAL* brs_1
           );

int fzyfsz (int   n,
            REAL* md,     REAL* ud1, REAL* ud2,
            REAL* cmd,
            REAL* cld_1,  REAL* cld_2,
            REAL* cld_l2, REAL* cld_l1,
            REAL* bud_1,  REAL* bud_2,
            REAL* brs_2,  REAL* brs_1
           );

int fzyfsl (int   n,
            REAL* rs,     REAL* x,
            REAL* cmd,
            REAL* cld_1,  REAL* cld_2,
            REAL* cld_l2, REAL* cld_l1,
            REAL* bud_1,  REAL* bud_2,
            REAL* brs_2,  REAL* brs_1
           );
#endif

/* --------------------------- END glsp.h --------------------------- */
