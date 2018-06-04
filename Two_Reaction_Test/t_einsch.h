/* -------------------- DECLARATIONS t_einsch.h --------------------- */

#ifndef T_EINSCH_H_INCLUDED
#define T_EINSCH_H_INCLUDED

typedef struct
{
  dglfnk y_strich;                         /* right hand side,        */
  char   *(*dgl_text)(void);               /* right hand side as text */
  REAL   (*y_exakt)(REAL x0,               /* and analytic solution   */
                    REAL y0,
                    REAL x
                   );
} bsptyp;

bsptyp *dgl_waehlen(int nummer);

#endif

/* ------------------------- END t_einsch.h ------------------------- */
