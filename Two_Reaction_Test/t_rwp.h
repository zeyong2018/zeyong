/* ---------------------- DECLARATIONS t_rwp.h ---------------------- */

#ifndef T_RWP_H_INCLUDED
#define T_RWP_H_INCLUDED

typedef struct
{
  int       n;                         /* number of equations,        */
  dglsysfnk rechte_seite;              /* right hand side,            */
  char      *(*dgl_text)(void);        /* right hand side as text     */
  void      (*rand)(REAL ya[],         /* boundary conditions         */
                    REAL yb[],
                    REAL f[]
                   );
  char      *(*randtext)(void);        /* boundary conditions as text */
} bsptyp;

bsptyp *rwp_waehlen(int nummer);

#endif

/* -------------------------- END t_rwp.h --------------------------- */
