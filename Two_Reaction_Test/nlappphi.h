/* -------------------- DECLARATIONS nlappphi.h --------------------- */

#ifndef NLAPPPHI_H_INCLUDED
#define NLAPPPHI_H_INCLUDED

typedef struct
{
  int       n;                      /* number of coefficients         */
                                    /* (0: arbitrarily many)          */
  approxfnk PHI;                    /* approximating function         */
  ableitfnk ABL;                    /* its partial derivatives        */
  char      *(*PHI_text)(void);     /* approximating function as text */
} bsptyp3;

bsptyp3 *nliansf_waehlen(unsigned int nummer);

#endif

/* ------------------------- END nlappphi.h ------------------------- */
