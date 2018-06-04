/* -------------------- DECLARATIONS nlglstst.h --------------------- */

#ifndef NLGLSTST_H_INCLUDED
#define NLGLSTST_H_INCLUDED

typedef struct
{
  int   n;                             /* number of equations         */
  nlgls fkt;                           /* one of the n functions      */
  char  *(*fkt_text)(void);            /* function definition as text */
} bsptyp1;

bsptyp1 *nlgls_waehlen(int nummer);

#endif

/* ------------------------- END nlglstst.h ------------------------- */
