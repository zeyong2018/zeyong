/* -------------------- DECLARATIONS happrans.h --------------------- */

#ifndef HAPPRANS_H_INCLUDED
#define HAPPRANS_H_INCLUDED

typedef struct
{
  int       n;                            /* number of test functions */
                                          /* (0: arbitrarily many)    */
  ansatzfnk phi;                          /* test function i          */
  char      *(*phi_text)(void);           /* model function as text   */
} bsptyp2;

bsptyp2 *linansf_waehlen(unsigned int nummer);

#endif

/* ------------------------- END happrans.h ------------------------- */
