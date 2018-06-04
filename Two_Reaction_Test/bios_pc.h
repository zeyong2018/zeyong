/* --------------------- DECLARATIONS bios_pc.h --------------------- */

/*
  Using Turbo C under MS-DOS bioskey() enables you to wait until
  a key is pressed without any influence by redirected standard input
  as with getch(). For Pure C bioskey() is here simulated by a macro.
  If the MiNT library is used, an additional macro for itoa() must be
  defined because itoa() is not declared anywhere else (but should be
  in stdlib.h according to ANSI C!).
*/

#ifndef BIOS_PC_H_INCLUDED
#define BIOS_PC_H_INCLUDED

#ifndef __MINT__            /* Pure C without additions?              */
#include <tos.h>            /*  for  Bconin                           */
#else                       /* Pure C expanded by MiNT library?       */
#include <osbind.h>         /*  for  Bconin                           */
#include <support.h>        /*  for  _itoa                            */
#define itoa  _itoa
#endif
#undef bioskey
#define bioskey(x)  Bconin(2)

#endif

/* ------------------------- END bios_pc.h -------------------------- */
