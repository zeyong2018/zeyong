/***********************************************************************
*                                                                      *
* imported macros:                                                     *
* ================                                                     *
* GRAFIKGEWOLLT        to be defined by the user                       *
* __DJGPP__            standard macro of DJGPP 1 and 2                 *
* __PUREC__            standard macro of Pure C 1.0 and 1.1            *
* _QC                  standard macro of QuickC 2.0, 2.5 and MC 6.0    *
* __TURBOC__           standard macro of all Turbo/Borland compilers   *
* __OS2__              standard macro of all OS/2 compilers            *
*                                                                      *
* exported macros:                                                     *
* ================                                                     *
* DJGPP2               for recognition of DJGPP 2                      *
* INIT_FONTS           for Initialization of text output with DJGPP 2  *
* GRAFIK_H             H file of the BGI-compatible graphics library   *
* BIOS_H               H file for keyboard input (only BGI graphics)   *
* BGIGRAFMOEGLICH      for recognition of BGI-compatible graphics      *
* MCGRAFMOEGLICH       for recognition of QuickC-compatible graphics   *
* GRAFIKMOEGLICH       for recognition of any kind of graphics         *
* MITGRAFIK            for recognition if graphics works and is wanted *
*                                                                      *
* Author: Juergen Dietel, Computing Center, RWTH Aachen                *
* Date:   3.17.1996 - 2.26.1997                                        *
***********************************************************************/

#ifdef __DJGPP__
#if __DJGPP__ >= 2
#define DJGPP2
#endif
#endif

#ifdef DJGPP2
#define GRAFIK_H  <libbcc.h>
#define INIT_FONTS                          \
  {                                         \
    registerbgifont(&_litt_font);           \
    settextstyle(SMALL_FONT, HORIZ_DIR, 4); \
  }
#else
#define GRAFIK_H  <graphics.h>
#define INIT_FONTS
#endif

#ifdef __PUREC__
#define BIOS_H    <bios_pc.h>
#else
#define BIOS_H    <bios.h>
#endif

#if (defined(__TURBOC__) && ! defined(__OS2__)) || defined(DJGPP2)
#define BGIGRAFMOEGLICH
#endif

#ifdef _QC
#define MCGRAFMOEGLICH
#endif

#if defined(BGIGRAFMOEGLICH) || defined(MCGRAFMOEGLICH)
#define GRAFIKMOEGLICH
#endif

#if defined(GRAFIKMOEGLICH) && defined(GRAFIKGEWOLLT)
#define MITGRAFIK
#endif
