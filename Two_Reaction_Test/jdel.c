#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/***********************************************************************
*                                                                      *
* This program deletes all files that are mentioned in the command line*
* of a program. These may include so-called reply files which are not  *
* deleted themselves, but rather those files named in the reply files. *
* In order to distinguish betwen standard files and reply files, the   *
* latter must carry the prefix `-f'.                                   *
*                                                                      *
* language: ANSI C                                                     *
* Compiler: Borland C++ 1.00 fuer OS/2                                 *
*           Borland C++ 2.0 oder 3.1                                   *
*           Pure C 1.1                                                 *
*           GNU CC 2.5.8 oder 2.6.3                                    *
*           IBM C++ Set 2.1                                            *
*                                                                      *
* Author: Juergen Dietel, Computer Center, RWTH Aachen                 *
* date:   3. 23. 1994, 8. 16. 1995                                     *
*                                                                      *
***********************************************************************/

#include <stdio.h>        /*  for  fprintf, stderr, FILE, fopen,      */
                          /*       NULL, fclose, fgets, FILENAME_MAX, */
                          /*       printf, size_t, remove, fgetc,     */
                          /*       fflush                             */
#include <string.h>       /*  for  strlen, memmove, strncmp           */
#include <stdlib.h>       /*  for  exit                               */

#if defined(__TURBOC__)
#if __TURBOC__ <= 0x0200                              /* Turbo C 2.0? */
#define FILENAME_MAX  80
#endif
#endif



/* ------------------------------------------------------------------ */

/***********************************************************************
* declare global variable for module                                   *
***********************************************************************/

static int
  Interaktiv = 0,            /* Flag, indicating that no file should  */
                             /* be deleted without an interactive     */
                             /* check                                 */
  Stumm      = 0;            /* Flag, indicating that no screen       */
                             /* output is needed                      */


/* ------------------------------------------------------------------ */

static void saeubern
    (
     char *text
    )

/***********************************************************************
* remove all space, tab or carriage return symbols at the start or     *
* finish of `text' statements.                                         *
*                                                                      *
* global variables :                                                   *
* =================                                                    *
* size_t, memmove, strlen                                              *
***********************************************************************/

{
  size_t laenge,        /* old and new length of `text'               */
         vorne,         /* number of such symbols at start of  `text' */
         hinten;        /* ditto at end of `text'                     */
  int    zeichen;       /* current symbol in `text'                   */

  laenge = strlen(text);
  vorne  = 0;
  while (vorne < laenge                     &&
         ((zeichen = text[vorne]) == ' '    ||
         zeichen                  == '\t'   ||
         zeichen                  == '\n'))
    vorne++;
  hinten = laenge;
  while (hinten > vorne                          &&
         ((zeichen = text[hinten - 1]) == ' '    ||
         zeichen                       == '\t'   ||
         zeichen                       == '\n'))
    hinten--;
  laenge -= vorne + (laenge - hinten);
  memmove(text, text + vorne, laenge);
  text[laenge] = '\0';
}



/* ------------------------------------------------------------------ */

static void loeschen
    (
     char *dateiname
    )

/***********************************************************************
*                                                                      *
* global variables:                                                    *
* ================                                                     *
* saeubern, strlen, Interaktiv, printf, fgetc, remove, Stumm,          *
* fprintf, stderr, fflush                                              *
***********************************************************************/

{
  int antwort,
      fehler   = 0;


  saeubern(dateiname);

  if (strlen(dateiname) > 0)         /* Does file name exist? */
  {
    if (! Stumm && ! Interaktiv)
      printf("Delete file %s \n", dateiname);

    if (Interaktiv)
    {
      printf("Delete file %s  (j/n)? ", dateiname);
      antwort = fgetc(stdin);
      fflush(stdin);
    }

    if (! Interaktiv || antwort != 'n')
      fehler = remove(dateiname);

    if (! Stumm && fehler == -1)
      fprintf(stderr, "file %s cannot be "
                      "deleted.\n",
                      dateiname);
  }
}



/* ------------------------------------------------------------------ */

static void dateien_bearbeiten
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
*                                                                      *
* global  variables:                                                   *
* =================                                                    *
* FILENAME_MAX, fopen, NULL, fprintf, stderr, fgets, loeschen          *
***********************************************************************/

{
  int  i;
  FILE *antwortdatei;
  char dateiname[FILENAME_MAX];


  for (i = 1; i < argc; i++)
    if (strncmp(argv[i], "-f", 2) == 0)
    {
      antwortdatei = fopen(argv[i] + 2, "r");
      if (antwortdatei == NULL)
        fprintf(stderr, "file %s cannot be "
                        "opened for "
                        "reading.\n",
                        argv[i] + 2);
      else
      {
        while (fgets(dateiname, (int)sizeof(dateiname),
               antwortdatei) != NULL)
          loeschen(dateiname);
        fclose(antwortdatei);
      }
    }
    else if (*argv[i] != '-')
      loeschen(argv[i]);
}



/* ------------------------------------------------------------------ */

static void hilfe_zeigen(void)

/***********************************************************************
*                                                                      *
* global variables :                                                   *
* =================                                                    *
* fprintf, stderr, exit                                                *
***********************************************************************/

{
  fprintf(stderr,
          "\n"
          "Purpose: delete all files, whose names appear in "
          "command line.\n"
          "        If file name has prefix -f, delete all files "
          "        listed rowwise in it instead.\n"
          "Use:    jdel [-i] [-s] [-f]filename [...]\n"
          "        -i  interactive\n"
          "        -s  silent\n"
         );
  exit(1);
}



/* ------------------------------------------------------------------ */

static void kommandozeile_analysieren
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
*                                                                      *
* global variables:                                                    *
* ================                                                     *
* hilfe_zeigen, strlen, Interaktiv, Stumm                              *
***********************************************************************/

{
  int i, j;


  if (argc <= 1)
    hilfe_zeigen();

  for (i = 1; i < argc; i++)
  {
    printf("argv[%d] = >>>%s<<<\n", i, argv[i]);
    if (*argv[i] == '-')
      if (strlen(argv[i]) >= 2)
        for (j = 1; j < strlen(argv[i]); )
          switch (argv[i][j])
          {
            case 'i': Interaktiv = 1;
                      j++;
                      break;
            case 's': Stumm = 1;
                      j++;
                      break;
            case 'f': j = (int)strlen(argv[i]);
                      break;
            default:  hilfe_zeigen();
          }
      else
        hilfe_zeigen();
  }
}



/* ------------------------------------------------------------------ */

int main
    (
     int  argc,
     char *argv[]
    )

/***********************************************************************
*                                                                      *
* globale variables:                                                   *
* ==================                                                   *
* kommandozeile_analysieren, dateien_bearbeiten                        *
***********************************************************************/

{
  kommandozeile_analysieren(argc, argv);

  dateien_bearbeiten(argc, argv);   /* work on all files in the       */
                                    /* command line                   */

  return 0;
}
