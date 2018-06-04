#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "u_proto.h"

static int isintab (char *name);

char *input = " < poly\\POLxx ";
/* UNIX
   char *input = " < poly/POLxx ";
 */
char output[20];
char prog[100];

char *testtab[] = { "tnewpol ",
                    "tmueller",
                    "tbauhube"
                  };

#define TESTNAMES (sizeof (testtab) / sizeof (testtab[0]))
#define MAXPOL  10

int main (int argc, char *argv[])
/*--------------------------------------------------------------------*/
/* Test program for  :            tnewpol                             */
/*                                tmueller                            */
/*                                tbauhube                            */
/* Call  : tprog 'progname'                                           */
/*                                                                    */
/* At most 10 polynomials are being tested, which reside in the       */
/* (relative) file  'poly' : POL01, ....., MAT10                      */
/*                                                                    */
/* The output of the tests is stored in the file  t'progname'         */
/*--------------------------------------------------------------------*/
{
  int i, numpol;
  char ti[3], * pos;

  if (argc < 2 || argv[1] == NULL)
  {
    printf ("Program name missing\n");
    return (1);
  }

  /* get the number of polynomials ...................................*/
  if (argv[2] == NULL) numpol = MAXPOL;
  else
    if (sscanf (argv[2], "%d", &numpol) <= 0)
      numpol = MAXPOL;

  numpol = (numpol > MAXPOL) ? MAXPOL : numpol;

  strcpy (prog, argv[1]);

  if (!isintab (argv[1]))
  {
    printf ("Inadmissable test procedure name\n");
    printf ("The following are inadmissable:\n");
    for (i = 0; i < TESTNAMES; i++)
      printf ("%s\n", testtab[i]);

    return (1);
  }

  strcpy (output, argv[1]);
  output[0] = 't';
  output[8] = '\0';

  /* delete output file ..............................................*/
  unlink (output);

  /* make full name ..................................................*/
  strcat (prog, input);
  strcat (prog, ">> ");
  strcat (prog, output);

  pos = strchr (prog, 'x');

  for (i = 1; i <= numpol; i++)  /* for all test polynomials .........*/
  {
    sprintf (ti, "%02d", i);
    pos[0] = ti[0];
    pos[1] = ti[1];
    printf ("%s\n", prog);
    if (system (prog) == -1) exit (1);
  }
  return (0);
}


static int isintab (char *name)
/*--------------------------------------------------------------------*/
/* Check whether 'name' is contained in testtab                       */
/* RETURN: 0 no matching entry, 1 matching entry                      */
/*--------------------------------------------------------------------*/
{
  int i;

  for (i = 0; i < TESTNAMES; i++)
    if (strcmp (name, testtab[i]) == 0) return (1);

  return (0);
}
