#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <string.h>

#include <basis.h>
#include <u_proto.h>

char *input = " < mat\\MATxx ";
char output[20];
char prog[100];

char *testtab[] = { "mgauss",
                    "mpivot",
                    "meigen"
                  };

#define TESTNAMES (sizeof (testtab) / sizeof (testtab[0]))

int main (int argc, char * argv[])
/*--------------------------------------------------------------------*/
/* Testing environment for the programs: mgauss                       */
/*                                       mpivot                       */
/*                                       meigen                       */
/*                                                                    */
/* Call: tprog 'progname'                                             */
/*                                                                    */
/* We test at most 47 matrices, which are stored in the file 'mat':   */
/* MAT01, ....., MAT47                                                */
/*                                                                    */
/* The output is stored in the file t'progname'                       */
/*--------------------------------------------------------------------*/
{
  int i, nummat;
  char ti[3], * pos;

  if (argc < 2 || argv[1] == NULL)
  {
    printf ("Program name missing\n");
    return (1);
  }

  /* find number of matrices ........................................ */
  if (argv[2] == NULL) nummat = 47;
    else nummat = atoi (argv[2]);

  nummat = (nummat > 47) ? 47 : nummat;

  strcpy (prog, argv[1]);
  if (!isintab (argv[1]))
  {
    printf ("Invalid test procedure name\n");
    printf ("valid names are:\n");
    for (i = 0; i < TESTNAMES; i++)
      printf ("%s\n", testtab[i]);

    return (1);
  }

  strcpy (output, argv[1]);
  output[0] = 't';
  output[8] = '\0';

  /* make full name ..................................................*/
  strcat (prog, input);
  strcat (prog, ">> ");
  strcat (prog, output);

  /* delete output file ..............................................*/
  unlink (output);
  pos = strchr (prog, 'x');

  for (i = 1; i <= nummat; i++)  /* for all test matrices ............*/
  {
    sprintf (ti, "%02d", i);
    pos[0] = ti[0];
    pos[1] = ti[1];
    printf ("%s\n", prog);
    if (system (prog) == -1) return (1);
  }
  return (0);
}


static int isintab (char *name)
/*--------------------------------------------------------------------*/
/* Checks whether 'name' exists inside testtab.                       */
/* Return value: 0 entry not found, 1 entry found                     */
/*--------------------------------------------------------------------*/
{
  int i;

  for (i = 0; i < TESTNAMES; i++)
    if (strcmp (name, testtab[i]) == 0) return (1);

  return (0);
}
