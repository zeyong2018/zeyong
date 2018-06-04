#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

/*--------------------------------------------------------------------*/
/* Generate a random matrix                                           */
/*--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main (void)
{
  time_t t;
  int i, j, n;

  if (scanf("%d",&n) != 1)
  {
    printf ("Input must be a number\n");
    return 1;
  }
  time (&t);
  srand((unsigned)t);

  printf("%d \n",n);

  for (i = 0;i < n; i++)
  { for (j = 0; j < n; j++)
     printf("%10.5f ",(double)rand() / (double)rand());
    printf("\n");
  }
  return (0);
}
