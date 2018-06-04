#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>         /*  for  umleiten, fprintf, stderr, scanf, */
                           /*       printf, REAL, LZS, fehler_melden, */
                           /*       LZP                               */
#include <t_einsch.h>      /*  for  bsptyp, dgl_waehlen               */
#include <einschr.h>       /*  for  dglesv                            */



int main(int argc, char *argv[])

{
  int    methode,
         intpol,
         rand,
         neu,
         weiter,
         bspnummer;
  REAL   x,
         y,
         h,
         xend,
         eps_abs,
         eps_rel,
         x0,                                      /* initial value    */
         y0;                                      /* y(x0) = y0       */
  int    fehler;
  bsptyp *beispiel;


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard ones           */
    return fehler;                                          /* 1 or 2 */

#ifdef INTERAKTIV
  fprintf(stderr, "\nExample: ");
#endif
  scanf("%d", &bspnummer);
  if ((beispiel = dgl_waehlen(bspnummer)) == NULL)
  {
    fehler_melden("non existing example number",
                  0, __FILE__, __LINE__);
    return 3;
  }

#ifdef INTERAKTIV
  fprintf(stderr, "\n\nMethod (2, 3, 5) = ");
#endif
  scanf("%d", &methode);

#ifdef INTERAKTIV
  fprintf(stderr, "\nbound (0, 1) = ");
#endif
  scanf("%d", &rand);


#ifdef INTERAKTIV
  fprintf(stderr, "\ninitial x-value  x = ");
#endif
  scanf("%"LZS"f", &x);
  x0 = x;

#ifdef INTERAKTIV
  fprintf(stderr, "\ninitial y-value  y = ");
#endif
  scanf("%"LZS"f", &y);
  y0 = y;

#ifdef INTERAKTIV
  fprintf(stderr, "\ninitial step size  h = ");
#endif
  scanf("%"LZS"f", &h);
#ifdef INTERAKTIV
  fprintf(stderr, "\n");
#endif

  printf("\n"
         "Solve an ordinary DE using a one-step method\n"
         "============================================\n\n"
         "DE:         y'  =  %s\n\n"
         "Method  = %22d\n"
         "bound   = %22d\n"
         "x       = %22.15"LZP"e\n"
         "y       = %22.15"LZP"e\n"
         "h       = %22.15"LZP"e\n",
         (*beispiel->dgl_text)(), methode, rand, x, y, h);

  do
  {
#ifdef INTERAKTIV
    fprintf(stderr, "\nfinal desired x-value  xend = ");
#endif
    scanf("%"LZS"f", &xend);

#ifdef INTERAKTIV
    fprintf(stderr, "\neps_abs = ");
#endif
    scanf("%"LZS"f", &eps_abs);

#ifdef INTERAKTIV
    fprintf(stderr, "\neps_rel = ");
#endif
    scanf("%"LZS"f", &eps_rel);

#ifdef INTERAKTIV
    fprintf(stderr, "\nFlag  intpol = ");
#endif
    scanf("%d", &intpol);

#ifdef INTERAKTIV
    fprintf(stderr, "\nusing old data?  "
                    "new (0, 1) = ");
#endif
    scanf("%d", &neu);

    printf("\n"
           "xend    = %22.15"LZP"e\n"
           "eps_abs = %22.15"LZP"e\n"
           "eps_rel = %22.15"LZP"e\n"
           "intpol  = %22d\n"
           "new     = %22d\n",
           xend, eps_abs, eps_rel, intpol, neu);

    fehler = dglesv(&x, &y, beispiel->y_strich, xend, &h, eps_abs,
                    eps_rel, intpol, methode, rand, neu);

    printf("\n"
           "error code: error             = %22d\n"
           "final x-value            x    = %22.15"LZP"e\n"
           "final y-value            y    = %22.15"LZP"e\n"
           "exact y-value            y    = %22.15"LZP"e\n"
           "final step size          h    = %22.15"LZP"e\n",
           fehler, x, y, (*beispiel->y_exakt)(x0, y0, x), h);

#ifdef INTERAKTIV
    fprintf(stderr, "\nContinue ? (1=yes, 0=no) ");   /* stop program */
#endif
    scanf("%d", &weiter);

  }
  while (weiter);
#ifdef INTERAKTIV
  fprintf(stderr, "\n");
#endif


  return 0;
}
