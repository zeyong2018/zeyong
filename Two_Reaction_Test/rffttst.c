#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>     /*  for  PI, TRUE, FALSE, umleiten, SIN, REAL, */
                       /*       printf, scanf, NULL, fprintf, stderr, */
                       /*       FABS, REAL, TWO, fehler_melden, LZP   */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR                                */
#include <fft.h>       /*  for  rfft                                  */



/***********************************************************************
* Test program for  rfft()  , the fast real Fourier transform (FFT)    *
*                                                                      *
* Scope of the program:                                                *
* =====================                                                *
* ("pi" always denotes the number 3.14 ... here.)                      *
* The program reads tau and k from the first file in the command line. *
* It then computes the functional values of y(x) = sin(k * x) at       *
* equidistant nodes in [0,2pi). Subsequently rfft() computes the       *
* Fourier coefficients of the trigonometric interpolating polynomial   *
* (Fourier analysis). Here the (2*k)th coefficient must be one with all*
* others equal to zero. A second call of rfft() (Fourier synthesis)    *
* then reverses the transformation and computes the functional values  *
* from the computed Fourier coefficients. These ought to be the given  *
* equidistant ones. If this is not so, we print out the reconstructed  *
* function value and its distance to the original one.                 *
* The output data is gathered in the second file mentioned in the      *
* command line. If this is missing, all output is dumped to the screen.*
* If the first file is missing, the input must come from the keyboard. *
* Several test data sets are available in rffttst.ei*.                 *
* The output should coincide with the files in rffttst.au* .           *
*                                                                      *
* Construction of an input file:                                       *
* ==============================                                       *
* tau    Exponent for the number  2 ^ tau of nodes                     *
* k      Frequency factor for  sin(k * x)                              *
***********************************************************************/

int main(int  argc,
         char *argv[]
        )

{
  int  N;             /* number of nodes or coefficients              */
  int  tau,           /* 2 ^ tau = N                                  */
       j,             /* Loop variable                                */
       k,             /* Frequency factor in  sin(k * x)              */
       synthese,      /* Flag for rfft(), determining direction of    */
                      /* transform                                    */
       fehler;        /* error code of rfft()                         */
  REAL *y,            /* function values or  Fourier coefficients     */
       *rett,         /* Vector of functional values for later        */
                      /* comparison                                   */
       faktor,        /* 2 * k * pi / N  for computing function values*/
                      /* of sin(k * x) at equidistant nodes in [0,2pi)*/
       eps;           /* accuracy bound for comparing original and    */
                      /* recomputed function values                   */
  void *vmblock;      /* List of dynamically allocated vectors        */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "tau (number N of nodes = 2^tau):       ");
#endif
  scanf("%d", &tau);
#ifdef INTERAKTIV
  fprintf(stderr, "k   (Frequency factor for  sin(k * x): ");
#endif
  scanf("%d", &k);


  /* -------------- initialize accuracy bound ----------------------- */

#if defined(__GNUC__) || defined(__PUREC__) || \
   (defined(__TURBOC__) && (defined(__OS2__) || defined(__WIN32__)))
#if defined(__GNUC__) && defined(__MINT__) && ! defined(FLOAT)
  eps = (REAL)1.1e6 * MACH_EPS;
#else
  eps = (REAL)55000.0 * MACH_EPS;
#endif
#else
  eps = (REAL)5500.0 * MACH_EPS;
#endif


  /* ------------------ print out input data ------------------------ */

  printf("\n"
         "fast real Fourier transform (FFT)\n"
         "=================================\n\n"
         "tau = %3d    (number of nodes = 2^tau)\n"
         "k   = %3d    (Frequency factor for function  sin(k * x))"
         "\n",
         tau, k
        );


  if (tau > 8 * sizeof(int) - 2)        /* tau too large? (overflow!) */
  {
    fehler_melden("tau zu gross", 0, __FILE__, __LINE__);
    return 3;
  }

  N = 1 << tau;                                         /* 2 ^ tau    */

  vmblock = vminit();                            /* allocate storage  */
  y    = (REAL *)vmalloc(vmblock, VEKTOR, N, 0);
  rett = (REAL *)vmalloc(vmblock, VEKTOR, N, 0);
  if (! vmcomplete(vmblock))                     /* no success ?      */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

  faktor = TWO * k * PI / N;              /* generate N equidistant   */
  for (j = 0; j < N; j++)                 /* nodes of    sin(k * x)   */
    rett[j] = y[j] = SIN(faktor * j);     /* in the interval [0,2pi)  */
                                          /* for x = j*2*pi/N, j= ..  */


  /* ------- compute coefficients of trig interpolating ------------- */
  /* ------- polynomial (Fourier analysis)              ------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "\n");
#endif
  fprintf(stderr, "real Fourier analysis...");
  fehler = rfft(tau, y, synthese = FALSE);
  fprintf(stderr, " done\n");

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("rfft(): tau < 1", 10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("rfft(): tau too large",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ----------- print results -------------------------------------- */

  printf("\n\n"
         "discrete Fourier analysis for  f(x) = sin(k * x)  in [0,2pi):"
         "\n"
         "-------------------------------------------------------------"
         "\n"
        );
  for (j = 0; j < N; j++)
    if (FABS(y[j]) > eps)
      printf("non-zero Fourier coefficient  y[%5d] = "
             "%- 14.5"LZP"e\n", j, y[j]);


  /* ------- recreate original function values from ----------------- */
  /* ------- Fourier coefficients (Fourier synthesis)  -------------- */

  fprintf(stderr, "real Fourier synthesis...");
  fehler = rfft(tau, y, synthese = TRUE);
  fprintf(stderr, " done\n");

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("rfft(): tau < 1", 10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("rfft(): tau too large",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ------------ report results of synthesis ----------------------- */

  printf("\n\n"
         "discrete Fourier synthesis for  f(x) = sin(k * x)  "
         "in [0,2pi):\n"
         "---------------------------------------------------"
         "-----------\n");
  for (j = 0; j < N; j++)
    if (FABS(y[j] - rett[j]) > eps)
      printf ("Deviation %12.3"LZP"e for the synthesized "
              "functional value  y[%5d] = %- 14.5"LZP"e\n",
              FABS(y[j] - rett[j]), j, y[j]);


   synthese = synthese;                              /* calm compiler */

   return 0;
}
