#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>     /*  for  PI, TRUE, FALSE, umleiten, SIN, COS,  */
                       /*       REAL, printf, scanf, NULL, fprintf,   */
                       /*       stderr, FABS, REAL, TWO,              */
                       /*       fehler_melden, LZP                    */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VVEKTOR                               */
#include <fft.h>       /*  for  complex, fftb                         */



/***********************************************************************
* This is a test program for the function fftb() for computing the     *
* fast complex Fourier transform (FFT) for an arbitrary number of      *
* nodes.                                                               *
*                                                                      *
* Mode of operation:                                                   *
*                                                                      *
* ( i denotes the imaginary unit, pi = 3.14... throughout.)            *
* Read the numbers N and k from the first mentioned file. Compute the  *
* values of  exp(k * i * x) at N equidistant nodes in [0,2pi[. Compute *
* Fourier coefficients of the trigonometric interpolating polynomial   *
* via fftb().                                                          *
* For this example the kth coefficient must be one and all others zero.*
* A second call of fftb() recreates functional values from the Fourier *
* coefficients, which must be close to the original input.             *
* If this is not the case, we compute the norm of the difference.      *
* The second mentioned file serves as output file. If there is only one*
* file, the output is sent to the screen. If even the input file name  *
* is missing, the input data must be supplied from the keyboard.       *
* Several test sets are included in the file  fftbtst.ei*. The corres- *
* ponding output can be found in fftbtst.au* .                         *
*                                                                      *
* Input file:                                                          *
* ==========                                                           *
* N      number of nodes                                               *
* k      exponent to be used in  exp(k * i * x)                        *
***********************************************************************/

int main(int  argc,
         char *argv[]
        )

{
  int     N;          /* number of nodes and Fourier coefficients     */
  int     j,          /* variable                                     */
          k,          /* Exponent used in  exp(k * i * x)             */
          synthese,   /* Flag for fftb(), to specify whether to       */
                      /* compute Fourier coefficients or function     */
                      /* values                                       */
          fehler;     /*  error code from fftb()                      */
  complex *y,         /* functional values or Fourier coefficients    */
          *rett;      /* Vector with original function values         */
  REAL    faktor,     /* step size 2 * k * pi / N  used to evaluate   */
                      /* exp(k * i * x)  at  N equidistant nodes in   */
                      /* [0,2pi[                                      */
          eps;        /* error bound                                  */
  void    *vmblock;   /* List of dynamically allocated vectors        */
#ifdef DEBUG
  REAL    maxabw;
#endif


  if ((fehler = umleiten(argc, argv))   /* assign potential input and */
      != 0)                             /* output files to standard   */
    return fehler;  /* 1 oder 2 */      /* ones.                      */


  /* -------------------- read input -------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "N (Number of nodes):             ");
#endif
  scanf("%d", &N);
#ifdef INTERAKTIV
  fprintf(stderr, "k   (Expontent used in  exp(k * i * x):  ");
#endif
  scanf("%d", &k);


  /* -------------- set error bound        -------------- */

#if defined(__PUREC__)
  eps = (REAL)110000.0 * MACH_EPS;
#elif defined(__GNUC__) && defined(atarist) && ! defined(FLOAT)
  eps = (REAL)400000000.0 * MACH_EPS;
#else
  eps = (REAL)40000.0 * MACH_EPS;
#endif


  /* ------------------ print input data ------------------- */

  printf("\n"
         "complex fast Fourier transform (FFT)"
         "for arbitrarily many node\n"
         "============================================="
         "================================\n\n"
         "N   = %3d    (number of nodes)\n"
         "k   = %3d    (exponent of exp(k * i * x))  "
         "\n",
         N, k
        );


  vmblock = vminit();                 /* initialize storage */
  y    = (complex *)vmalloc(vmblock, VVEKTOR, N, sizeof(*y));
  rett = (complex *)vmalloc(vmblock, VVEKTOR, N, sizeof(*rett));
  if (! vmcomplete(vmblock))   /* storage demand not  */
  {                            /*     satisfied?              */
    fehler_melden("lack of available storage", 0, __FILE__, __LINE__);
    return 3;
  }

  faktor = TWO * k * PI / N;              /* N equidistant nodes in   */
  for (j = 0; j < N; j++)                 /* [0,2pi[; compute function*/
    rett[j].x = y[j].x = COS(faktor * j), /* values for exp(k*i*x)    */
    rett[j].y = y[j].y = SIN(faktor * j); /*  (x = j*2*pi/N)          */


  /* ------- compute the coeficients of the trigonomotric      ------ */
  /* ------- interpolating polynomial (Fourier analysis)       ------ */

#ifdef INTERAKTIV
  fprintf(stderr, "\n");
#endif
  fprintf(stderr, "complex Fourier analysis "
                  "for arbitrarily many nodes");
  fehler = fftb(N, y, synthese = FALSE) ;
  fprintf(stderr, " done\n");

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("fftb(): N < 2", 10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("fftb(): N too large",
                    10 + fehler, __FILE__, __LINE__);
    case 3:
      fehler_melden("fftb(): memory exhausted",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ----------- put out results                          ----------- */

  printf("\n\n"
         "discrete Fourier analysis of  f(x) = exp(k * i * x)  "
         "in [0,2pi]:\n"
         "----------------------------------------------------"
         "-----------\n");
  for (j = 0; j < N; j++)
    if (FABS(y[j].x) + FABS(y[j].y) > eps)
      printf("non vanishing Fourier koefficient  y[%5d] = "
             "(%- 14.5"LZP"e,%- 14.5"LZP"e)\n", j, y[j].x, y[j].y);


  /* ------- reconstruct functional values from Fourier       -------*/
  /* ------- coefficients (Fourier synthesis)                ------- */

  fprintf(stderr, "complex Fourier synthesis for "
                  "arbitrarily many coefficients...");
  fehler = fftb(N, y, synthese = TRUE);
  fprintf(stderr, " done\n");

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("fftb(): N < 1", 10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("fftb(): N too large",
                    10 + fehler, __FILE__, __LINE__);
    case 3:
      fehler_melden("fftb(): lack of memory",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ----------- put out results from synthesis           ---------- */

  printf("\n\n"
         "discrete Fourier synthesis of  f(x) = exp(k * i * x)  "
         "in [0,2pi]:\n"
         "-----------------------------------------------------"
         "-----------\n");
#ifdef DEBUG
  maxabw = ZERO;
#endif
  for (j = 0; j < N; j++)
  {
    if (FABS(y[j].x - rett[j].x) + FABS(y[j].y - rett[j].y) > TEN * eps)
      printf ("Deviation of %12.3"LZP"e when  synthesizing the "
              "function value y[%5d] = (%- 14.5"LZP"e,%- 14.5"LZP"e)\n",
              FABS(y[j].x - rett[j].x) + FABS(y[j].y - rett[j].y), j,
              y[j].x, y[j].y);
#ifdef DEBUG
    if (FABS(y[j].x - rett[j].x) + FABS(y[j].y - rett[j].y) > maxabw)
      maxabw = FABS(y[j].x - rett[j].x) + FABS(y[j].y - rett[j].y);
#endif
  }

#ifdef DEBUG
  printf("maximale Abweichung = %12.3"LZP"e\n", maxabw);
#endif


  synthese = synthese;                           /* calm the compiler */

  return 0;
}
