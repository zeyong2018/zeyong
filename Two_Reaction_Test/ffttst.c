#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <basis.h>     /*  for  PI, TRUE, FALSE, umleiten, SIN, COS,  */
                       /*       REAL, printf, scanf, NULL, fprintf,   */
                       /*       stderr, FABS, REAL, TWO,              */
                       /*       fehler_melden, LZP                    */
#include <vmblock.h>   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VVEKTOR                               */
#include <fft.h>       /*  for  complex, fft                          */



/***********************************************************************
* This is a test program for  fft()  for the fast complex Fourier      *
* transform (FFT).                                                     *
*                                                                      *
* Scope of the program :                                               *
* ======================                                               *
* ("i" shall denote the imaginary unit and "pi" the number 3.14...)    *
* The program reads tau and k from the first file in the command line. *
* It then computes the functional values of y(x) = exp(k * i * x) at   *
* equidistant nodes in [0,2pi). Subsequently  fft() computes the       *
* Fourier coefficients of the trigonometric interpolating polynomial   *
* (Fourier analysis). Here the kth coefficient must be one with all    *
* others equal to zero. A second call of fft() (Fourier synthesis)     *
* then reverses the transformation and computes the functional values  *
* from the computed Fourier coefficients. These ought to be the given  *
* equidistant ones. If this is not so, we print out the reconstructed  *
* function value and its distance to the original one.                 *
* The output data is gathered in the econd file mentioned in the       *
* command line. If this is missing, all output is dumped to the screen.*
* If the first file is missing, the input must come from the keyboard. *
* Several test data sets are available in  ffttst.ei*.                 *
* The output should coincide with the files in  ffttst.au* .           *
*                                                                      *
* Construction of an input file:                                       *
* ==============================                                       *
* tau    Exponent for the number  2 ^ tau of nodes                     *
* k      Frequency factor for exp(k * i * x)                           *
***********************************************************************/

int main(int  argc,
         char *argv[]
        )

{
  int     N;          /* number of nodes or coefficients              */
  int     tau,        /* 2 ^ tau = N                                  */
          j,          /* Loop variable                                */
          k,          /* Frequency factor in  exp(k * i * x)          */
          synthese,   /* Flag for fft(), determining direction of     */
                      /* transform                                    */
          fehler;     /* error code of  fft()                         */
  complex *y,         /* function values or  Fourier coefficients     */
          *rett;      /* Vector of functional values for later        */
                      /* comparison                                   */
  REAL    faktor,     /* 2 * k * pi / N  for computing function values*/
                      /* of  exp(k * i * x)  at  equidistant nodes    */
                      /* in [0,2pi)                                   */
          eps;        /* accuracy bound for comparing original and    */
                      /* recomputed function values                   */
  void    *vmblock;   /* List of dynmically allocated vectors         */


  if ((fehler = umleiten(argc, argv))   /* assign input/output files  */
      != 0)                             /* to standard files          */
    return fehler;  /* 1 or 2 */


  /* -------------------- read input data --------------------------- */

#ifdef INTERAKTIV
  fprintf(stderr, "tau (number N of nodes: N = 2^tau):        ");
#endif
  scanf("%d", &tau);
#ifdef INTERAKTIV
  fprintf(stderr, "k   (Frequency factor for  exp(k * i * x): ");
#endif
  scanf("%d", &k);


  /* -------------- set accuracy bound ------------------------------ */

#if defined(__GNUC__) || defined(__IBMC__) || defined(__IBMCPP__) || \
    defined(__PUREC__) || (defined(__TURBOC__) &&                    \
    (defined(__OS2__) || defined(__WIN32__)))
  eps = (REAL)200000.0 * MACH_EPS;
#else
  eps = (REAL)11000.0 * MACH_EPS;
#endif


  /* ------------------ print out input ---------------------------- */

  printf("\n"
         "fast complex Fourier transform (FFT)\n"
         "====================================\n\n"
         "tau = %3d    (Number N of nodes: N = 2^tau)\n"
         "k   = %3d    (Frequency factor for the function  "
         "exp(k * i * x))"
         "\n",
         tau, k
        );


  if (tau > 8 * sizeof(int) - 2)      /* tau too large?               */
  {                                   /* (overflow for 2^tau!)        */
    fehler_melden("tau zu gross", 0, __FILE__, __LINE__);
    return 3;
  }

  N = 1 << tau;                                         /* 2 ^ tau    */

  vmblock = vminit();                            /* allocate storage  */
  y    = (complex *)vmalloc(vmblock, VVEKTOR, N, sizeof(*y));
  rett = (complex *)vmalloc(vmblock, VVEKTOR, N, sizeof(*rett));
  if (! vmcomplete(vmblock))                     /* lack of memory    */
  {
    fehler_melden("lack of memory", 0, __FILE__, __LINE__);
    return 4;
  }

  faktor = TWO * k * PI / N;              /* generate N equidistant   */
  for (j = 0; j < N; j++)                 /* nodes for exp(k*i*x)     */
    rett[j].x = y[j].x = COS(faktor * j), /* in the interval [0,2pi)  */
    rett[j].y = y[j].y = SIN(faktor * j); /* for x = j*2*pi/N, j=..   */


  /* ------- compute coefficients of the trigonometric -------------- */
  /* ------- interpolating polynomial (Fourier analysis) ------------ */

#ifdef INTERAKTIV
  fprintf(stderr, "\n");
#endif
  fprintf(stderr, "complex Fourier analysis...");
  fehler = fft(tau, y, synthese = FALSE) ;
  fprintf(stderr, " done\n");

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("fft(): tau < 1", 10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("fft(): tau too large",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ---------------------- print results --------------------------- */

  printf("\n\n"
         "discrete Fourier analysis of  f(x) = exp(k * i * x)  "
         "in [0,2pi):\n"
         "-----------------------------------------------------"
         "-----------\n");
  for (j = 0; j < N; j++)
    if (FABS(y[j].x) + FABS(y[j].y) > eps)
      printf("non vanishing Fourier coefficient  y[%5d] = "
             "(%- 14.5"LZP"e,%- 14.5"LZP"e)\n", j, y[j].x, y[j].y);


  /* ------- reconstruct original function values from the ---------- */
  /* ------- Fourier coefficients (Fourier synthesis)  -------------- */

  fprintf(stderr, "complex Fourier synthesis...");
  fehler = fft(tau, y, synthese = TRUE);
  fprintf(stderr, " done\n");

  switch (fehler)
  {
    case 0:
      break;
    case 1:
      fehler_melden("fft(): tau < 1", 10 + fehler, __FILE__, __LINE__);
      break;
    case 2:
      fehler_melden("fft(): tau too large",
                    10 + fehler, __FILE__, __LINE__);
  }

  if (fehler != 0)
    return 10 + fehler;


  /* ----------- report results of synthesis ------------------------ */

  printf("\n\n"
         "discrete Fourier synthesis for  f(x) = exp(k * i * x)  "
         "in [0,2pi):\n"
         "-------------------------------------------------------"
         "-----------\n");
  for (j = 0; j < N; j++)
    if (FABS(y[j].x - rett[j].x) + FABS(y[j].y - rett[j].y) > eps)
      printf ("Deviation %12.3"LZP"e when synthesizing "
              "function value  y[%5d] = "
              "(%- 14.5"LZP"e,%- 14.5"LZP"e)\n",
              FABS(y[j].x - rett[j].x) + FABS(y[j].y - rett[j].y), j,
              y[j].x, y[j].y);


   synthese = synthese;                          /* calm the compiler */

   return 0;
}
