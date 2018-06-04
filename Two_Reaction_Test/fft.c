#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/*.BA*/



/*.FE{C 8.1.5.3}{Complex Discrete Fourier Transformation (FFT)}
                {Complex Discrete Fourier Transformation (FFT)}*/

/*.BE*/
/* -------------------------- MODULE fft.c -------------------------- */

/***********************************************************************
*                                                                      *
* functions used for the fast Fourier transformation:                  *
* ---------------------------------------------------                  *
* - real:                                          rfft()              *
* - complex:                                       fft()               *
* - complex or real with                                               *
*       arbitrary number of nodes or coefficients: fftb()              *
*                                                                      *
* Programing language: ANSI C                                          *
* Compiler:            Borland C++ 2.0                                 *
* Computer:            IBM PS/2 70 with 80387                          *
* Author:              Klaus Niederdrenk (FORTRAN)                     *
* Adaptation:          Juergen Dietel, Computer Center, RWTH Aachen    *
* Source:              [NIED84]                                        *
* Date:                9.28.1992                                       *
*                                                                      *
***********************************************************************/

#include <basis.h>      /*  for  PI, REAL, ONE, ZERO, TWO, HALF, SIN, */
                        /*       COS                                  */
#include <vmblock.h>    /*  for  vminit, vmalloc, MATRIX, vmfree      */
#include <u_proto.h>    /*  for  house                                */
#include <fft.h>        /*  for  rfft, fft, fftb, complex             */



/* ------------------------------------------------------------------ */
/*.BA*/

int rfft        /* fast real Fourier transform .......................*/
/*.IX{rfft}*/
        (
         int  tau,        /* 2^tau = number of nodes .................*/
         REAL y[],        /* nodes / Fourier coefficients ............*/
         int  synthese    /* direction of transform ..................*/
        )                 /* error code ..............................*/

/***********************************************************************
* If synthese = 0, this function computes the discrete Fourier         *
* coefficients                                                         *
*          a[0], ...., a[N/2]  and  b[1], ...., b[N/2 - 1]             *
* of the corresponding discrete partial Fourier sum                    *
*    a[0] + Sum (k=1,2,...,N/2-1)  (a[k] * cos(k * omega * x) +        *
*                                    b[k] * sin(k * omega * x))        *
*         + a[N/2]*cos(N/2*omega*x)                                    *
* for  N = 2^tau  given real functional values y[0], ...., y[N-1].     *
* Here  omega  =  2 * PI / L  (L = length of period).                  *
* If synthese = 1, we compute the inverse transform. (Fourier          *
* synthesis)                                                           *
* The (inverse) transform is executed via fast Fourier Transform, FFT, *
* for half of the total length.                                        *
.BE*)
* This follows the ideas in:                                           *
*      K. Niederdrenk: Die endliche Fourier- und Walsh-Transformation  *
*                      mit einer Einfuehrung in die Bildverarbeitung,  *
*                      2nd vol. 1984, Wiesbaden.                       *
* This book also contains a detailed derivation of the algorithm.      *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* tau:      log 2 of the length of the data, i.e., the number of nodes *
*           is N = 2^tau. (tau >= 1)                                   *
* y:        [0..N-1] vector, depending on synth :                      *
*             synthese = 0: y denotes the function values              *
*             synthese = 1: y denotes the discrete Fourier coefficients*
*                             y[0] = a[0],                             *
*                             y[k] = a[(k+1)/2] for k=1,3, ..,N-1      *
*                             y[k] = b[k/2]     for k=2,4, ..,N-2,     *
*                             i.e., y contains the values              *
*                             a[0], a[1], b[1], a[2], b[2], ...        *
* synthese: Flag to control the direction of the transform             *
*             synthese = 0: compute the discrete Fourier coefficients  *
*                           (Fourier analysis)                         *
*             synthese = 1: Compute the functional values              *
*                           (Fourier synthesis)                        *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* y: [0..N-1] vector, depending on   synthese :                        *
*      synthese = 0: y contains the discrete Fourier coefficients      *
*                      y[0] = a[0],                                    *
*                      y[2*k-1] = a[k] for k=1,2, ..,N/2               *
*                      y[2*k] = b[k]   for k=1,2, ..,N/2-1,            *
*                      i.e., y contains the string of values           *
*                      a[0], a[1], b[1], a[2], b[2], ...               *
*      synthese = 1: y contains the functional values                  *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: all is ok                                                         *
* 1: tau < 1                                                           *
* 2: tau too large (overflow for 2^tau)                                *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, PI, SIN, COS, ONE, ZERO, TWO, HALF                             *
.BA*)
***********************************************************************/
/*.BE*/

{

  int  N,           /* number of nodes                                */
       Nd2,         /* N / 2                                          */
       Nd4,         /* N / 4                                          */
       sigma,       /* bit reversal of the tau-1 binary bits of j     */
       min_n,       /* 2 ^ (tau - 1 - n)                              */
       n_min_0,     /* 2 ^ n                                          */
       n_min_1,     /* 2 ^ (n - 1)                                    */
       ind1,        /* aux indices                                    */
       ind2,        /*                                                */
       k, j,        /* Loop indices                                   */
       n, l;        /*                                                */
  REAL faktor,      /* Normalizing factor 2/N for Fourier analysis;   */
                    /* set to 1 for synthesis                         */
       arg,         /* argument of (wr,wi)                            */
       arg_m,       /* argument of the Nth roots of unity             */
       arg_md2,     /* same for (N/2)th roots of unity                */
       ew_r,        /* Real and imaginary parts of                    */
       ew_i,        /* the primitive nth root of unity                */
       eps_r,       /* Real and imaginary parts of the primitive      */
       eps_i,       /* (Nth root of unity) ^ k   or of the            */
                    /* ((N/2)th root of unity) ^ (l * 2^min_n)        */
       ur, ui,      /* aux variable                                   */
       wr, wi,      /* Real and imaginary parts of the primitive      */
                    /* ((N/2)th root of unity) ^ (2^min_n)            */
       rett,        /* aux variables                                  */
       yhilf,       /*                                                */
       hilf1,       /*                                                */
       hilf2,       /*                                                */
       hilf3,       /*                                                */
       hilf4;       /*                                                */


  if (tau < 1)                        /* tau too small?               */
    return 1;

  if (tau > 8 * sizeof(int) - 2)      /* tau too large ?              */
    return 2;                         /* (overflow for  2^tau!)       */

  N       = 1 << tau;                 /* N  =  2 ^ tau                */
  Nd2     = N / 2;
  Nd4     = Nd2 / 2;
  faktor  = ONE / Nd2;
  arg_md2 = TWO * PI * faktor;
  arg_m   = HALF * arg_md2;

  if (synthese)
    faktor = ONE;


  if (synthese)               /* combine real data for an FFT of      */
  {                           /* half the length                      */
    yhilf =  y[1];
    y[1]  =  y[0] - y[N - 1];
    y[0]  += y[N - 1];

    ew_r = COS(arg_m);       /* (ew_r,ew_i) = Nth root of unity       */
    ew_i = SIN(arg_m);
    eps_r = ONE;             /* (eps_r,eps_i) = (Nth root of unity)^k */
    eps_i = ZERO;

    for (k = 1; k < Nd4; k++)
    {
      ind1 = 2 * k;
      ind2 = N - ind1;
      rett = eps_r;
      eps_r = rett * ew_r - eps_i * ew_i;
      eps_i = rett * ew_i + eps_i * ew_r;
      hilf1 = HALF * (eps_r * (yhilf   - y[ind2 - 1]) +
                     eps_i * (y[ind1] + y[ind2]));
      hilf2 = HALF * (eps_i * (yhilf   - y[ind2 - 1]) -
                      eps_r * (y[ind1] + y[ind2]));
      hilf3 = HALF * (yhilf   + y[ind2 - 1]);
      hilf4 = HALF * (y[ind1] - y[ind2]);
      yhilf = y[ind1 + 1];
      y[ind1]     = hilf3 - hilf2;
      y[ind1 + 1] = hilf1 - hilf4;
      y[ind2]     = hilf2 + hilf3;
      y[ind2 + 1] = hilf1 + hilf4;
    }
    y[Nd2 + 1] = y[Nd2];
    y[Nd2]     = yhilf;
  }


  for (j = 0; j < Nd2; j++)                       /* use bit reversal,*/
  {                                               /* in case of       */
    for (k = j, n = 1, sigma = 0; n < tau; n++)   /* Fourier analysis */
      sigma <<= 1,                                /* normalize as well*/
      sigma |=  k & 1,
      k     >>= 1;
    if (j <= sigma)
      ind1 = 2 * j,
      ind2 = 2 * sigma,
      ur = y[ind1],
      ui = y[ind1 + 1],
      y[ind1]     = y[ind2]     * faktor,
      y[ind1 + 1] = y[ind2 + 1] * faktor,
      y[ind2]     = ur          * faktor,
      y[ind2 + 1] = ui          * faktor;
  }


  for (min_n = Nd2, n_min_1 = 1, n = 1; n < tau; n++)   /* perform FFT*/
  {                                                     /* of half the*/
    min_n   /= 2;                                       /* length     */
    n_min_0 =  2 * n_min_1;
    arg =  arg_md2 * min_n;
    wr  =  COS(arg);
    wi  =  synthese ? ONE : -ONE;
    wi  *= SIN(arg);
    eps_r = ONE;                         /* (eps_r,eps_i) =           */
    eps_i = ZERO;                        /* ((N/2)th root of unity) ^ */
                                         /* (l * 2^min_n)             */
    for (l = 0; l < n_min_1; l++)
    {
      for (j = 0; j <= Nd2 - n_min_0; j += n_min_0)
      {
        ind1 = (j + l) * 2;
        ind2 = ind1 + n_min_0;
        ur = y[ind2] * eps_r - y[ind2 + 1] * eps_i;
        ui = y[ind2] * eps_i + y[ind2 + 1] * eps_r;
        y[ind2]     =  y[ind1]     - ur;
        y[ind2 + 1] =  y[ind1 + 1] - ui;
        y[ind1]     += ur;
        y[ind1 + 1] += ui;
      }
      rett  = eps_r;
      eps_r = rett * wr - eps_i * wi;
      eps_i = rett * wi + eps_i * wr;
    }
    n_min_1 = n_min_0;
  }


  if (! synthese)                       /* separate the jointly trans-*/
  {                                     /* formed data in case of     */
    yhilf = y[N - 1];                   /* Fourier analysis           */
    y[N - 1] = HALF * (y[0] - y[1]);
    y[0]     = HALF * (y[0] + y[1]);

    ew_r =  COS(arg_m);
    ew_i = -SIN(arg_m);

    eps_r = ONE;                           /* (eps_r,eps_i) =         */
    eps_i = ZERO;                          /* (Nth root of unity) ^ k */

    for (k = 1; k < Nd4; k++)
    {
      rett = eps_r;
      eps_r = rett * ew_r - eps_i * ew_i;
      eps_i = rett * ew_i + eps_i * ew_r;
      ind1 = k * 2;
      ind2 = N - ind1;
      hilf1 = HALF * (eps_i * (y[ind1]     - y[ind2]) +
                      eps_r * (y[ind1 + 1] + yhilf));
      hilf2 = HALF * (eps_r * (y[ind1]     - y[ind2]) -
                      eps_i * (y[ind1 + 1] + yhilf));
      hilf3 = HALF * (y[ind1]     + y[ind2]);
      hilf4 = HALF * (y[ind1 + 1] - yhilf);
      yhilf = y[ind2 - 1];
      y[ind1 - 1] = hilf1 + hilf3;
      y[ind1]     = hilf2 - hilf4;
      y[ind2 - 1] = hilf3 - hilf1;
      y[ind2]     = hilf2 + hilf4;
    }
    y[Nd2 - 1] = y[Nd2];
    y[Nd2]     = yhilf;
  }


  return 0;
}



/* --------------- complex multiplication:  c = a * b --------------- */

#define COMMUL(c, a, b)         \
/*.IX{COMMUL}*/                 \
{                               \
  REAL re, im;                  \
  re  = a.x * b.x - a.y * b.y;  \
  im  = a.x * b.y + a.y * b.x;  \
  c.x = re;                     \
  c.y = im;                     \
}



/* ------------------------------------------------------------------ */
/*.BA*/

int fft         /* fast complex Fourier transform ....................*/
/*.IX{fft}*/
       (
        int     tau,      /* 2^tau = number of nodes .................*/
        complex y[],      /* node vector or Fourier coefficients .....*/
        int     synthese  /* direction of transform ..................*/
       )                  /* error code ..............................*/

/***********************************************************************
* If synthese = 0, this function computes the discrete Fourier         *
* coefficients  c(-N/2),....,c(N/2-1)  of the corresponding discrete   *
* partial Fourier sum                                                  *
*        Sum (k=-N/2,...,N/2-1) (c(k) * exp(i * k * omega * x))        *
* for  N = 2^tau  given real or complex functional values y[0], ....,  *
* y[N-1].  Here  omega  =  2 * PI / L  (L = length of period).         *
* If synthese = 1, we compute the inverse transform. (Fourier          *
* synthesis)                                                           *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* tau:      log 2 of the number N of data points; N = 2^tau for tau > 0*
* y:        [0..N-1,0..1] vector for N complex numbers, dependng on    *
*           the value of  synthese :                                   *
*             synthese = 0: y contains the function values y(i).       *
*             synthese = 1: y contains the discrete Fourier            *
*                           coefficients c(i), namely those for        *
*                           i=0,...,N/2-1 in y[i] and those for        *
*                           i=-N/2,...,-1 in y[i+N].                   *
* synthese: Flag governing the direction of the transform :            *
*             synthese = 0: compute the discrete Fourier coefficients  *
*                           (Fourieranalysis)                          *
*             synthese = 1: compute functional values                  *
*                           (Fouriersynthesis)                         *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y: [0..N-1] vector for N complex numbers, depending on synthese :    *
*      synthese = 0: the discrete  Fourier coefficients c(i), namely   *
*                    those with indices i=0,...,N/2-1 in y[i] and      *
*                    those for i=-N/2,...,-1 in y[i+N];                *
*      synthese = 1: the function values  y(i)                         *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: all is ok                                                         *
* 1: tau < 1                                                           *
* 2: tau too large (overflow in 2^tau)                                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* complex, REAL, PI, SIN, COS, ONE, ZERO, TWO, COMMUL                  *
.BA*)
***********************************************************************/
/*.BE*/

{
  int     N,        /* number of data points                          */
          n_min_0,  /* 2 ^ n                                          */
          n_min_1,  /* 2 ^ (n - 1)                                    */
          sigma,    /* bit reversal of the tau bits of j              */
          ind1,     /* aux variables                                  */
          ind2,     /*                                                */
          j, n, l;  /* Loop indices                                   */
  REAL    faktor,   /* normalizing factor 1/N for Fourier coefficients*/
                    /* in case of F. analysis; set equal to 1 in case */
                    /* of synthesis of the function values            */
          ewphi;    /* argument of the Nth root of unity              */
  complex ew,       /* Nth root of unity                              */
          w,        /* ew ^ (2^(tau - n))                             */
          eps,      /* ew ^ (l * 2^(tau - n))                         */
          u;        /* aux complex variable                           */


  if (tau < 1)                        /* tau too small ?              */
    return 1;

  if (tau > 8 * sizeof(int) - 2)      /* tau too large ?              */
    return 2;                         /* (overflow for  2^tau!)       */


  N      = 1 << tau;                  /* N  =  2 ^ tau                */
  faktor = ONE / N;
  ewphi  = -TWO * PI * faktor;
  if (synthese)
    ewphi = -ewphi;
  ew.x = COS(ewphi);
  ew.y = SIN(ewphi);

  if (synthese)
    faktor = ONE;


  for (j = 0; j < N; j++)                         /* Bit reversal     */
  {                                               /* and normalize in */
    for (l = j, sigma = 0, n = tau; n > 0; n--)   /* Fourier analysis */
      sigma <<= 1,                                /* case             */
      sigma |=  l & 1,
      l     >>= 1;
    if (j <= sigma)
      u          = y[j],
      y[j].x     = y[sigma].x * faktor,
      y[j].y     = y[sigma].y * faktor,
      y[sigma].x = u.x        * faktor,
      y[sigma].y = u.y        * faktor;
  }


  /* ------ perform the transformation (Analysis or Synthesis) ------ */

  for (n_min_1 = 1, n = 1; n <= tau; n++, n_min_1 = n_min_0)
  {
    for (w = ew, l = tau - n; l != 0; l--)    /* w = ew ^ (2^(tau-n)) */
      COMMUL(w, w, w);

    eps.x = ONE;                          /* initialize eps = 1       */
    eps.y = ZERO;
    for (n_min_0 = n_min_1 + n_min_1, l = 0; l < n_min_1; l++)
    {
      for (j = 0; j <= N - n_min_0; j += n_min_0)
      {
        ind1      =  j + l;
        ind2      =  ind1 + n_min_1;
        COMMUL(u, y[ind2], eps);         /* u        =  y(ind2) * eps */
        y[ind2].x =  y[ind1].x - u.x;    /* y(ind2)  =  y(ind1) - u   */
        y[ind2].y =  y[ind1].y - u.y;
        y[ind1].x += u.x;                /* y(ind1)  += u             */
        y[ind1].y += u.y;
      }
      COMMUL(eps, eps, w);               /* eps      *= w             */
    }
  }


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int fftb        /* complex FFT for arbitary many nodes ...............*/
/*.IX{fftb}*/
        (
         int     N,        /* number of nodes ........................*/
         complex y[],      /* node vector or Fourier coefficients ....*/
         int     synthese  /* direction for transformation ...........*/
        )                  /* error code .............................*/

/***********************************************************************
* For  synthese = 0  and an arbitrary number N of real or complex      *
* funktional values y(0),...,y(N-1), this function determines the      *
* discrete Fourier coefficients c(k) of the discrete partial  Fourier  *
* sum                                                                  *
*     Sum (k=-N/2,...,N/2-1) (c(k) * exp(i * k * omega * x)),          *
* if  N is even, or                                                    *
*     Sum (k=-(N-1)/2,...,(N-1)/2) (c(k) * exp(i * k * omega * x)),    *
* if  N is odd; here :                                                 *
*         omega  =  2 * PI / L            (L = Period length),         *
*         i      =  imaginary unit  (0,1).                             *
* For  synthese = 1  we perform the inverse transformation (Fourier    *
* synthesis).                                                          *
* The function uses the FFT for powers of 2 via discrete convolutions. *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* N         number of nodes or Fourier coefficients                    *
*           Note that both of the aux vectors f1 and g below need to   *
*           accomodate  2^tau complex numbers where N <= (2^tau + 1)/3.*
* y:        [0..N-1] vector for N complex numbers depending on the     *
*           value of synthese :                                        *
*             synthese = 0: y contains the function values  y(k)       *
*             synthese = 1: y contains the discrete Fourier            *
*                           coefficients c(k), namely c(k) for         *
*                           k=0,...,N/2-1 (N even) or c(k) for         *
*                           k=0,...,(N-1)/2 (N odd) in y[k] and c(k)   *
*                           for k=-N/2,...,-1 (N even) or c(k)         *
*                           for k=-(N-1)/2,...,-1 (N odd) in  y[k+N].  *
* synthese: Flag to control direction of transformation:               *
*             synthese = 0: compute discrete Fourier coefficients      *
*                           (Fourieranalysis)                          *
*             synthese = 1: compute function values                    *
*                           (Fouriersynthesis)                         *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* y: [0..N-1] vector for N complex numbers depending on synthese :     *
*             synthese = 1: y contains the function values  y(k)       *
*             synthese = 0: y contains the discrete Fourier            *
*                           coefficients c(k), namely c(k) for         *
*                           k=0,...,N/2-1 (N even) or c(k) for         *
*                           k=0,...,(N-1)/2 (N odd) in y[k] and c(k)   *
*                           for k=-N/2,...,-1 (N even) or c(k)         *
*                           for k=-(N-1)/2,...,-1 (N odd) in  y[k+N].  *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: all is ok                                                         *
* 1: N < 1                                                             *
* 2: N too larges (overflow in tau = 2^N)                              *
* 3: lack of memory                                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* complex, fft, REAL, PI, SIN, COS, ONE, ZERO, TWO, vminit, vmalloc,   *
* MATRIX, vmcomplete, vmfree, COMMUL                                   *
.BA*)
***********************************************************************/
/*.BE*/

{
  int     tau,       /* [ln(3*N-2)/ln(2)+1]                           */
          l,         /* 2 ^ tau                                       */
          k,         /* Loop variable                                 */
          fehler;    /* error code von fft()                          */
  complex ew1,       /* Nth root of unity                             */
          ew2,       /* ew2  =  ew1 ^ 2                               */
          ew3,       /* aux variables                                 */
          ew4,       /*                                               */
          ewk,       /* ewk  =  ew1 ^ (k^2),                          */
          *f1,       /* [0..2^tau-1] aux vectors                      */
          *g;        /*                                               */
  REAL    faktor,    /* nornmalizing factor  1/N for the Fourier      */
                     /* coefficients in case of Fourier analysis;     */
                     /* otherwise we set factor  =  1                 */
          faktl,     /* faktor * l                                    */
          ew1phi;    /* argument of the complex number ew1            */
  void    *vmblock;  /* List of dynamically allocated vectors         */


  if (N < 2)                            /* N too small ?              */
    return 1;


  /* ------------ find suitable power of 2 for the length ----------- */
  /* ------------ of the auxiliary vectors f1 and g ----------------- */

  tau = (int)(LOG(THREE * (REAL)N - TWO) / LOG(TWO)) + 1;

  if (tau > 8 * sizeof(int) - 2)      /* tau too large ?              */
    return 2;                         /* (overflow for 2^tau!)        */

  l = 1 << tau;                                 /* l  =  2 ^ tau      */
  if (l / 2 >= 3 * N - 2)                       /* l too large ?      */
    l /= 2,                                     /* halve l and reduce */
    tau--;                                      /* tau by 1           */


  /* -------------- allocate aux vectors ---------------------------- */

  vmblock = vminit();
  f1 = (complex *)vmalloc(vmblock, VVEKTOR, l, sizeof(*f1));
  g  = (complex *)vmalloc(vmblock, VVEKTOR, l, sizeof(*g));
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 3;
  }


  faktor = ONE / N;
  ew1phi = -PI * faktor;
  if (synthese)
    ew1phi = -ew1phi;

  ew1.x = COS(ew1phi);                     /* ew1  =  exp(i * ew1phi) */
  ew1.y = SIN(ew1phi);

  COMMUL(ew2, ew1, ew1);

  if (synthese)
    faktor = ONE;


  /* ------------- initialize aux vectors f1 and g ------------------ */

  for (k = 0; k < l; k++)
    f1[k].x = f1[k].y = g[k].x = g[k].y = ZERO;
  f1[0]      = y[0];
  g[N - 1].x = ONE;
  g[N - 1].y = ZERO;


  ewk = ew1;
  ew3 = ew1;

  for (k = 1; k < N; k++)
  {
    COMMUL(f1[k], y[k], ewk);              /* f1(k)     =  y(k) * ewk */
    g[N - 1 + k].x = ewk.x;                /* g(N-1+k)  =  CONJG(ewk) */
    g[N - 1 + k].y = -ewk.y;
    g[N - 1 - k]   = g[N - 1 + k];
    COMMUL(ew3, ew3, ew2);                 /* ew3       *= ew2        */
    COMMUL(ewk, ewk, ew3);                 /* ewk       *= ew3        */
  }


  /* ----------- dicrete convolution of the vectors f1 and g -------- */
  /* ----------- using the FFT                            ----------- */

  fehler = fft(tau, f1, FALSE);
  if (fehler)                          /* this should not occur !     */
  {
    vmfree(vmblock);
    return fehler;
  }

  fehler = fft(tau, g,  FALSE);
  if (fehler)                          /* ditto !                     */
  {
    vmfree(vmblock);
    return fehler;
  }

  for (k = 0; k < l; k++)
    COMMUL(f1[k], f1[k], g[k]);                    /* f1(k)  *=  g(k) */

  fehler = fft(tau, f1, TRUE);
  if (fehler)                          /* impossible again !          */
  {
    vmfree(vmblock);
    return fehler;
  }


  /* ------- store necessary data in the vector y ------------------- */

  faktl  = faktor * l;

  y[0].x = f1[N - 1].x * faktl;     /* y(0)  =  f1(N-1) * faktor * l  */
  y[0].y = f1[N - 1].y * faktl;
  ewk    = ew1;
  ew3    = ew1;

  for (k = 1; k < N; k++)
  {
    ew4.x = ewk.x * faktl;          /* y(k) = f1(k+N-1) * ewk * faktl */
    ew4.y = ewk.y * faktl;
    COMMUL(y[k], f1[k + N - 1], ew4);

    COMMUL(ew3, ew3, ew2);                            /* ew3  *=  ew2 */
    COMMUL(ewk, ewk, ew3);                            /* ewk  *=  ew3 */
  }


  vmfree(vmblock);
  return 0;
}

/* ------------------------------------------------------------------ */

static int ist_zweierpotenz    /* Is m a power of 2? .................*/
    (
     int m
    )

/***********************************************************************
* Find out whether m is a power of two; return t with 2^t = m if so,   *
* otherwise return 0.                                                  *
***********************************************************************/

{
  int potenz;
  int exponent;

  for (potenz = 1, exponent = 0; potenz != 0; potenz <<= 1, exponent++)
    if (potenz == m)
      return exponent;

  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int fdicht /* trigonometric interpolating polynomial at shifted nodes */
    (
     int     M,            /* number of nodes          ...............*/
     complex F[],          /* table of values for F       ............*/
     REAL    p,            /* Period of F   ..........................*/
     REAL    theta         /* shift                   ................*/
    )

/***********************************************************************
*  This function  computes the values of the trigonometric             *
*  interpolating polynomial (i.e. the discrete partial Fourier         *
*  sum) using the Fast Fourier Transform (FFT) for a given set         *
*  of function values F[0]..F[M-1] of a p-periodic                     *
*  function for equidistant nodes t[0]..t[M-1],                        *
*  t[j] = j*p/M  at a set of shifted nodes t[j] + theta  for           *
*  j=0..M-1 ("Increase of number of nodes").                           *
.BE*)
*                                                                      *
*  Input parameters:                                                   *
*  =================                                                   *
*  M      Number of original nodes.                                    *
*         If M is a power of two (M = 2^tau  for a positive integer    *
*         tau), we make use of the function fft(), which runs for all  *
*         values of M. If M is not a power of two, then we call on     *
*         fftb().                                                      *
*  F      [0..M-1] vector for M complex functional values              *
*  p      Period interval for the underlying function                  *
*  theta  Shift parameter: the newly computed values correspond to     *
*         arguments shifted by theta.                                  *
*                                                                      *
*  Output parameters:                                                  *
*  ==================                                                  *
*  F      [0..M-1] vector for M complex function values of the trig.   *
*         interpolating polynomial at the shifted equidistant nodes    *
*                                                                      *
*  Return value:                                                       *
*  =============                                                       *
*  0: no error                                                         *
*  1: M < 2                                                            *
*  2: tau too large for fft() or fftb(), i.e., M is too large          *
*  3: memory exceeded                                                  *
*                                                                      *
*  Global names used:                                                  *
*  ==================                                                  *
*  complex, REAL, ist_zweierpotenz, fft, fftb, TWO, PI, COS, SIN,      *
*  COMMUL                                                              *
*                                                                      *
************************************************************************
*                                                                      *
*  Author:             Klaus Niederdrenk (FORTRAN 77)                  *
*  Date:               06.30.1994                                      *
*  Adaptation:         J. Dietel                                       *
*  Date:               8.2.1995                                        *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/

{
  int     tau;                    /* log two ofM          */
  int     k;                      /* variable                       */
  int     fehler;                 /* eror code from fft() or fftb() */
  REAL    faktor;
  REAL    tetpi;
  complex ek;
  complex ekk;                    /* complex conjugate of `ek'      */
  complex h;


  if (M < 2)                            /* M too small ?    */
    return 1;

  tau = ist_zweierpotenz(M);    /* Is M a power of two? */

  /* Compute  discrete Fourier transform                   */

  if (tau != 0)                          /* M =  2^tau? */
  {
    fehler = fft(tau, F, 0);
    if (fehler)
      return fehler;
  }
  else                                   /* M not a power of two?     */
  {
    fehler = fftb(M, F, 0);
    if (fehler)
      return fehler;
  }

  /* Find values of the discrete Fourier transform at shifted nodes   */

  faktor = TWO * PI * theta / p;
  tetpi  = (REAL)(-M / 2) * faktor;
  ek.x   = COS(tetpi);
  ek.y   = SIN(tetpi);
  COMMUL(h, F[M / 2], ek);
  for (k = -(M - 1) / 2; k <= -1; k++)
  {
    tetpi = (REAL)k * faktor;
    ek.x  = COS(tetpi);
    ek.y  = SIN(tetpi);
    ekk   = ek;
    ekk.y = -ekk.y;
    COMMUL(F[k + M], F[k + M], ek);
    COMMUL(F[-k],    F[-k],    ekk);
  }
  F[M / 2] = h;

  /* Find functional values                                           */

  if (tau)                               /* M = 2^tau?? */
    return fft(tau, F, 1);
  else                                   /* M not a power of two?     */
    return fftb(M, F, 1);
}



/* ------------------------------------------------------------------ */
/*.BA*/

int fourn     /* Fourier transform for nonperiodic functions          */
    (
     int     M,            /* number of nodes          ...............*/
     complex F[],          /* table of values for f       ............*/
     REAL    a,            /* left endpoint of the interval of support*/
     REAL    deltax        /* step size            ...................*/
    )                      /* erorr code .............................*/

/***********************************************************************
*  This function uses the fast Fourier transform (FFT) to compute      *
*  approximate values of the Fourier transform                         *
*                                                                      *
*       f^(tj) = (Integral of) f(x) * exp(-i*tj*x) dx                  *
*                                                                      *
*  (i: imaginary unit; i^2 = -1 ) for j=-M/2..M/2-1.                   *
*  Here the nonperiodic function  f is known by its functional         *
*  values F[0]..F[M-1] at equidistant nodes with                       *
*  uniform distance  deltax. These nodes lie in the support of f.      *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
*  M       Number of nodes; M must be even.                            *
*          If M is a power of two (M = 2^tau  for a positive integer   *
*          tau), we make use of the function fft(), which runs for all *
*          values of M. If M is not a power of two, then we call on    *
*          fftb().                                                     *
*  F       [0..M-1] vector for the M complex functional values of f at *
*          the nodes                                                   *
*                  x[j] = a + j * deltax, j=0..M-1.                    *
*  a       starting point of the functional values.                    *
*  deltax  uniform distance of the nodes.                              *
*                                                                      *
*  Output parameters:                                                  *
*  ==================                                                  *
*  F       [0..M-1] vector with M complex funktion values for the      *
*          Fourier transform f^, as follows:                           *
*                  F[k+M]  contains  f^(tk)  for  k = -M/2..-1,        *
*                  F[k]    contains  f^(tk)  for  k = 0..-M/2-1,       *
*          where  tk = k / (M * deltax).                               *
*                                                                      *
*  Return value:                                                       *
*  =============                                                       *
*  0: no error                                                         *
*  1: M < 2                                                            *
*  2: tau too large in fft() or fftb(), i.e. M is too large            *
*  3: not enough memory                                                *
*                                                                      *
*  Global names used:                                                  *
*  ==================                                                  *
*  complex, REAL, ist_zweierpotenz, fft, fftb, TWO, COS, SIN,          *
*  COMMUL                                                              *
*                                                                      *
************************************************************************
*                                                                      *
*  Author:             Klaus Niederdrenk (FORTRAN 77)                  *
*  Date:               06.30.1994                                      *
*  Adaptation:         J. Dietel                                       *
*  Date:               8.3.1995                                        *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/

{
  int     tau;                    /* log two of M               */
  int     k;                      /* variable                   */
  int     fehler;                 /* error code of fft() or fftb() */
  REAL    x;
  REAL    faktor;
  REAL    arg;
  complex ek;
  complex ekk;                    /* conjugate of `ek'          */


  if (M < 2)                                /* M too small ?    */
    return 1;

  tau = ist_zweierpotenz(M);    /* decide whether M is a power of two */

  /* Determine discrete Fourier coefficients                    */

  if (tau != 0)                                  /* M =  2^tau? */
  {
    fehler = fft(tau, F, 0);
    if (fehler)
      return fehler;
  }
  else                                 /* M not a power of two? */
  {
    fehler = fftb(M, F, 0);
    if (fehler)
      return fehler;
  }

  /*      Transform to the nonperiodic case:                         */
  /*      Adjust values to those of the Fourier transform            */

  x      = (REAL)M * deltax;
  faktor = TWO * PI * a / x;
   for (k = -M / 2 + 1; k < 0; k++)
  {
    arg   = (REAL)k * faktor;
    ek.x  = x * COS(arg);
    ek.y  = x * SIN(arg);
    ekk   = ek;
    ekk.y = -ekk.y;
    COMMUL(F[k + M], F[k + M], ekk);
    COMMUL(F[-k],    F[-k],    ek);
  }
  F[0].x *= x;
  F[0].y *= x;
  arg    =  (REAL)(M / 2) * faktor;
  ek.x   =  x * COS(arg);
  ek.y   =  x * SIN(arg);
  COMMUL(F[M / 2], F[M / 2], ek);


  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int ffako     /* discrete cyclic periodic convolution and correlation */
    (
     int     M,            /* number of nodes          ...............*/
     complex F[],          /* table of function values for f .........*/
     complex H[]           /* table of function values for h .........*/
    )                      /* error code .............................*/

/***********************************************************************
*  This function uses the fast Fourier transform (FFT) to compute      *
*  the discrete values of the convolution                              *
*                                                                      *
*       Falt[j] = 1/M * ((sum: k=0..M-1) F[j-k] * H[k])                *
*                                                                      *
*  and of the discrete cyclic correlation                              *
*                                                                      *
*       Korr[j] = 1/M * ((sum: k=0..M-1) F[j+k] * conjug(H[k]))        *
*                                                                      *
*  of f and h at the given complex functional values                   *
*  F[0]..F[M-1]  und  H[0]..H[M-1]  at                                 *
*  equidistant nodes in the period interval                            *
*  for  j=0..M-1.                                                      *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
*  M  Number of nodes.                                                 *
*     If M is a power of two (M = 2^tau  for a positive integer tau),  *
*     we make use of the function fft(), which runs for all values of  *
*     M. If M  is not a power of two, then we call on fftb().          *
*  F  [0..M-1] vector with the complex function values for f.          *
*  H  ditto for h.                                                     *
*                                                                      *
*  Output parameters:                                                  *
*  ==================                                                  *
*  F  [0..M-1] vector with the complex values of the discrete cyclic   *
*     convolution                                                      *
*  H  [0..M-1] vector with the complex values of the discrete cyclic   *
*     correlation                                                      *
*                                                                      *
*  Return value:                                                       *
*  =============                                                       *
*  0: no error                                                         *
*  1: M < 2                                                            *
*  2: tau too large in fft() or fftb(), i.e., M too large              *
*  3: out of memory                                                    *
*                                                                      *
*  Global names used:                                                  *
*  ==================                                                  *
*  complex, ist_zweierpotenz, fft, fftb, COMMUL                        *
*                                                                      *
************************************************************************
*                                                                      *
*  Author:             Klaus Niederdrenk (FORTRAN 77)                  *
*  Date:               9.8.1994                                        *
*  Adaptation:         Juergen Dietel, Rechenzentrum, RWTH Aachen      *
*  Date:               8.2.1995                                        *
.BA*)
***********************************************************************/
/*.BE*/

{
  int     tau;                    /* log two of  M               */
  int     k;                      /* variable                    */
  int     fehler;                 /* error code  fft() or fftb() */
  complex h1;
  complex hkk;                    /* complex conjugate of`H[k]'  */


  if (M < 2)                                 /* M too small ?    */
    return 1;

  tau = ist_zweierpotenz(M);              /* M a power of two ?  */

  /* Determine discrete Fourier transform                        */

  if (tau != 0)                                    /* M = 2^tau? */
  {
    fehler = fft(tau, F, 0)  |
             fft(tau, H, 0);
    if (fehler)
      return fehler;
  }
  else                                 /* M not a powewr of two? */
  {
    fehler = fftb(M, F, 0)  |
             fftb(M, H, 0);
    if (fehler)
      return fehler;
  }

  /* Determine the Fourier transform of the discrete cyclic      */
  /* convolution and of the discrete cyclic correlation          */

  for ( k = 0; k < M; k++)
  {
    COMMUL(h1,   F[k], H[k]);
    hkk   = H[k];

    hkk.y = -hkk.y;
    COMMUL(H[k], F[k], hkk);
    F[k]  = h1;
  }

  /* Determine functional values                                 */

  if (tau)                                         /* M = 2^tau? */
    return fft(tau, F, 1)  |
           fft(tau, H, 1);
  else                                  /* M not a power of two ?*/
    return fftb(M, F, 1)   |
           fftb(M, H, 1);
}



/* ------------------------------------------------------------------ */
/*.BA*/

int ffakon /* discrete cyclic nonperiodic convolution and correlation */
    (
     int     M,            /* number of nodes for f minus one         */
     complex F[],          /* table of values of        f ............*/
     int     N,            /* number of nodes for h minus one  .......*/
     complex H[],          /* table of values for h       ............*/
     int     tau,          /* L = 2^tau = Length of  F and H .........*/
     REAL    deltax,       /* step size            ...................*/
     int     faltung       /* flag selecting convolution / correlation*/
    )                      /* error code .............................*/

/***********************************************************************
*  This function uses the fast Fourier transform (FFT) to compute      *
*  approximations for the convolution                                  *
*                                                                      *
*       Falt[j] = (Integral of) f(xj-t) * h(t) dt                      *
*                                                                      *
*  or approximations for the correlation                               *
*                                                                      *
*       Korr[j] = (Integral of) f(xj+t) * conjug(h(t)) dt              *
*                                                                      *
*  for  j=0..M+N.                                                      *
*  Here f and h are two nonperiodic functions given by their           *
*  funktional values  F[0]..F[M]  und  H[0]..H[N]                      *
*  at equidistant nodes with uniform distance  deltax.                 *
*  These nodes lie in the support of f and h.                          *
.BE*)
*                                                                      *
*  Input parameters:                                                   *
*  =================                                                   *
*  M        M+1  functional values of f                                *
*  F        [0..L-1] vector with M+1 complex function values for f     *
*           (L = 2^tau)                                                *
*  N        N+1  functional values of h                                *
*  H        [0..L-1] vector with N+1 complex function values for h     *
*           (L = 2^tau)                                                *
*  tau      log two of the length of F and H;                          *
*           tau >= 1  and  L = 2^tau >= M+N+1.                         *
*  deltax   step size                                                  *
*  faltung  flag: = 0: correlation; ~= 0: convolution                  *
*                                                                      *
*  Output parameters:                                                  *
*  ==================                                                  *
*  F        [0..L-1] vector with M+N+1 complex approximate values      *
*           of the convolution Falt[j] (flag `faltung' set) or of      *
*           correlation Korr[j] (Flag `faltung' cleared) for j=0..M+N  *
*  H        [0..L-1] vector of discrete Fourier coefficients of h      *
*           (auxiliary vector)                                         *
*                                                                      *
*  Return value:                                                       *
*  =============                                                       *
*  0: no error                                                         *
*  1: M < 1  or  N < 1                                                 *
*  2: tau < 1  or  tau too large                                       *
*                                                                      *
*  Global names used:                                                  *
*  ==================                                                  *
*  complex, REAL, fft, ZERO, COMMUL                                    *
*                                                                      *
************************************************************************
*                                                                      *
*  Author:             Klaus Niederdrenk (FORTRAN 77)                  *
*  Date:               9.8.1994                                        *
*  Adaptation:         Juergen Dietel, Rechenzentrum, RWTH Aachen      *
*  Date:               8.3 1995                                        *
*                                                                      *
.BA*)
***********************************************************************/
/*.BE*/

{
  int     L;                                  /* 2 ^ tau              */
  int     j;                                  /* variable             */
  int     k;                                  /* variable             */
  int     fehler;                             /* error code of  fft() */
  complex Hkk;
  REAL    x;


  /* Check input                                                      */

  if (M < 1 || N < 1)                            /* M or N too small? */
    return 1;

  if (tau < 1)                                 /* tau too small?      */
    return 2;

  if (tau > 8 * sizeof(int) - 2)               /* tau too large?      */
    return 2;                              /* (overflow for 2^tau!)   */

  L = 1 << tau;                                        /* L  =  2^tau */

  /* Initialize expanded vectors                                      */

  if (faltung)
    for (j = M + 1; j < L; j++)
      F[j].x = F[j].y = ZERO;
  else
  {
    for (j = M + N + 1; j < L; j++)
      F[j].x = F[j].y = ZERO;
    for (j = M + N; j >= N; j--)
      F[j] = F[j - N];
    for (j = 0; j < N; j++)
      F[j].x = F[j].y = ZERO;
  }
  for (j = N + 1; j < L; j++)
    H[j].x = H[j].y = ZERO;

  /* Determine the  discrete Fourier coefficients via FFT for         */
  /* powers of two                                                    */

  fehler = fft(tau, F, 0)  |
           fft(tau, H, 0);
  if (fehler)
    return fehler;

  if (faltung)
    for (k = 0; k < L; k++)
    {
      COMMUL(F[k], F[k], H[k]);
    }
  else
    for (k = 0; k < L; k++)
    {
      Hkk   = H[k];
      Hkk.y = -Hkk.y;
      COMMUL(F[k], F[k], Hkk);
    }

  /* Determine function values                                        */

  fehler = fft(tau, F, 1);
  if (fehler)
    return fehler;

  x = (REAL)L * deltax;
  for (j = 0; j <= M + N; j++)
    F[j].x *= x,
    F[j].y *= x;
  for (j = M + N + 1; j < L; j++)
    F[j].x = F[j].y = ZERO;


  return 0;
}
/*.BA*/



/*.FE{C 8.2.2}
     {Nonlinear Root-Mean-Square Fitting}
     {Nonlinear Root-Mean-Square Fitting}*/

/*.BE*/
/* ------------------------------------------------------------------ */

static REAL quadsum
/*.IX{quadsum}*/
                   (
                    int  n,
                    REAL v[]
                   )

/***********************************************************************
* Compute the euclidean inner product of the elements in the vector    *
* v and return its value.                                              *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO, sqr                                                      *
***********************************************************************/

{
  REAL quadratsumme;
  int  i;

  for (quadratsumme = ZERO, i = n; i != 0; i--, v++)
    quadratsumme += sqr(*v);

  return quadratsumme;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int nli_happr    /* non-linear least squares .........................*/
/*.IX{nli\unt happr}*/
             (
              int       m,         /* number of nodes ................*/
              int       n,         /* number of coefficients - 1 .....*/
              REAL      x[],       /* x-values .......................*/
              REAL      y[],       /* y-values .......................*/
              REAL      w[],       /* positive weights ...............*/
              approxfnk PHI,       /* approximating function .........*/
              int       ablOK,     /* flag for derivatives via ABL ...*/
              ableitfnk ABL,       /* partial derivative of PHI ......*/
              int       *maxIt,    /* max/current number of iterations*/
              REAL      RelEps,    /* relative error bound ...........*/
              REAL      c[],       /* Starting/optimal coefficient ...*/
              REAL      *MiQuFe    /* mean square error ..............*/
             )                     /* error code .....................*/

/***********************************************************************
* Compute the coefficients c for the generally nonlinear function      *
*                       PHI(c[0],...,c[n],x)                           *
* so that                                                              *
*   (*)   (y[0] - PHI(c,x[0]))^2 + ... + (y[m-1] - PHI(c,x[m-1]))^2    *
* is minimal, i.e., so that the graph of PHI approximates the m nodes  *
* (x[i],y[i]), i=0, ..., m-1 with m >= n optimally wrt. the mean square*
* error.                                                               *
*                                                                      *
* From a given starting solution we use the damped Newton method for   *
* non-linear systems of equations to compute the optimal parameters    *
* c[k] of the least squares approximation  PHI.                        *
* Here the linear minimization problem that arises in each iteration   *
* step is solved using a  Householder transform.                       *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* m:      number of nodes                                              *
* n:      number of coefficients of the model functions - 1            *
* x:      [0..m-1] x-values                                            *
* y:      [0..m-1] y-values                                            *
* w:      [0..m-1] weight vector                                       *
* PHI:    points to a function which supplies the value of the approx. *
*         function with coefficients c[0], ..., c[n] at x.             *
*         PHI has the form :                                           *
*               REAL PHI(REAL c[], REAL x)                             *
*               {                                                      *
*                 return <value of function at  x>                     *
*               }                                                      *
* ablOK:  Flag designating whether the partial derivatives of the model*
*         function wrt. its coefficients can be computed via the       *
*         function  ABL() (TRUE), or whether they must be computed via *
*         central differencs quotients (FALSE).                        *
* ABL:    Pointer for a function which computes the partial derivatives*
*         f[0],...,f[n] of the model function wrt. its coefficients    *
*         c[0],...,c[n] at a. This function must have the form:        *
*               void ABL(REAL x, REAL c[], REAL f[])                   *
*               {                                                      *
*                 f[0] = <derivative of function wrt. c[0]>            *
*                 ...                                                  *
*                 f[n] = <derivative of function wrt. c[n]>            *
*               }                                                      *
*         If the formal derivatives are unknown, set the flag ablOK to *
*         zero (= FALSE): In this case the partial derivatives are     *
*         approximated by central differenz quotients.                 *
* maxIt:  maximal number of iterations                                 *
* RelEps: relative error bound for optimal coefficients                *
* c:      [0..n] vector with starting coefficients                     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* c:      [0..n] optimal coefficient vector                            *
* maxIt   actual number of iterations                                  *
* MiQuFe: mean square error                                            *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* 0: no error                                                          *
* 1: m, maxIt or RelEps too small                                      *
* 2: error in Householder transformation                               *
* 3: lack of storage                                                   *
* 4: after maxIt iterations the desired accuracy has not been reached; *
*    i.e.,  RelEps is too small or the iterations converge too slowly  *
*    due to an inappropriate start.                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* quadsum, REAL, approxfnk, ableitfnk, MACH_EPS, FABS, SQRT, POW, ONE, *
* vminit, vmalloc, vmcomplete, vmfree, VEKTOR, MATRIX, house, ZERO,    *
* TWO, HALF, THREE                                                     *
.BA*)
***********************************************************************/
/*.BE*/

{
  void *vmblock;      /* List of dyn. allocated vectors and matrices  */
  int  Iteration,     /* iteration counter for Newton method          */
       i, k;          /* Loop indices                                 */
  REAL **Ab,          /* [0..m-1,0..n] system matrix                  */
       *b,            /* [0..m-1] right hand side vector              */
       *s,            /* [0..n] step counter for Newton               */
       *cneu,         /* [0..n] vector with new approximation for     */
                      /* optimal coefficients                         */
       *g,            /* [0..m-1] vector with square roots of weights */
       NeuerFehler,   /* most recent mean square error                */
       Daempfung,     /* Damping factor                               */
       Faktor,        /* aux variables to approximate Jacobi matrix   */
       ck,
       Hk,
       ZHk,
       DQ;


  if (m <= n || *maxIt <= 0 || RelEps <= ZERO)
    return 1;


  vmblock = vminit();                             /* allocate storage */
  Ab   = (REAL **)vmalloc(vmblock, MATRIX, m, n + 1);
  s    = (REAL *) vmalloc(vmblock, VEKTOR, n + 1, 0);
  cneu = (REAL *) vmalloc(vmblock, VEKTOR, n + 1, 0);
  g    = (REAL *) vmalloc(vmblock, VEKTOR, m,     0);
  b    = (REAL *) vmalloc(vmblock, VEKTOR, m,     0);
  if (! vmcomplete(vmblock))                      /* lack of memory ? */
    return 3;


  for (i = 0; i < m; i++)                 /* put roots of weights     */
    g[i] = SQRT(w[i]);                    /* into the vector g        */

  Faktor = POW(MACH_EPS, ONE / THREE);

  /* --- form weighted differences of y-values and those of --------- */
  /* --- the current approximation, store in right hand side b ------ */

  for (i = 0; i < m; i++)
    b[i] = (y[i] - (*PHI)(c, x[i])) * g[i];

  *MiQuFe = quadsum(m, b);            /*  Norm of the right hand side */


  for (Iteration = 1; ; Iteration++)  /* apply Newton method to       */
  {                                   /*  f(c) =                      */
                                      /*  ( y[0]   - PHI(c,x[0])   )  */
                                      /*  (         ...            )  */
                                      /*  ( y[m-1] - PHI(c,x[m-1]) )  */

    if (ablOK)                    /* form the Jacobi matrix via ABL or*/
      for (i = 0; i < m; i++)     /* via central differences          */
        (*ABL)(x[i], c, Ab[i]);
    else
      for (k = 0; k <= n; k++)
      {
        ck = c[k];
        Hk = Faktor;
        if (ck != ZERO)
          Hk *= FABS(ck);
        ZHk = HALF / Hk;
        for (i = 0; i < m; i++)
          c[k]     = ck + Hk,
          DQ       = (*PHI)(c, x[i]),
          c[k]     = ck - Hk,
          Ab[i][k] = (DQ - (*PHI)(c, x[i])) * ZHk;
        c[k] = ck;
      }

    for (i = 0; i < m; i++)              /* adjust rows of matrix by  */
      for (k = 0; k <= n; k++)           /* the weights               */
        Ab[i][k] *= g[i];

    if (house(m, n + 1, Ab, b))    /* compute the direction s for     */
    {                              /* Newton as the solution of the   */
      vmfree(vmblock);             /* linear minimization via a       */
      return 2;                    /* Householder transformation      */
    }
    for (i = 0; i <= n; i++)       /* copy solution in b onto s       */
      s[i] = b[i];

    for (Daempfung = ONE, k = 0;        /* try to find a suitable     */
         k <= 10;                       /* damping factor             */
         k++, Daempfung /= TWO)
    {
      for (i = 0; i <= n; i++)               /* compute a new approx. */
        cneu[i] = c[i] + s[i] * Daempfung;   /* cneu                  */
      for (i = 0; i < m; i++)                /* reevaluate value of   */
        b[i] = (y[i] - (*PHI)(cneu, x[i]))   /* the new function at   */
               * g[i];                       /* cneu                  */
      NeuerFehler = quadsum(m, b);           /* find its norm         */
      if (NeuerFehler <= *MiQuFe)       /* damping factor suitable ?  */
      {
        for (i = 0; i <= n; i++)     /* damp the direction vector     */
          s[i] = cneu[i] - c[i],
          c[i] = cneu[i];            /* store newest approximation in */
        break;                       /* c and end loop                */
      }
    }

    if (k > 10)                       /* damping no good ?            */
    {                                 /* => work without damping !    */
      for (i = 0; i <= n; i++)
        c[i] += s[i];
      for (i = 0; i < m; i++)     /* value of the minimizing function */
        b[i] = (y[i] - (*PHI)(c, x[i])) * g[i];    /* at c            */
      NeuerFehler = quadsum(m, b);      /* Norm of the function value */
    }

    if (quadsum(n + 1, s) <             /* desired accuracy reached ? */
        RelEps * quadsum(n + 1, c))
      break;                               /* stop iteration          */

    if (Iteration > *maxIt)                       /* maxIt exceeded ? */
    {                                             /* stop iteration   */
      vmfree(vmblock);
      return 4;
    }

    *MiQuFe = NeuerFehler;
  }


  *MiQuFe = SQRT(NeuerFehler);                  /* mean square error  */

  *maxIt = Iteration;                 /* actual number of iterations  */


  vmfree(vmblock);

  return 0;
}

/* --------------------------- END  fft.c --------------------------- */
