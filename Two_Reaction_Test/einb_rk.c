#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------ MODULE einb_rk.c ------------------------ */

/***********************************************************************
*                                                                      *
* Solve a first order ordinary differential equation system using the  *
* -------------------------------------------------------------------  *
* Runge-Kutta embedding formulas                                       *
* ------------------------------                                       *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Volker Krueger (FORTRAN 77)                    *
* Adaption:             Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing FORTRAN 77 codes                      *
* Date:                 3.26.1993; 10.30. 1995                         *
*                                                                      *
*                                                                      *
***********************************************************************/

#include <basis.h>    /*  for  REAL, TWO, THREE, FOUR, FIVE, SIX,     */
                      /*       EIGHT, HALF, ONE, ZERO, NINE, SQRT,    */
                      /*       dglsysfnk, norm_max, FALSE, POW, FABS, */
                      /*       TRUE, copy_vector, max, SIGN, LOG,     */
                      /*       TEN, MACH_EPS, POSMAX                  */
#include <vmblock.h>  /*  for  vminit, vmalloc, MATRIX, VEKTOR,       */
                      /*       vmcomplete, vmfree                     */
#include <einb_rk.h>  /*  for  einb_rk                                */



/* ------------------------------------------------------------------ */

#define NE  23          /* total number of embedding formulas         */

typedef struct          /* Structure used to describe each embedding  */
{                       /* formula completely                         */
  int  stufenzahl;      /* level number  m for the embedding formula  */
  int  k1_km;           /* new k1 vector = old km vector after a      */
                        /* successful step                            */
  REAL qg;              /* global error order of lower order method   */
  long maxschritt;      /* maximal number of Runge-Kutta integrations */
  REAL *As;             /* [0..m-1] vector with coefficients for A~   */
  REAL *A;              /* [0..m-1] vector with coefficients for A    */
  REAL *a;              /* [0..m-2] vector with coefficients for a    */
  REAL *b;              /* [0..m*(m-1)/2] vector with the coefficient */
                        /* matrix b rowwise                           */
} koefftyp;
/*.IX{koefftyp}*/

/* ------------------------------------------------------------------ */

static void init_koeff     /* find data of one embedding formula .....*/
/*.IX{init\unt koeff}*/
                      (
                       int      neinb,  /* Number of formula .........*/
                       long     fmax,   /* max. # of calls of dgl() ..*/
                       koefftyp *koeff  /* Data of formula ...........*/
                      )

/***********************************************************************
* Find the level number, the maximal number of integration steps, the  *
* global error order of the lower order method as well as all          *
* coefficients of an embedding formula from its number.                *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* neinb  Number of embedding formula chosen  (0, 1,..., NE -1)         *
* fmax   upper bound for number of function evaluations of the right   *
*        hand side dgl() (needed to compute the maximal number of      *
*        Runge-Kutta integrations)                                     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* koeff  Structure with the essential data for the embedding formula   *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* koefftyp, REAL, TWO, THREE, FOUR, FIVE, SIX, EIGHT, HALF, ONE, ZERO, *
* NINE, SQRT, TEN, NE                                                  *
***********************************************************************/

{
  typedef REAL R;      /* abbreviation for type change                */

  REAL  h;             /* correction factor for embedding formulas    */
                       /*  16 and 22                                  */

  static
    int ok16 = FALSE;  /* Flag to indicate whether special procedure  */
                       /* for formula 16 has been executed            */
  static
    int ok22 = FALSE;  /* Flag to indicate whether special procedure  */
                       /* for formula 22 has been executed            */


  /* -------- combine all coefficients of A~, A, a, b for all ------- */
  /* -------  of the NE embedding formulas into one vector  --------- */

  static REAL koeffizienten[     11 +
                                 24 +
                            3 *  32 +
                            4 *  41 +
                            6 *  51 +
                            3 *  62 +
                                 74 +
                            3 * 116 +
                                167
                           ] =
  {

                  /* ----- 0th embedding formula: method RK3(2)  ---- */
/* A~ */
                  ONE /              SIX,
                  TWO /            THREE,
                  ONE /              SIX,
/* A */
                 ZERO,
                  ONE,
                 ZERO,
/* a */
                 HALF,
                  ONE,
/* b */
                 HALF,
                 -ONE,
                  TWO,

                   /* ---- 1st embedding formula: method RKF4(3) ---- */
/* A~ */
              (R)79.0 /         (R)490.0,
                 ZERO,
            (R)2175.0 /        (R)3626.0,
            (R)2166.0 /        (R)9065.0,
                 ZERO,
/* A */
             (R)229.0 /        (R)1470.0,
                 ZERO,
            (R)1125.0 /        (R)1813.0,
           (R)13718.0 /       (R)81585.0,
                  ONE /          (R)18.0,
/* a */
                  TWO /           (R)7.0,
               (R)7.0 /          (R)15.0,
              (R)35.0 /          (R)38.0,
                  ONE,
/* b */
                  TWO /           (R)7.0,
              (R)77.0 /         (R)900.0,
             (R)343.0 /         (R)900.0,
             (R)805.0 /        (R)1444.0,
          (R)-77175.0 /       (R)54872.0,
           (R)97125.0 /       (R)54872.0,
              (R)79.0 /         (R)490.0,
                 ZERO,
            (R)2175.0 /        (R)3626.0,
            (R)2166.0 /        (R)9065.0,

                  /* ---- 2nd embedding formula: method RKF5(4)  ---- */
/* A~ */
              (R)16.0 /         (R)135.0,
                 ZERO,
            (R)6656.0 /       (R)12825.0,
           (R)28561.0 /       (R)56430.0,
                -NINE /          (R)50.0,
                  TWO /          (R)55.0,
/* A */
              (R)25.0 /         (R)216.0,
                 ZERO,
            (R)1408.0 /        (R)2565.0,
            (R)2197.0 /        (R)4104.0,
                 -ONE /             FIVE,
                 ZERO,
/* a */
                  ONE /             FOUR,
                THREE /            EIGHT,
              (R)12.0 /          (R)13.0,
                  ONE,
                 HALF,
/* b */
                  ONE /             FOUR,
                THREE /          (R)32.0,
                 NINE /          (R)32.0,
            (R)1932.0 /        (R)2197.0,
           (R)-7200.0 /        (R)2197.0,
            (R)7296.0 /        (R)2197.0,
             (R)439.0 /         (R)216.0,
               -EIGHT,
            (R)3680.0 /         (R)513.0,
            (R)-845.0 /        (R)4104.0,
               -EIGHT /          (R)27.0,
                  TWO,
           (R)-3544.0 /        (R)2565.0,
            (R)1859.0 /        (R)4104.0,
             (R)-11.0 /          (R)40.0,

                   /* ---- 3rd embedding formula: method RK5(4)6M --- */
/* A~ */
              (R)19.0 /         (R)216.0,
                 ZERO,
            (R)1000.0 /        (R)2079.0,
            (R)-125.0 /         (R)216.0,
              (R)81.0 /          (R)88.0,
                 FIVE /          (R)56.0,
/* A */
              (R)31.0 /         (R)540.0,
                 ZERO,
             (R)190.0 /         (R)297.0,
            (R)-145.0 /         (R)108.0,
             (R)351.0 /         (R)220.0,
                  ONE /          (R)20.0,
/* a */
                  ONE /             FIVE,
                THREE /              TEN,
                THREE /             FIVE,
                  TWO /            THREE,
                  ONE,
/* b */
                  ONE /             FIVE,
                THREE /          (R)40.0,
                 NINE /          (R)40.0,
                THREE /              TEN,
                -NINE /              TEN,
                  SIX /             FIVE,
             (R)226.0 /         (R)729.0,
             (R)-25.0 /          (R)27.0,
             (R)880.0 /         (R)729.0,
              (R)55.0 /         (R)729.0,
            (R)-181.0 /         (R)270.0,
                 FIVE /              TWO,
            (R)-266.0 /         (R)297.0,
             (R)-91.0 /          (R)27.0,
             (R)189.0 /          (R)55.0,

                  /* ---- 4th embedding formula: method RKE5(4)  ---- */
/* A~ */
              (R)14.0 /         (R)336.0,
                 ZERO,
                 ZERO,
              (R)35.0 /         (R)336.0,
             (R)162.0 /         (R)336.0,
             (R)125.0 /         (R)336.0,
/* A */
                  ONE /              SIX,
                 ZERO,
                 FOUR /              SIX,
                  ONE /              SIX,
                 ZERO,
                 ZERO,
/* a */
                 HALF,
                 HALF,
                  ONE,
                  TWO /            THREE,
                  ONE /             FIVE,
/* b */
                 HALF,
                  ONE /             FOUR,
                  ONE /             FOUR,
                 ZERO,
                 -ONE,
                  TWO,
               (R)7.0 /          (R)27.0,
                  TEN /          (R)27.0,
                 ZERO,
                  ONE /          (R)27.0,
              (R)28.0 /         (R)625.0,
            (R)-125.0 /         (R)625.0,
             (R)546.0 /         (R)625.0,
              (R)54.0 /         (R)625.0,
            (R)-378.0 /         (R)625.0,

  /* ---- 5th embedding formula; method  HIHA5                  ---- */
/* A~ */
                  ONE /          (R)12.0,
                 ZERO,
              (R)27.0 /          (R)32.0,
                -FOUR /            THREE,
             (R)125.0 /          (R)96.0,
                 FIVE /          (R)48.0,
                 ZERO,
/* A */
                  TWO /          (R)15.0,
                 ZERO,
              (R)27.0 /          (R)80.0,
                 -TWO /          (R)15.0,
              (R)25.0 /          (R)48.0,
                  ONE /          (R)24.0,
                  ONE /              TEN,
/* a */
                  TWO /             NINE,
                  ONE /            THREE,
                 HALF,
                THREE /             FIVE,
                  ONE,
                  ONE,
/* b */
                  TWO /             NINE,
                  ONE /          (R)12.0,
                  ONE /             FOUR,
                  ONE /            EIGHT,
                 ZERO,
                THREE /            EIGHT,
              (R)91.0 /         (R)500.0,
             (R)-27.0 /         (R)100.0,
              (R)78.0 /         (R)125.0,
                EIGHT /         (R)125.0,
             (R)-11.0 /          (R)20.0,
              (R)27.0 /          (R)20.0,
              (R)12.0 /             FIVE,
             (R)-36.0 /             FIVE,
                 FIVE,
                  ONE /          (R)12.0,
                 ZERO,
              (R)27.0 /          (R)32.0,
                -FOUR /            THREE,
             (R)125.0 /          (R)96.0,
                 FIVE /          (R)48.0,

                  /* ---- 6th embedding formula: method RK5(4)7S  --- */
/* A~ */
              (R)19.0 /         (R)200.0,
                 ZERO,
                THREE /             FIVE,
            (R)-243.0 /         (R)400.0,
              (R)33.0 /          (R)40.0,
               (R)7.0 /          (R)80.0,
                 ZERO,
/* A */
             (R)431.0 /        (R)5000.0,
                 ZERO,
             (R)333.0 /         (R)500.0,
           (R)-7857.0 /       (R)10000.0,
             (R)957.0 /        (R)1000.0,
             (R)193.0 /        (R)2000.0,
                 -ONE /          (R)50.0,
/* a */
                  TWO /             NINE,
                  ONE /            THREE,
                 FIVE /             NINE,
                  TWO /            THREE,
                  ONE,
                  ONE,
/* b */
                  TWO /             NINE,
                  ONE /          (R)12.0,
                  ONE /             FOUR,
              (R)55.0 /         (R)324.0,
             (R)-25.0 /         (R)108.0,
              (R)50.0 /          (R)81.0,
              (R)83.0 /         (R)330.0,
             (R)-13.0 /          (R)22.0,
              (R)61.0 /          (R)66.0,
                 NINE /         (R)110.0,
             (R)-19.0 /          (R)28.0,
                 NINE /             FOUR,
                  ONE /           (R)7.0,
             (R)-27.0 /           (R)7.0,
              (R)22.0 /           (R)7.0,
              (R)19.0 /         (R)200.0,
                 ZERO,
                THREE /             FIVE,
            (R)-243.0 /         (R)400.0,
              (R)33.0 /          (R)40.0,
               (R)7.0 /          (R)80.0,

                  /* ---- 7th embedding formula: method RK5(4)7M  --- */
/* A~ */
              (R)35.0 /         (R)384.0,
                 ZERO,
             (R)500.0 /        (R)1113.0,
             (R)125.0 /         (R)192.0,
           (R)-2187.0 /        (R)6784.0,
              (R)11.0 /          (R)84.0,
                 ZERO,
/* A */
            (R)5179.0 /       (R)57600.0,
                 ZERO,
            (R)7571.0 /       (R)16695.0,
             (R)393.0 /         (R)640.0,
          (R)-92097.0 /      (R)339200.0,
             (R)187.0 /        (R)2100.0,
                  ONE /          (R)40.0,
/* a */
                  ONE /             FIVE,
                THREE /              TEN,
                 FOUR /             FIVE,
                EIGHT /             NINE,
                  ONE,
                  ONE,
/* b */
                  ONE /             FIVE,
                THREE /          (R)40.0,
                 NINE /          (R)40.0,
              (R)44.0 /          (R)45.0,
             (R)-56.0 /          (R)15.0,
              (R)32.0 /             NINE,
           (R)19372.0 /        (R)6561.0,
          (R)-25360.0 /        (R)2187.0,
           (R)64448.0 /        (R)6561.0,
            (R)-212.0 /         (R)729.0,
            (R)9017.0 /        (R)3168.0,
            (R)-355.0 /          (R)33.0,
           (R)46732.0 /        (R)5247.0,
              (R)49.0 /         (R)176.0,
           (R)-5103.0 /       (R)18656.0,
              (R)35.0 /         (R)384.0,
                 ZERO,
             (R)500.0 /        (R)1113.0,
             (R)125.0 /         (R)192.0,
           (R)-2187.0 /        (R)6784.0,
              (R)11.0 /          (R)84.0,

  /* ---- 8th embedding formula; method  RK5(4)7C          ---- */
/* A~ */
              (R)35.0 /         (R)432.0,
                 ZERO,
            (R)8500.0 /       (R)14553.0,
          (R)-28561.0 /       (R)84672.0,
             (R)405.0 /         (R)704.0,
              (R)19.0 /         (R)196.0,
                 ZERO,
/* A */
              (R)11.0 /         (R)108.0,
                 ZERO,
            (R)6250.0 /       (R)14553.0,
           (R)-2197.0 /       (R)21168.0,
              (R)81.0 /         (R)176.0,
             (R)171.0 /        (R)1960.0,
                  ONE /          (R)40.0,
/* a */
                  ONE /             FIVE,
                THREE /              TEN,
                  SIX /          (R)13.0,
                  TWO /            THREE,
                  ONE,
                  ONE,
/* b */
                  ONE /             FIVE,
                THREE /          (R)40.0,
                 NINE /          (R)40.0,
             (R)264.0 /        (R)2197.0,
             (R)-90.0 /        (R)2197.0,
             (R)840.0 /        (R)2197.0,
             (R)932.0 /        (R)3645.0,
             (R)-14.0 /          (R)27.0,
            (R)3256.0 /        (R)5103.0,
            (R)7436.0 /       (R)25515.0,
            (R)-367.0 /         (R)513.0,
              (R)30.0 /          (R)19.0,
            (R)9940.0 /        (R)5643.0,
          (R)-29575.0 /        (R)8208.0,
            (R)6615.0 /        (R)3344.0,
              (R)35.0 /         (R)432.0,
                 ZERO,
            (R)8500.0 /       (R)14553.0,
          (R)-28561.0 /       (R)84672.0,
             (R)405.0 /         (R)704.0,
              (R)19.0 /         (R)196.0,

                  /* ---- 9th embedding formula: method RK6(5)8M  --- */
/* A~ */
              (R)61.0 /         (R)864.0,
                 ZERO,
           (R)98415.0 /      (R)321776.0,
           (R)16807.0 /      (R)146016.0,
            (R)1375.0 /        (R)7344.0,
            (R)1375.0 /        (R)5408.0,
             (R)-37.0 /        (R)1120.0,
                  ONE /              TEN,
/* A */
             (R)821.0 /       (R)10800.0,
                 ZERO,
           (R)19683.0 /       (R)71825.0,
          (R)175273.0 /      (R)912600.0,
             (R)395.0 /        (R)3672.0,
             (R)785.0 /        (R)2704.0,
                THREE /          (R)50.0,
                 ZERO,
/* a */
                  ONE /              TEN,
                  TWO /             NINE,
                THREE /           (R)7.0,
                THREE /             FIVE,
                 FOUR /             FIVE,
                  ONE,
                  ONE,
/* b */
                  ONE /              TEN,
                 -TWO /          (R)81.0,
              (R)20.0 /          (R)81.0,
             (R)615.0 /        (R)1372.0,
            (R)-270.0 /         (R)343.0,
            (R)1053.0 /        (R)1372.0,
            (R)3243.0 /        (R)5500.0,
             (R)-54.0 /          (R)55.0,
           (R)50949.0 /       (R)71500.0,
            (R)4998.0 /       (R)17875.0,
          (R)-26492.0 /       (R)37125.0,
              (R)72.0 /          (R)55.0,
            (R)2808.0 /       (R)23375.0,
          (R)-24206.0 /       (R)37125.0,
             (R)338.0 /         (R)459.0,
            (R)5561.0 /        (R)2376.0,
             (R)-35.0 /          (R)11.0,
          (R)-24117.0 /       (R)31603.0,
          (R)899983.0 /      (R)200772.0,
           (R)-5225.0 /        (R)1836.0,
            (R)3925.0 /        (R)4056.0,
          (R)465467.0 /      (R)266112.0,
           (R)-2945.0 /        (R)1232.0,
        (R)-5610201.0 /    (R)14158144.0,
        (R)10513573.0 /     (R)3212352.0,
         (R)-424325.0 /      (R)205632.0,
          (R)376225.0 /      (R)454272.0,
                 ZERO,

  /* ---- 10th embedding formula; method  RK6(5)8S        ---- */
/* A~ */
              (R)29.0 /         (R)324.0,
                 ZERO,
            (R)3400.0 /        (R)7371.0,
          (R)-16807.0 /       (R)25272.0,
            (R)-125.0 /        (R)1944.0,
              (R)25.0 /          (R)24.0,
                  ONE /          (R)84.0,
                  ONE /            EIGHT,
/* A */
            (R)2041.0 /       (R)21600.0,
                 ZERO,
             (R)748.0 /        (R)1755.0,
           (R)-2401.0 /       (R)46800.0,
              (R)11.0 /         (R)108.0,
              (R)59.0 /         (R)160.0,
                THREE /          (R)50.0,
                 ZERO,
/* a */
                  ONE /             FOUR,
                THREE /              TEN,
                  SIX /           (R)7.0,
                THREE /             FIVE,
                 FOUR /             FIVE,
                  ONE,
                  ONE,
/* b */
                  ONE /             FOUR,
                THREE /          (R)25.0,
                 NINE /          (R)50.0,
             (R)102.0 /         (R)343.0,
           (R)-1368.0 /         (R)343.0,
            (R)1560.0 /         (R)343.0,
               -THREE /         (R)100.0,
              (R)36.0 /          (R)25.0,
             (R)-12.0 /          (R)13.0,
             (R)147.0 /        (R)1300.0,
              (R)37.0 /         (R)225.0,
             (R)-48.0 /          (R)25.0,
             (R)872.0 /         (R)351.0,
              (R)49.0 /        (R)1053.0,
                  TWO /          (R)81.0,
              (R)11.0 /         (R)648.0,
              (R)14.0 /            THREE,
          (R)-10193.0 /        (R)2106.0,
          (R)-30331.0 /       (R)50544.0,
            (R)1025.0 /        (R)1944.0,
              (R)59.0 /          (R)48.0,
             (R)796.0 /        (R)1701.0,
            (R)-352.0 /          (R)63.0,
          (R)134096.0 /       (R)22113.0,
          (R)-78281.0 /       (R)75816.0,
           (R)-9425.0 /       (R)20412.0,
             (R)781.0 /         (R)504.0,
                 ZERO,

  /* ---- 11th embedding formula; method  RK6(5)8C         ---- */
/* A~ */
                  ONE /          (R)12.0,
                 ZERO,
            (R)-216.0 /        (R)1235.0,
            (R)6561.0 /       (R)12376.0,
            (R)1375.0 /        (R)5304.0,
            (R)1375.0 /        (R)5928.0,
                -FIVE /         (R)168.0,
                  ONE /              TEN,
/* A */
             (R)163.0 /        (R)1440.0,
                 ZERO,
           (R)-2628.0 /        (R)6175.0,
           (R)13851.0 /       (R)17680.0,
            (R)1525.0 /        (R)7956.0,
            (R)6575.0 /       (R)23712.0,
                THREE /          (R)50.0,
                 ZERO,
/* a */
                  ONE /              TEN,
                  ONE /              SIX,
                  TWO /             NINE,
                THREE /             FIVE,
                 FOUR /             FIVE,
                  ONE,
                  ONE,
/* b */
                  ONE /              TEN,
                  ONE /          (R)36.0,
                 FIVE /          (R)36.0,
                  TEN /         (R)243.0,
              (R)20.0 /         (R)243.0,
                EIGHT /          (R)81.0,
            (R)4047.0 /        (R)5500.0,
             (R)-18.0 /          (R)55.0,
           (R)-4212.0 /        (R)1375.0,
           (R)17901.0 /        (R)5500.0,
           (R)-5587.0 /        (R)4125.0,
              (R)24.0 /          (R)55.0,
            (R)9576.0 /        (R)1375.0,
         (R)-140049.0 /       (R)23375.0,
              (R)38.0 /          (R)51.0,
           (R)12961.0 /        (R)2376.0,
             (R)-35.0 /          (R)33.0,
         (R)-160845.0 /        (R)5434.0,
         (R)1067565.0 /       (R)38896.0,
         (R)-103375.0 /       (R)47736.0,
           (R)32875.0 /       (R)35568.0,
          (R)702799.0 /      (R)199584.0,
           (R)-1865.0 /        (R)2772.0,
        (R)-2891375.0 /      (R)152152.0,
        (R)19332955.0 /     (R)1089088.0,
        (R)-5356375.0 /     (R)4009824.0,
         (R)2207875.0 /     (R)2987712.0,
                 ZERO,

                 /* ---- 12th embedding formula: method RKV6(5)  ---- */
/* A~ */
              (R)57.0 /         (R)640.0,
                 ZERO,
             (R)-16.0 /          (R)65.0,
            (R)1377.0 /        (R)2240.0,
             (R)121.0 /         (R)320.0,
                 ZERO,
             (R)891.0 /        (R)8320.0,
                  TWO /          (R)35.0,
/* A */
                THREE /          (R)80.0,
                 ZERO,
                 FOUR /          (R)25.0,
             (R)243.0 /        (R)1120.0,
              (R)77.0 /         (R)160.0,
              (R)73.0 /         (R)700.0,
                 ZERO,
                 ZERO,
/* a */
                  ONE /          (R)18.0,
                  ONE /              SIX,
                  TWO /             NINE,
                  TWO /            THREE,
                  ONE,
                EIGHT /             NINE,
                  ONE,
/* b */
                  ONE /          (R)18.0,
                 -ONE /          (R)12.0,
                  ONE /             FOUR,
                 -TWO /          (R)81.0,
                 FOUR /          (R)27.0,
                EIGHT /          (R)81.0,
              (R)40.0 /          (R)33.0,
                -FOUR /          (R)11.0,
             (R)-56.0 /          (R)11.0,
              (R)54.0 /          (R)11.0,
            (R)-369.0 /          (R)73.0,
              (R)72.0 /          (R)73.0,
            (R)5380.0 /         (R)219.0,
          (R)-12285.0 /         (R)584.0,
            (R)2695.0 /        (R)1752.0,
           (R)-8716.0 /         (R)891.0,
             (R)656.0 /         (R)297.0,
           (R)39520.0 /         (R)891.0,
            (R)-416.0 /          (R)11.0,
              (R)52.0 /          (R)27.0,
                 ZERO,
            (R)3015.0 /         (R)256.0,
                -NINE /             FOUR,
           (R)-4219.0 /          (R)78.0,
            (R)5985.0 /         (R)128.0,
            (R)-539.0 /         (R)384.0,
                 ZERO,
             (R)693.0 /        (R)3328.0,

                 /* ---- 13th embedding formula: method RKF6(5)  ---- */
/* A~ */
               (R)7.0 /        (R)1408.0,
                 ZERO,
            (R)1125.0 /        (R)2816.0,
                 NINE /          (R)32.0,
             (R)125.0 /         (R)768.0,
                 ZERO,
                 FIVE /          (R)66.0,
                 FIVE /          (R)66.0,
/* A */
              (R)31.0 /         (R)384.0,
                 ZERO,
            (R)1125.0 /        (R)2816.0,
                 NINE /          (R)32.0,
             (R)125.0 /         (R)768.0,
                 FIVE /          (R)66.0,
                 ZERO,
                 ZERO,
/* a */
                  ONE /              SIX,
                 FOUR /          (R)15.0,
                  TWO /            THREE,
                 FOUR /             FIVE,
                  ONE,
                 ZERO,
                  ONE,
/* b */
                  ONE /              SIX,
                 FOUR /          (R)75.0,
              (R)16.0 /          (R)75.0,
                 FIVE /              SIX,
               -EIGHT /            THREE,
                 FIVE /              TWO,
               -EIGHT /             FIVE,
             (R)144.0 /          (R)25.0,
                -FOUR,
              (R)16.0 /          (R)25.0,
             (R)361.0 /         (R)320.0,
             (R)-18.0 /             FIVE,
             (R)407.0 /         (R)128.0,
             (R)-11.0 /          (R)80.0,
              (R)55.0 /         (R)128.0,
             (R)-11.0 /         (R)640.0,
                 ZERO,
              (R)11.0 /         (R)256.0,
             (R)-11.0 /         (R)160.0,
              (R)11.0 /         (R)256.0,
                 ZERO,
              (R)93.0 /         (R)640.0,
             (R)-18.0 /             FIVE,
             (R)803.0 /         (R)256.0,
             (R)-11.0 /         (R)160.0,
              (R)99.0 /         (R)256.0,
                 ZERO,
                  ONE,

  /* ---- 14th embedding formula; method  RKF6(5)B         ---- */
/* A~ */
             (R)167.0 /        (R)2080.0,
                 ZERO,
           (R)91375.0 /      (R)229152.0,
              (R)41.0 /         (R)144.0,
        (R)43046721.0 /   (R)269010560.0,
                 ZERO,
             (R)125.0 /      (R)150192.0,
             (R)743.0 /        (R)9856.0,
/* A */
              (R)21.0 /         (R)260.0,
                 ZERO,
            (R)7625.0 /       (R)19096.0,
              (R)25.0 /          (R)88.0,
         (R)1594323.0 /     (R)9929920.0,
              (R)53.0 /         (R)704.0,
                 ZERO,
                 ZERO,
/* a */
                  ONE /              SIX,
                 FOUR /          (R)15.0,
                  TWO /            THREE,
              (R)65.0 /          (R)81.0,
                  ONE,
                  ONE /          (R)15.0,
                  ONE,
/* b */
                  ONE /              SIX,
                 FOUR /          (R)75.0,
              (R)16.0 /          (R)75.0,
                 FIVE /              SIX,
               -EIGHT /            THREE,
                 FIVE /              TWO,
          (R)-43745.0 /       (R)26244.0,
          (R)353860.0 /       (R)59049.0,
         (R)-493675.0 /      (R)118098.0,
          (R)155155.0 /      (R)236196.0,
           (R)16539.0 /       (R)13780.0,
            (R)-204.0 /          (R)53.0,
          (R)232595.0 /       (R)69006.0,
             (R)-91.0 /         (R)636.0,
          (R)314928.0 /      (R)747565.0,
         (R)-409063.0 /      (R)195000.0,
             (R)124.0 /          (R)75.0,
           (R)25357.0 /        (R)8680.0,
           (R)-8928.0 /        (R)1375.0,
        (R)79184709.0 /    (R)19394375.0,
                 ZERO,
          (R)231213.0 /      (R)193180.0,
           (R)-2820.0 /         (R)743.0,
         (R)9512525.0 /     (R)2902158.0,
           (R)-5113.0 /       (R)80244.0,
       (R)584348904.0 /  (R)1561522235.0,
                 ZERO,
           (R)30800.0 /     (R)2989089.0,

  /* ---- 15th embedding formula; method RKC6(5)     ---- */
/* A~ */
        (R)17572349.0 /   (R)289262523.0,
                 ZERO,
        (R)57513011.0 /   (R)201864250.0,
        (R)15587306.0 /   (R)354501571.0,
        (R)71783021.0 /   (R)234982865.0,
        (R)29672000.0 /   (R)180480167.0,
        (R)65567621.0 /   (R)127060952.0,
       (R)-79074570.0 /   (R)210557597.0,
                 ZERO,
/* A */
        (R)15231665.0 /   (R)510830334.0,
                 ZERO,
        (R)59452991.0 /   (R)116050448.0,
       (R)-28398517.0 /   (R)122437738.0,
        (R)56673824.0 /   (R)137010559.0,
        (R)68003849.0 /   (R)426673583.0,
         (R)7097631.0 /    (R)37564021.0,
       (R)-71226429.0 /   (R)583093742.0,
                  ONE /          (R)20.0,
/* a */
                  TWO /          (R)15.0,
                  ONE /             FIVE,
                THREE /              TEN,
                 14.0 /          (R)25.0,
                 19.0 /          (R)25.0,
        (R)35226607.0 /    (R)35688279.0,
                  ONE,
                  ONE,
/* b */
                  TWO /          (R)15.0,
                  ONE /          (R)20.0,
                THREE /          (R)20.0,
                THREE /          (R)40.0,
                 ZERO,
                 NINE /          (R)40.0,
        (R)86727015.0 /   (R)196851553.0,
       (R)-60129073.0 /    (R)52624712.0,
       (R)957436434.0 /  (R)1378352377.0,
        (R)83886832.0 /   (R)147842441.0,
       (R)-86860849.0 /    (R)45628967.0,
       (R)111022885.0 /    (R)25716487.0,
       (R)108046682.0 /   (R)101167669.0,
      (R)-141756746.0 /    (R)36005461.0,
        (R)73139862.0 /    (R)60170633.0,
        (R)77759591.0 /    (R)16096467.0,
       (R)-49252809.0 /     (R)6452555.0,
      (R)-381680111.0 /    (R)51572984.0,
       (R)879269579.0 /    (R)66788831.0,
       (R)-90453121.0 /    (R)33722162.0,
       (R)111179552.0 /   (R)157155827.0,
       (R)237564263.0 /    (R)39280295.0,
      (R)-100523239.0 /    (R)10677940.0,
      (R)-265574846.0 /    (R)27330247.0,
       (R)317978411.0 /    (R)18988713.0,
      (R)-124494385.0 /    (R)35453627.0,
        (R)86822444.0 /   (R)100138635.0,
       (R)-12873523.0 /   (R)724232625.0,
        (R)17572349.0 /   (R)289262523.0,
                 ZERO,
        (R)57513011.0 /   (R)201864250.0,
        (R)15587306.0 /   (R)354501571.0,
        (R)71783021.0 /   (R)234982865.0,
        (R)29672000.0 /   (R)180480167.0,
        (R)65567621.0 /   (R)127060952.0,
       (R)-79074570.0 /   (R)210557597.0,

  /* ---- 16th embedding formula; method  RKV6(5)9A        ---- */
/* A~ */
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
/* A */
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
/* a */
                  ONE /            EIGHT,
                 ZERO,
                 ZERO,
                 NINE /          (R)16.0,
                  ONE /              TWO,
                 NINE /              TEN,
                  ONE,
                  ONE,
/* b */
                  ONE / EIGHT,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,

  /* ---- 17th embedding formula; method  RKV6(5)9B        ---- */
/* A~ */
              (R)11.0 /         (R)144.0,
                 ZERO,
                 ZERO,
             (R)256.0 /         (R)693.0,
                 ZERO,
             (R)125.0 /         (R)504.0,
             (R)125.0 /         (R)528.0,
                 FIVE /          (R)72.0,
                 ZERO,
/* A */
                  ONE /          (R)18.0,
                 ZERO,
                 ZERO,
              (R)32.0 /          (R)63.0,
                 -TWO /            THREE,
             (R)125.0 /         (R)126.0,
                 ZERO,
                -FIVE /          (R)63.0,
                 FOUR /          (R)21.0,
/* a */
                  ONE /            EIGHT,
                  ONE /              SIX,
                  ONE /             FOUR,
                  ONE /              TWO,
                THREE /             FIVE,
                 FOUR /             FIVE,
                  ONE,
                  ONE,
/* b */
                  ONE /            EIGHT,
                  ONE /          (R)18.0,
                  ONE /             NINE,
                  ONE /          (R)16.0,
                 ZERO,
                THREE /          (R)16.0,
                  ONE /             FOUR,
                 ZERO,
               -THREE /             FOUR,
                  ONE,
             (R)134.0 /         (R)625.0,
                 ZERO,
            (R)-333.0 /         (R)625.0,
             (R)476.0 /         (R)625.0,
              (R)98.0 /         (R)625.0,
             (R)-98.0 /        (R)1875.0,
                 ZERO,
              (R)12.0 /         (R)625.0,
           (R)10736.0 /       (R)13125.0,
           (R)-1936.0 /        (R)1875.0,
              (R)22.0 /          (R)21.0,
                 NINE /          (R)50.0,
                 ZERO,
              (R)21.0 /          (R)25.0,
           (R)-2924.0 /        (R)1925.0,
              (R)74.0 /          (R)25.0,
             (R)-15.0 /           (R)7.0,
              (R)15.0 /          (R)22.0,
              (R)11.0 /         (R)144.0,
                 ZERO,
                 ZERO,
             (R)256.0 /         (R)693.0,
                 ZERO,
             (R)125.0 /         (R)504.0,
             (R)125.0 /         (R)528.0,
                 FIVE /          (R)72.0,

                  /* ---- 18th embedding formula: method RKV7(6)  --- */
/* A~ */
            (R)2881.0 /       (R)40320.0,
                 ZERO,
                 ZERO,
            (R)1216.0 /        (R)2961.0,
           (R)-2624.0 /        (R)4095.0,
        (R)24137569.0 /    (R)57482880.0,
                -FOUR /          (R)21.0,
                 ZERO,
            (R)4131.0 /        (R)3920.0,
            (R)-157.0 /        (R)1260.0,
/* A */
               (R)7.0 /          (R)90.0,
                 ZERO,
                 ZERO,
              (R)16.0 /          (R)45.0,
              (R)16.0 /          (R)45.0,
                 ZERO,
                  TWO /          (R)15.0,
               (R)7.0 /          (R)90.0,
                 ZERO,
                 ZERO,
/* a */
                  ONE /          (R)12.0,
                  ONE /              SIX,
                  ONE /             FOUR,
                THREE /             FOUR,
              (R)16.0 /          (R)17.0,
                 HALF,
                  ONE,
                  TWO /            THREE,
                  ONE,
/* b */
                  ONE /          (R)12.0,
                 ZERO,
                  ONE /              SIX,
                  ONE /          (R)16.0,
                 ZERO,
                THREE /          (R)16.0,
              (R)21.0 /          (R)16.0,
                 ZERO,
             (R)-81.0 /          (R)16.0,
                 NINE /              TWO,
         (R)1344688.0 /      (R)250563.0,
                 ZERO,
        (R)-1709184.0 /       (R)83521.0,
         (R)1365632.0 /       (R)83521.0,
          (R)-78208.0 /      (R)250563.0,
            (R)-559.0 /         (R)384.0,
                 ZERO,
                  SIX,
            (R)-204.0 /          (R)47.0,
              (R)14.0 /          (R)39.0,
           (R)-4913.0 /       (R)78208.0,
            (R)-625.0 /         (R)224.0,
                 ZERO,
              (R)12.0,
            (R)-456.0 /          (R)47.0,
              (R)48.0 /          (R)91.0,
           (R)14739.0 /      (R)136864.0,
                  SIX /           (R)7.0,
          (R)-12253.0 /       (R)99144.0,
                 ZERO,
              (R)16.0 /          (R)27.0,
              (R)16.0 /         (R)459.0,
           (R)29072.0 /      (R)161109.0,
           (R)-2023.0 /       (R)75816.0,
             (R)112.0 /       (R)12393.0,
                 ZERO,
           (R)30517.0 /        (R)2512.0,
                 ZERO,
           (R)-7296.0 /         (R)157.0,
          (R)268728.0 /        (R)7379.0,
            (R)2472.0 /        (R)2041.0,
        (R)-3522621.0 /    (R)10743824.0,
             (R)132.0 /         (R)157.0,
                 ZERO,
          (R)-12393.0 /        (R)4396.0,

                  /* --- 19th embedding formula: method RK8(7)13M --- */
/* A~ */
        (R)14005451.0 /   (R)335480064.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
       (R)-59238493.0 /  (R)1068277825.0,
       (R)181606767.0 /   (R)758867731.0,
       (R)561292985.0 /   (R)797845732.0,
     (R)-1041891430.0 /  (R)1371343529.0,
       (R)760417239.0 /  (R)1151165299.0,
       (R)118820643.0 /   (R)751138087.0,
      (R)-528747749.0 /  (R)2220607170.0,
                  ONE /             FOUR,
/* A */
        (R)13451932.0 /   (R)455176623.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
      (R)-808719846.0 /   (R)976000145.0,
      (R)1757004468.0 /  (R)5645159321.0,
       (R)656045339.0 /   (R)265891186.0,
     (R)-3867574721.0 /  (R)1518517206.0,
       (R)465885868.0 /   (R)322736535.0,
        (R)53011238.0 /   (R)667516719.0,
                  TWO /          (R)45.0,
                 ZERO,
/* a */
                  ONE /          (R)18.0,
                  ONE /          (R)12.0,
                  ONE /            EIGHT,
                 FIVE /          (R)16.0,
                THREE /            EIGHT,
              (R)59.0 /         (R)400.0,
              (R)93.0 /         (R)200.0,
      (R)5490023248.0 /  (R)9719169821.0,
              (R)13.0 /          (R)20.0,
      (R)1201146811.0 /  (R)1299019798.0,
                  ONE,
                  ONE,
/* b */
                  ONE /          (R)18.0,
                  ONE /          (R)48.0,
                  ONE /          (R)16.0,
                  ONE /          (R)32.0,
                 ZERO,
                THREE /          (R)32.0,
                 FIVE /          (R)16.0,
                 ZERO,
             (R)-75.0 /          (R)64.0,
              (R)75.0 /          (R)64.0,
                THREE /          (R)80.0,
                 ZERO,
                 ZERO,
                THREE /          (R)16.0,
                THREE /          (R)20.0,
        (R)29443841.0 /   (R)614563906.0,
                 ZERO,
                 ZERO,
        (R)77736538.0 /   (R)692538347.0,
       (R)-28693883.0 /  (R)1125000000.0,
        (R)23124283.0 /  (R)1800000000.0,
        (R)16016141.0 /   (R)946692911.0,
                 ZERO,
                 ZERO,
        (R)61564180.0 /   (R)158732637.0,
        (R)22789713.0 /   (R)633445777.0,
       (R)545815736.0 /  (R)2771057229.0,
      (R)-180193667.0 /  (R)1043307555.0,
        (R)39632708.0 /   (R)573591083.0,
                 ZERO,
                 ZERO,
      (R)-433636366.0 /   (R)683701615.0,
      (R)-421739975.0 /  (R)2616292301.0,
       (R)100302831.0 /   (R)723423059.0,
       (R)790204164.0 /   (R)839813087.0,
       (R)800635310.0 /  (R)3783071287.0,
       (R)246121993.0 /  (R)1340847787.0,
                 ZERO,
                 ZERO,
    (R)-37695042795.0 / (R)15268766246.0,
      (R)-309121744.0 /  (R)1061227803.0,
       (R)-12992083.0 /   (R)490766935.0,
      (R)6005943493.0 /  (R)2108947869.0,
       (R)393006217.0 /  (R)1396673457.0,
       (R)123872331.0 /  (R)1001029789.0,
     (R)-1028468189.0 /   (R)846180014.0,
                 ZERO,
                 ZERO,
      (R)8478235783.0 /   (R)508512852.0,
      (R)1311729495.0 /  (R)1432422823.0,
    (R)-10304129995.0 /  (R)1701304382.0,
    (R)-48777925059.0 /  (R)3047939560.0,
     (R)15336726248.0 /  (R)1032824649.0,
    (R)-45442868181.0 /  (R)3398467696.0,
      (R)3065993473.0 /   (R)597172653.0,
       (R)185892177.0 /   (R)718116043.0,
                 ZERO,
                 ZERO,
     (R)-3185094517.0 /   (R)667107341.0,
      (R)-477755414.0 /  (R)1098053517.0,
      (R)-703635378.0 /   (R)230739211.0,
      (R)5731566787.0 /  (R)1027545527.0,
      (R)5232866602.0 /   (R)850066563.0,
     (R)-4093664535.0 /   (R)808688257.0,
      (R)3962137247.0 /  (R)1805957418.0,
        (R)65686358.0 /   (R)487910083.0,
       (R)403863854.0 /   (R)491063109.0,
                 ZERO,
                 ZERO,
     (R)-5068492393.0 /   (R)434740067.0,
      (R)-411421997.0 /   (R)543043805.0,
       (R)652783627.0 /   (R)914296604.0,
     (R)11173962825.0 /   (R)925320556.0,
    (R)-13158990841.0 /  (R)6184727034.0,
      (R)3936647629.0 /  (R)1978049680.0,
      (R)-160528059.0 /   (R)685178525.0,
       (R)248638103.0 /  (R)1413531060.0,
                 ZERO,

                   /* ---- 20th embedding formula: method RKF8(7) --- */
/* A~ */
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
              (R)34.0 /         (R)105.0,
                 NINE /          (R)35.0,
                 NINE /          (R)35.0,
                 NINE /         (R)280.0,
                 NINE /         (R)280.0,
                 ZERO,
              (R)41.0 /         (R)840.0,
              (R)41.0 /         (R)840.0,
/* A */
              (R)41.0 /         (R)840.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
              (R)34.0 /         (R)105.0,
                 NINE /          (R)35.0,
                 NINE /          (R)35.0,
                 NINE /         (R)280.0,
                 NINE /         (R)280.0,
              (R)41.0 /         (R)840.0,
                 ZERO,
                 ZERO,
/* a */
                  TWO /          (R)27.0,
                  ONE /             NINE,
                  ONE /              SIX,
                 FIVE /          (R)12.0,
                 HALF,
                 FIVE /              SIX,
                  ONE /              SIX,
                  TWO /            THREE,
                  ONE /            THREE,
                  ONE,
                 ZERO,
                  ONE,
/* b */
                  TWO /          (R)27.0,
                  ONE /          (R)36.0,
                  ONE /          (R)12.0,
                  ONE /          (R)24.0,
                 ZERO,
                  ONE /            EIGHT,
                 FIVE /          (R)12.0,
                 ZERO,
             (R)-25.0 /          (R)16.0,
              (R)25.0 /          (R)16.0,
                  ONE /          (R)20.0,
                 ZERO,
                 ZERO,
                  ONE /             FOUR,
                  ONE /             FIVE,
             (R)-25.0 /         (R)108.0,
                 ZERO,
                 ZERO,
             (R)125.0 /         (R)108.0,
             (R)-65.0 /          (R)27.0,
             (R)125.0 /          (R)54.0,
              (R)31.0 /         (R)300.0,
                 ZERO,
                 ZERO,
                 ZERO,
              (R)61.0 /         (R)225.0,
                 -TWO /             NINE,
              (R)13.0 /         (R)900.0,
                  TWO,
                 ZERO,
                 ZERO,
             (R)-53.0 /              SIX,
             (R)704.0 /          (R)45.0,
            (R)-107.0 /             NINE,
              (R)67.0 /          (R)90.0,
                THREE,
             (R)-91.0 /         (R)108.0,
                 ZERO,
                 ZERO,
              (R)23.0 /         (R)108.0,
            (R)-976.0 /         (R)135.0,
             (R)311.0 /          (R)54.0,
             (R)-19.0 /          (R)60.0,
              (R)17.0 /              SIX,
                 -ONE /          (R)12.0,
            (R)2383.0 /        (R)4100.0,
                 ZERO,
                 ZERO,
            (R)-341.0 /         (R)164.0,
            (R)4496.0 /        (R)1025.0,
            (R)-301.0 /          (R)82.0,
            (R)2133.0 /        (R)4100.0,
              (R)45.0 /          (R)82.0,
              (R)45.0 /         (R)164.0,
              (R)18.0 /          (R)41.0,
                THREE /         (R)205.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 -SIX /          (R)41.0,
               -THREE /         (R)205.0,
               -THREE /          (R)41.0,
                THREE /          (R)41.0,
                  SIX /          (R)41.0,
                 ZERO,
           (R)-1777.0 /        (R)4100.0,
                 ZERO,
                 ZERO,
            (R)-341.0 /         (R)164.0,
            (R)4496.0 /        (R)1025.0,
            (R)-289.0 /          (R)82.0,
            (R)2193.0 /        (R)4100.0,
              (R)51.0 /          (R)82.0,
              (R)33.0 /         (R)164.0,
              (R)12.0 /          (R)41.0,
                 ZERO,
                  ONE,

                  /* ---- 21st embedding formula: method RKV8(7)  --- */
/* A~ */
              (R)31.0 /         (R)720.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
              (R)16.0 /          (R)75.0,
           (R)16807.0 /       (R)79200.0,
           (R)16807.0 /       (R)79200.0,
             (R)243.0 /        (R)1760.0,
                 ZERO,
                 ZERO,
             (R)243.0 /        (R)1760.0,
              (R)31.0 /         (R)720.0,
/* A */
              (R)13.0 /         (R)288.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
              (R)32.0 /         (R)125.0,
           (R)31213.0 /      (R)144000.0,
            (R)2401.0 /       (R)12375.0,
            (R)1701.0 /       (R)14080.0,
            (R)2401.0 /       (R)19200.0,
              (R)19.0 /         (R)450.0,
                 ZERO,
                 ZERO,
/* a */
                  ONE /             FOUR,
                  ONE /          (R)12.0,
                  ONE /            EIGHT,
                  TWO /             FIVE,
                 HALF,
                  SIX /           (R)7.0,
                  ONE /           (R)7.0,
                  TWO /            THREE,
                  TWO /           (R)7.0,
                  ONE,
                  ONE /            THREE,
                  ONE,
/* b */
                  ONE /             FOUR,
                 FIVE /          (R)72.0,
                  ONE /          (R)72.0,
                  ONE /          (R)32.0,
                 ZERO,
                THREE /          (R)32.0,
             (R)106.0 /         (R)125.0,
                 ZERO,
            (R)-408.0 /         (R)125.0,
             (R)352.0 /         (R)125.0,
                  ONE /          (R)48.0,
                 ZERO,
                 ZERO,
                EIGHT /          (R)33.0,
             (R)125.0 /         (R)528.0,
           (R)-1263.0 /        (R)2401.0,
                 ZERO,
                 ZERO,
           (R)39936.0 /       (R)26411.0,
          (R)-64125.0 /       (R)26411.0,
            (R)5520.0 /        (R)2401.0,
              (R)37.0 /         (R)392.0,
                 ZERO,
                 ZERO,
                 ZERO,
            (R)1625.0 /        (R)9408.0,
                 -TWO /          (R)15.0,
              (R)61.0 /        (R)6720.0,
           (R)17176.0 /       (R)25515.0,
                 ZERO,
                 ZERO,
          (R)-47104.0 /       (R)25515.0,
            (R)1325.0 /         (R)504.0,
          (R)-41792.0 /       (R)25515.0,
           (R)20237.0 /      (R)145800.0,
            (R)4312.0 /        (R)6075.0,
          (R)-23834.0 /      (R)180075.0,
                 ZERO,
                 ZERO,
          (R)-77824.0 /     (R)1980825.0,
         (R)-636635.0 /      (R)633864.0,
          (R)254048.0 /      (R)300125.0,
            (R)-183.0 /        (R)7000.0,
                EIGHT /          (R)11.0,
            (R)-324.0 /        (R)3773.0,
           (R)12733.0 /        (R)7600.0,
                 ZERO,
                 ZERO,
          (R)-20032.0 /        (R)5225.0,
          (R)456485.0 /       (R)80256.0,
          (R)-42599.0 /        (R)7125.0,
          (R)339227.0 /      (R)912000.0,
           (R)-1029.0 /        (R)4180.0,
            (R)1701.0 /        (R)1408.0,
            (R)5145.0 /        (R)2432.0,
          (R)-27061.0 /      (R)204120.0,
                 ZERO,
                 ZERO,
           (R)40448.0 /      (R)280665.0,
        (R)-1353775.0 /     (R)1197504.0,
           (R)17662.0 /       (R)25515.0,
          (R)-71687.0 /     (R)1166400.0,
              (R)98.0 /         (R)225.0,
                  ONE /          (R)16.0,
            (R)3773.0 /       (R)11664.0,
                 ZERO,
           (R)11203.0 /        (R)8680.0,
                 ZERO,
                 ZERO,
          (R)-38144.0 /       (R)11935.0,
         (R)2354425.0 /      (R)458304.0,
          (R)-84046.0 /       (R)16275.0,
          (R)673309.0 /     (R)1636800.0,
            (R)4704.0 /        (R)8525.0,
            (R)9477.0 /       (R)10912.0,
           (R)-1029.0 /         (R)992.0,
                 ZERO,
             (R)729.0 /         (R)341.0,

                  /* ---- 22nd embedding formula: method RKV9(8)  --- */
/* A~ */
              (R)23.0 /         (R)525.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
             (R)171.0 /        (R)1400.0,
              (R)86.0 /         (R)525.0,
              (R)93.0 /         (R)280.0,
           (R)-2048.0 /        (R)6825.0,
               -THREE /       (R)18200.0,
              (R)39.0 /         (R)175.0,
                 ZERO,
                 NINE /          (R)25.0,
             (R)233.0 /        (R)4200.0,
/* A */
             (R)103.0 /        (R)1680.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
             (R)-27.0 /         (R)140.0,
              (R)76.0 /         (R)105.0,
            (R)-201.0 /         (R)280.0,
            (R)1024.0 /        (R)1365.0,
                THREE /        (R)7280.0,
              (R)12.0 /          (R)35.0,
                 NINE /         (R)280.0,
                 ZERO,
                 ZERO,
/* a */
                  ONE /          (R)12.0,
                  ONE /             NINE,
                  ONE /              SIX,
                 ZERO,
                 ZERO,
                 ZERO,
                  TWO /            THREE,
                 HALF,
                  ONE /            THREE,
                  ONE /             FOUR,
                 FOUR /            THREE,
                 FIVE /              SIX,
                  ONE,
                  ONE /              SIX,
                  ONE,
/* b */
                  ONE /          (R)12.0,
                  ONE /          (R)27.0,
                  TWO /          (R)27.0,
                  ONE /          (R)24.0,
                 ZERO,
                  ONE /            EIGHT,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                  TWO /          (R)27.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
              (R)19.0 /         (R)256.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                -NINE /         (R)256.0,
              (R)11.0 /         (R)144.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 -ONE /          (R)16.0,
               -EIGHT /          (R)27.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
          (R)-26624.0 /          (R)81.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
            (R)2624.0 /        (R)1053.0,
                THREE /        (R)1664.0,
            (R)-137.0 /        (R)1296.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
            (R)-299.0 /          (R)48.0,
             (R)184.0 /          (R)81.0,
             (R)-44.0 /             NINE,
           (R)-5120.0 /        (R)1053.0,
             (R)-11.0 /         (R)468.0,
              (R)16.0 /             NINE,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
             (R)320.0 /         (R)567.0,
                 -ONE /        (R)1920.0,
                 FOUR /         (R)105.0,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
         (R)-169984.0 /        (R)9087.0,
             (R)-87.0 /       (R)30290.0,
             (R)492.0 /        (R)1165.0,
                 ZERO,
            (R)1260.0 /         (R)233.0
  };


  /* Note all level numbers, all k1_km flags, all global error orders */
  /* of low order and all addresses of their coefficients for the     */
  /* 0, 1, ..., NE - 1 embedding formulas.                            */

  /* These makros are designed to shorten the calling sequences of    */
  /* addresses.                                                       */

  #define m0    3            /* assign each level number a name       */
  #define m1    5
  #define m2    6
  #define m3    6
  #define m4    6
  #define m5    7
  #define m6    7
  #define m7    7
  #define m8    7
  #define m9    8
  #define m10   8
  #define m11   8
  #define m12   8
  #define m13   8
  #define m14   8
  #define m15   9
  #define m16   9
  #define m17   9
  #define m18  10
  #define m19  13
  #define m20  13
  #define m21  13
  #define m22  16

  #define O0             0   /* define an index Oi for each formula i;*/
  #define O1   nO(O0,  m0)   /* define `coefficients' for the vector, */
  #define O2   nO(O1,  m1)   /* so that Oi denotes the spot in        */
  #define O3   nO(O2,  m2)   /* `coefficients', where the coefficients*/
  #define O4   nO(O3,  m3)   /* A~, A, a, b of formula i start:       */
  #define O5   nO(O4,  m4)
  #define O6   nO(O5,  m5)   /* `coefficients[O3]  ==A~[0]'           */
  #define O7   nO(O6,  m6)   /* `coefficients[O3+1]==A~[1]', etc.     */
  #define O8   nO(O7,  m7)
  #define O9   nO(O8,  m8)
  #define O10  nO(O9,  m9)
  #define O11  nO(O10, m10)
  #define O12  nO(O11, m11)
  #define O13  nO(O12, m12)
  #define O14  nO(O13, m13)
  #define O15  nO(O14, m14)
  #define O16  nO(O15, m15)
  #define O17  nO(O16, m16)
  #define O18  nO(O17, m17)
  #define O19  nO(O18, m18)
  #define O20  nO(O19, m19)
  #define O21  nO(O20, m20)
  #define O22  nO(O21, m21)

  /* In order to get from one O value to the next (nO denotes next    */
  /* Offset), we need to realize that a method of level m has         */
  /* (m+m+m-1)+(1+2+...+m-1) = m*(m+5)/2-1  coeficients.              */
  #define nO(VorigerOffset,m)  (VorigerOffset + m * (m + 5) / 2 - 1)

  /* find the coefficient addresses depending on the m and O value    */
  #define MachKoeff(Offset,m)                        \
      /* A~= */  koeffizienten + Offset,             \
      /* A = */  koeffizienten + Offset +     m,     \
      /* a = */  koeffizienten + Offset + 2 * m,     \
      /* b = */  koeffizienten + Offset + 3 * m - 1

  static koefftyp formel[NE] =
  {
    { m0,  FALSE, TWO,    0, MachKoeff(O0, m0)       /*  0: RK3(2)    */
    },
    { m1,  TRUE,  THREE,  0, MachKoeff(O1, m1)       /*  1: RKF4(3)   */
    },
    { m2,  FALSE, FOUR,   0, MachKoeff(O2, m2)       /*  2: RKF5(4)   */
    },
    { m3,  FALSE, FOUR,   0, MachKoeff(O3, m3)       /*  3: RK5(4)6M  */
    },
    { m4,  FALSE, FOUR,   0, MachKoeff(O4, m4)       /*  4: RKE5(4)   */
    },
    { m5,  TRUE,  FOUR,   0, MachKoeff(O5, m5)       /*  5: HIHA5     */
    },
    { m6,  TRUE,  FOUR,   0, MachKoeff(O6, m6)       /*  6: RK5(4)7S  */
    },
    { m7,  TRUE,  FOUR,   0, MachKoeff(O7, m7)       /*  7: RK5(4)7M  */
    },
    { m8,  TRUE,  FOUR,   0, MachKoeff(O8, m8)       /*  8: RK5(4)7C  */
    },
    { m9,  FALSE, FIVE,   0, MachKoeff(O9, m9)       /*  9: RK6(5)8M  */
    },
    { m10, FALSE, FIVE,   0, MachKoeff(O10, m10)     /* 10: RK6(5)8S  */
    },
    { m11, FALSE, FIVE,   0, MachKoeff(O11, m11)     /* 11: RK6(5)8C  */
    },
    { m12, FALSE, FIVE,   0, MachKoeff(O12, m12)     /* 12: RKV6(5)   */
    },
    { m13, FALSE, FIVE,   0, MachKoeff(O13, m13)     /* 13: RKF6(5)A  */
    },
    { m13, FALSE, FIVE,   0, MachKoeff(O14, m14)     /* 14: RKF6(5)B  */
    },
    { m15, TRUE,  FIVE,   0, MachKoeff(O15, m15)     /* 15: RKC6(5)   */
    },
    { m16, TRUE,  FIVE,   0, MachKoeff(O16, m16)     /* 16: RKV6(5)9A */
    },
    { m17, TRUE,  FIVE,   0, MachKoeff(O17, m17)     /* 17: RKV6(5)9B */
    },
    { m18, FALSE, SIX,    0, MachKoeff(O18, m18)     /* 18: RKV7(6)   */
    },
    { m19, FALSE, (R)7.0, 0, MachKoeff(O19, m19)     /* 19: RK8(7)13M */
    },
    { m20, FALSE, (R)7.0, 0, MachKoeff(O20, m20)     /* 20: RKF8(7)   */
    },
    { m21, FALSE, (R)7.0, 0, MachKoeff(O21, m21)     /* 21: RKV8(7)   */
    },
    { m22, FALSE, EIGHT,  0, MachKoeff(O22, m22)     /* 22: RKV9(8)   */
    }
  };
  #undef m0
  #undef m1
  #undef m2
  #undef m3
  #undef m4
  #undef m5
  #undef m6
  #undef m7
  #undef m8
  #undef m9
  #undef m10
  #undef m11
  #undef m12
  #undef m13
  #undef m14
  #undef m15
  #undef m16
  #undef m17
  #undef m18
  #undef m19
  #undef m20
  #undef m21
  #undef m22
  #undef O0
  #undef O1
  #undef O2
  #undef O3
  #undef O4
  #undef O5
  #undef O6
  #undef O7
  #undef O8
  #undef O9
  #undef O10
  #undef O11
  #undef O12
  #undef O13
  #undef O14
  #undef O15
  #undef O16
  #undef O17
  #undef O18
  #undef O19
  #undef O20
  #undef O21
  #undef O22
  #undef nO
  #undef MachKoeff


  *koeff            = formel[neinb];
  koeff->maxschritt = fmax / koeff->stufenzahl;

  #define As  koeff->As
  #define A   koeff->A
  #define a   koeff->a
  #define b   koeff->b
  #define m   koeff->stufenzahl

  if (neinb == 16 &&       /* special case embedding formula 16?      */
      ! ok16)              /* correction not made?                    */
  {                        /* correct the following coefficients by   */
    ok16   = TRUE;         /* the factor h:                           */
    h      = SQRT(TEN);
    As[0]  =               (R)31.0      /           (R)324.0 -
                           (R)37.0      /          (R)4860.0 * h;
    As[3]  = (          (R)37435.0 -
                         (R)3235.0 * h) /         (R)69228.0;
    As[4]  =         (R)-1245184.0      /       (R)1090341.0 +
                      (R)9699328.0      /      (R)16355115.0 * h;
    As[5]  =               (R)71.0      /            (R)54.0 -
                           (R)74.0      /           (R)135.0 * h;
    As[6]  =              (R)625.0      /           (R)486.0 -
                          (R)250.0      /           (R)729.0 * h;
    As[7]  =              (R)-23.0      /            (R)21.0 +
                           (R)37.0      /           (R)105.0 * h;
    A[0]   =                  FIVE      /            (R)54.0 -
                               TWO      /           (R)135.0 * h;
    A[3]   = (           (R)2390.0 +
                         (R)2290.0 * h) /         (R)17307.0;
    A[4]   =            (R)40960.0      /        (R)121149.0 +
                       (R)262144.0      /        (R)605745.0 * h;
    A[5]   =                   TWO      /            (R)27.0 -
                           (R)64.0      /           (R)135.0 * h;
    A[7]   =           (R)150029.0      /        (R)443709.0 -
                       (R)236267.0      /       (R)2218545.0 * h;
    A[8]   =             (R)2411.0      /        (R)126774.0 +
                         (R)1921.0      /         (R)63387.0 * h;
    a[1]   =                  FOUR      /               NINE -
                              FOUR      /            (R)45.0 * h;
    a[2]   =                   TWO      /              THREE -
                               TWO      /            (R)15.0 * h;
    b[1]   = (           (R)-268.0 +
                           (R)92.0 * h) /           (R)405.0;
    b[2]   = (            (R)448.0 -
                          (R)128.0 * h) /           (R)405.0;
    b[3]   =                   ONE      /                SIX -
                               ONE      /            (R)30.0 * h;
    b[5]   =                   ONE      /                TWO -
                               ONE      /                TEN * h;
    b[6]   =            (R)11547.0      /         (R)32768.0 +
                          (R)405.0      /         (R)16384.0 * h;
    b[8]   =           (R)-18225.0      /         (R)32768.0 -
                         (R)5103.0      /         (R)16384.0 * h;
    b[9]   =            (R)12555.0      /         (R)16384.0 +
                         (R)2349.0      /          (R)8192.0 * h;
    b[10]  =         (R)19662371.0      /      (R)51149376.0 +
                       (R)441281.0      /      (R)12787344.0 * h;
    b[12]  =         (R)-3786045.0      /       (R)5683264.0 -
                       (R)252663.0      /        (R)710408.0 * h;
    b[13]  =       (R)1570556745.0      /    (R)1821486112.0 +
                    (R)290041461.0      /     (R)910743056.0 * h;
    b[14]  = (      (R)-41227072.0 +
                      (R)1374464.0 * h) /     (R)512292969.0;
    b[15]  =       (R)-154207593.0      /     (R)369412160.0 -
                   (R)1829424339.0      /   (R)11544130000.0 * h;
    b[17]  =       (R)2659895739.0      /    (R)1847060800.0 +
                    (R)653855409.0      /    (R)1154413000.0 * h;
    b[18]  =    -(R)349492176711.0      /  (R)591982986400.0 -
                 (R)359784638379.0      / (R)1479957466000.0 * h;
    b[19]  =     (R)153920585664.0      /   (R)92497341625.0 +
                 (R)311066673408.0      /  (R)462486708125.0 * h;
    b[20]  =            (R)-1944.0      /          (R)1625.0 -
                         (R)6804.0      /          (R)8125.0 * h;
    b[21]  =      (R)70594945601.0      /   (R)21406013856.0 +
                  (R)21473424323.0      /   (R)21406013856.0 * h;
    b[23]  =       (R)-794525145.0      /      (R)88090592.0 -
                    (R)249156075.0      /      (R)88090592.0 *  h;
    b[24]  =     (R)866290968775.0      /  (R)254097312624.0 +
                 (R)256998959765.0      /  (R)254097312240.0 * h;
    b[25]  = ((R)-15964196472448.0 -
                (R)5039429245312.0 * h) / (R)1286367645159.0;
    b[26]  = (          (R)17017.0 +
                         (R)5075.0 * h) /          (R)1116.0;
    b[27]  = (          (R)42875.0 +
                        (R)16625.0 * h) /         (R)90396.0;
    b[28]  =               (R)31.0      /           (R)324.0 -
                           (R)37.0      /          (R)4860.0 * h;
    b[31]  = (          (R)37435.0 -
                         (R)3235.0 * h) /         (R)69228.0;
    b[32]  =         (R)-1245184.0      /       (R)1090341.0 +
                      (R)9699328.0      /      (R)16355115.0 * h;
    b[33]  =               (R)71.0      /            (R)54.0 -
                           (R)74.0      /           (R)135.0 * h;
    b[34]  =              (R)625.0      /           (R)486.0 -
                          (R)250.0      /           (R)729.0 * h;
    b[35]  =              (R)-23.0      /            (R)21.0 +
                           (R)37.0      /           (R)105.0 * h;
  }

  if (neinb == 22 &&       /* special case embedding formula 22?      */
      ! ok22)              /* correction not made?                    */
  {                        /* correct the following coefficients using*/
    ok22   = TRUE;         /* the factor h:                           */
    h      = SQRT(SIX);
    a[3]   = (         TWO + h *         TWO) /     (R)15.0;
    a[4]   = (         SIX + h              ) /     (R)15.0;
    a[5]   = (         SIX - h              ) /     (R)15.0;
    b[6]   = (        FOUR + h *     (R)94.0) /    (R)375.0;
    b[8]   = (    (R)-94.0 - h *     (R)84.0) /    (R)125.0;
    b[9]   = (    (R)328.0 + h *    (R)208.0) /    (R)375.0;
    b[10]  = (        NINE - h              ) /    (R)150.0;
    b[13]  = (    (R)312.0 + h *     (R)32.0) /   (R)1425.0;
    b[14]  = (     (R)69.0 + h *     (R)29.0) /    (R)570.0;
    b[15]  = (    (R)927.0 - h *    (R)347.0) /   (R)1250.0;
    b[18]  = ( (R)-16248.0 + h *   (R)7328.0) /   (R)9375.0;
    b[19]  = (   (R)-489.0 + h *    (R)179.0) /   (R)3750.0;
    b[20]  = (  (R)14268.0 - h *   (R)5798.0) /   (R)9375.0;
    b[26]  = (     (R)16.0 - h              ) /     (R)54.0;
    b[27]  = (     (R)16.0 + h              ) /     (R)54.0;
    b[33]  = (    (R)118.0 - h *     (R)23.0) /    (R)512.0;
    b[34]  = (    (R)118.0 + h *     (R)23.0) /    (R)512.0;
    b[41]  = (    (R)266.0 - h              ) /    (R)864.0;
    b[42]  = (    (R)266.0 + h              ) /    (R)864.0;
    b[45]  = (   (R)5034.0 - h *    (R)271.0) /  (R)61440.0;
    b[51]  = (   (R)7859.0 - h *   (R)1626.0) /  (R)10240.0;
    b[52]  = (  (R)-2232.0 + h *    (R)813.0) /  (R)20480.0;
    b[53]  = (   (R)-594.0 + h *    (R)271.0) /    (R)960.0;
    b[54]  = (    (R)657.0 - h *    (R)813.0) /   (R)5120.0;
    b[55]  = (   (R)5996.0 - h *   (R)3794.0) /    (R)405.0;
    b[60]  = (  (R)-4342.0 - h *    (R)338.0) /        NINE;
    b[61]  = ( (R)154922.0 - h *  (R)40458.0) /    (R)135.0;
    b[62]  = (  (R)-4176.0 + h *   (R)3794.0) /     (R)45.0;
    b[63]  = ((R)-340864.0 + h * (R)242816.0) /    (R)405.0;
    b[64]  = (  (R)26304.0 - h *  (R)15176.0) /     (R)45.0;
    b[66]  = (   (R)3793.0 + h *   (R)2168.0) / (R)103680.0;
    b[71]  = (   (R)4042.0 + h *   (R)2263.0) /  (R)13824.0;
    b[72]  = ((R)-231278.0 + h *  (R)40717.0) /  (R)69120.0;
    b[73]  = (   (R)7947.0 - h *   (R)2168.0) /  (R)11520.0;
    b[74]  = (   (R)1048.0 - h *    (R)542.0) /    (R)405.0;
    b[75]  = (  (R)-1383.0 + h *    (R)542.0) /    (R)720.0;
    b[83]  = (   (R)5642.0 - h *    (R)337.0) /    (R)864.0;
    b[84]  = (   (R)5642.0 + h *    (R)337.0) /    (R)864.0;
    b[91]  = (  (R)33617.0 - h *   (R)2168.0) / (R)518400.0;
    b[96]  = (  (R)-3846.0 + h *     (R)31.0) /  (R)13824.0;
    b[97]  = ( (R)155338.0 - h *  (R)52807.0) / (R)345600.0;
    b[98]  = ( (R)-12537.0 + h *   (R)2168.0) /  (R)57600.0;
    b[99]  = (     (R)92.0 + h *    (R)542.0) /   (R)2025.0;
    b[100] = (  (R)-1797.0 - h *    (R)542.0) /   (R)3600.0;
    b[105] = ( (R)-36487.0 - h *  (R)30352.0) / (R)279600.0;
    b[110] = ( (R)-29666.0 - h *   (R)4499.0) /   (R)7456.0;
    b[111] = ((R)2779182.0 - h * (R)615973.0) / (R)186400.0;
    b[112] = ( (R)-94329.0 + h *  (R)91056.0) /  (R)93200.0;
    b[113] = ((R)-232192.0 + h * (R)121408.0) /  (R)17475.0;
    b[114] = ( (R)101226.0 - h *  (R)22764.0) /   (R)5825.0;
  }
#ifdef DEBUG

  {
  int i;
  printf("Coefficients of the embedding formula %2d:\n", neinb);
  printf("A~ = \n");
  for (i = 0; i < m; i++)
    printf("    %15.8"LZP"E\n", As[i]);
  printf("A = \n");
  for (i = 0; i < m; i++)
    printf("    %15.8"LZP"E\n", A[i]);
  printf("a = \n");
  for (i = 0; i < m - 1; i++)
    printf("    %15.8"LZP"E\n", a[i]);
  printf("b = \n");
  for (i = 0; i < m - 1; i++)
  {
    int j;
    printf("  ");
    for (j = 0; j <= i; j++)
      printf("  %15.8"LZP"E", b[i * (i + 1) / 2 + j]);
    printf("\n");
  }
  }
#endif
  #undef As
  #undef A
  #undef a
  #undef b
  #undef m
}



/* ------------------------------------------------------------------ */

static int rk_schritt     /* perform one Runge-Kutta step ............*/
/*.IX{rk\unt schritt}*/
                     (
                      REAL      x,         /* initial x-value ........*/
                      REAL      h,         /* Step size ..............*/
                      REAL      y[],       /* y-value at x ...........*/
                      int       n,         /* number of DE system ....*/
                      REAL      *k[],      /* Runge-Kutta vectors ki .*/
                      dglsysfnk dgl,       /* right hand side for DE .*/
                      REAL      yhigh[],   /* solution at x+h, h. ord.*/
                      REAL      ylow[],    /* solution at x+h, l. ord.*/
                      koefftyp  *koeff,    /* data foe method ........*/
                      int       newstp,    /* new step ? .............*/
                      int       firstp,    /* first step? ............*/
                      REAL      xzi,       /* overflow bound .........*/
                      REAL      yhilf[]    /* aux vector .............*/
                     )                     /* error code .............*/


/***********************************************************************
* Perform one integration step using the Runge-Kutta embedding formula *
* determined in koeff. This produces two approximate solutions at x+h, *
* one of higher order, one of lower error order.                       *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x       initial x-value                                              *
* h       step size                                                    *
* y       [0..n-1] y-vector of the initial solution at x               *
* n       number of DEs                                                *
* k       only for cleared newstp: [0..m-1,0..n-1] matrix, whose first *
*                                  row contains the k1 vector of the   *
*                                  previous try, where m = level       *
*                                  number of method.                   *
* dgl     right hand side of the system of DEs                         *
* koeff   Structure, which desribes all relevant parameters and        *
*         coefficients of the chosen Runge-Kutta embedding formula     *
* newstp  determines whether a new step should be performed (set) or   *
*         a step ought to be repeated with a smaller step size         *
*         (cleared). For a repeated step there is no need to           *
*         recompute k1.                                                *
* firstp  indicates that this is the first step; hence the fact that   *
*         "k1=km" cannot be used in the embedding formula.             *
* xzi     overflow bound                                               *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* k       [0..m-1,0..n-1] matrix, with the vectors ki of the           *
*         integration as rows. (m = level number)                      *
* yhigh   [0..n-1] vector with y-values of higher order approximation  *
* ylow    [0..n-1] vector, ditto for low order                         *
* yhilf   [0..n-1] aux vector for k                                    *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
*  = 0: all ok                                                         *
*  = 1: overflow                                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, dglsysfnk, koefftyp, norm_max, ZERO                            *
***********************************************************************/

{
  REAL *b;              /* pointer to one row of the matrix of        */
                        /* b-coefficients of the embedding formula    */
  REAL xhilf;           /* aux variable                               */
  REAL summe;           /* Sum used for  ki vector                    */
  REAL summe1;          /* Sum used to compute high order             */
                        /* approximation of y(x+h)                    */
  REAL summe2;          /* ditto for low order one                    */
  int  m;               /* level number of embedding formula          */
  int  i;               /* loop counters                              */
  int  j;
  int  l;


  m = koeff->stufenzahl;

  if (newstp)                                /* new step?             */
  {
    if (! firstp &&                          /* not the first step?   */
        koeff->k1_km)                        /* new k1 = old km?      */
      copy_vector(k[0], k[m - 1], n);        /* copy                  */
    else                                     /* not a k1_km formula   */
    {                                        /* or first step?        */
      dgl(x, y, k[0]);                       /* compute new k1        */
      if (norm_max(k[0], n) > xzi)           /* danger of overflow?   */
        return 1;                            /* report error          */
    }
  }

  /* ---- compute remaining ki vectors and watch for overflow    ---- */

  for (b = koeff->b, i = 1; i < m; b += i, i++)
  {
    xhilf = x + koeff->a[i - 1] * h;
    for (l = 0; l < n; l++)
    {
      for (summe = ZERO, j = 0; j < i; j++)
        summe += b[j] * k[j][l];
      yhilf[l] = y[l] + summe * h;
    }
    dgl(xhilf, yhilf, k[i]);
    if (norm_max(k[i], n) > xzi)
      return 1;
  }

  for (l = 0; l < n; l++)   /* compute both approximations for y(x+h) */
  {
    for (summe1 = ZERO, j = 0; j < m; j++)
      summe1 += koeff->As[j] * k[j][l];
    for (summe2 = ZERO, j = 0; j < m; j++)
      summe2 += koeff->A[j] * k[j][l];

    yhigh[l] = y[l] + summe1 * h;         /* high order approximation */
    ylow[l]  = y[l] + summe2 * h;         /* low order approximation  */
  }


  return 0;
}



/* ------------------------------------------------------------------ */

static int awpl               /* use desired embedding formula .......*/
/*.IX{awpl}*/
               (
                REAL      *x,          /* starting/final x-value .....*/
                REAL      *h,          /* starting/final step size ...*/
                REAL      beta,        /* desired final value ........*/
                REAL      abserr,      /* absolute error bound .......*/
                REAL      relerr,      /* relative error bound .......*/
                int       n,           /* number of DEs in system ....*/
                dglsysfnk dgl,         /* right hand side of DE system*/
                REAL      y[],         /* solution at  x .............*/
                int       hullstp,     /* step size control according */
                                       /* to Hull? ...................*/
                REAL      eps,         /* 100 * machine constant .....*/
                REAL      xzi,         /* overflow bound .............*/
                koefftyp  *koeff,      /* Data for formula ...........*/
                long      *aufrufe     /* actual # of calls of dgl() .*/
               )                       /* error code .................*/

/***********************************************************************
* Solve a system of first order ordinary differential equations using a*
* Runge-Kutta embedding formula for the interval [x0,beta]. The step   *
* size is adjusted following Hull or general initial value solvers with*
* automatic step size control.                                         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        initial x-value                                             *
* h        initial step size                                           *
* beta     desired final x-value                                       *
* abserr   absolute error bound (abserr >= 0). If abserr = 0  we check *
*          relative error only.                                        *
* relerr   relative error bound (relerr >= 0). If relerr = 0  we check *
*          the absolute error only.                                    *
* n        number of DEs                                               *
* dgl      right hand side of DE system                                *
* y        [0..n-1] vevtor with initial y-values at x                  *
* hullstp  step size control according to Hull or according to general *
*          initial value solvers                                       *
* eps      100 * machine constant                                      *
* xzi      overflow bound                                              *
* koeff    Structure with hte essential data and coefficients for the  *
*          chosen Runge-Kutta embedding formula                        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value (normally x = beta)                           *
* h        final step size                                             *
* y        [0..n-1] vector with the solution at x                      *
* aufrufe  counter for # of calls of dgl()                             *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
*   =   0: all ok                                                      *
*   =  -1: desired relative accuracy is below 100 * machine            *
*          constant in parts of the interval of integration.           *
*          In these regions we compute with 100 * machine constant     *
*          as absolute error bound.                                    *
*   =  -2: number of allowed evaluations has been reached              *
*   =  -3: possible overflow                                           *
*   =  -4: optimal step size became too small                          *
*   =  -5: lack of memory                                              *
*                                                                      *
* global names used:                                                   *
* REAL, dglsysfnk, koefftyp, vminit, vmalloc, MATRIX, VEKTOR,          *
* vmcomplete, vmfree, FALSE, ONE, rk_schritt, norm_max, TWO, POW,      *
* FABS, TRUE, copy_vector, max, ZERO                                   *
***********************************************************************/

{
  REAL *yhigh;          /* [0..n-1] high order y-value approximation  */
  REAL *ylow;           /* ditto for low order                        */
  REAL **k;             /* [0..m-1,0..n-1] matrix (see rk_schritt())  */
  REAL *ydiff;          /* [0..n-1] aux vector                        */
  REAL *yhilf;          /* [0..n-1] aux vector for rk_schritt()       */
  REAL temp;            /* aux storage for latest optimal step size   */
  REAL expo;            /* aux variable                               */
  REAL delta;           /* estimate for local error                   */
  REAL epslon;          /* Tolerance for local error                  */
  REAL xhilf;           /* aux variables                              */
  REAL s = ZERO;
  REAL xend;            /* aux variable for interval end point        */
  REAL qg;              /* global error order for low order method    */
                        /* (+ 1 for Hull step sizes)                  */
  int  newstp;          /* see rk_schritt()                           */
  int  m;               /* level number of embedding formula          */
  int  laststp;         /* Flag indicating the curent step will reach */
                        /* the end point of the interval              */
  int  firstp;          /* Flag for the first step                    */
  long step;            /* Step counter                               */
  int  j;               /* loop counter                               */
  int  rkerr;           /* error code from rk_schritt()               */
  int  fehler;          /* error code                                 */
  void *vmblock;        /* List of dynamically allocated vectors and  */
                        /* matrices                                   */


  /* ------------ allocate dynamic vectors and matrices ------------- */

  vmblock = vminit();
  k       = (REAL **)vmalloc(vmblock, MATRIX, koeff->stufenzahl, n);
  yhigh   = (REAL *) vmalloc(vmblock, VEKTOR, n,                 0);
  ylow    = (REAL *) vmalloc(vmblock, VEKTOR, n,                 0);
  ydiff   = (REAL *) vmalloc(vmblock, VEKTOR, n,                 0);
  yhilf   = (REAL *) vmalloc(vmblock, VEKTOR, n,                 0);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return -5;
  }


  fehler  = -2;                          /* prepare iteration     */
  firstp  = TRUE;
  newstp  = TRUE;
  laststp = FALSE;

  if (hullstp)
    qg = koeff->qg + ONE;
  else
    qg = koeff->qg;
  expo = ONE / qg;
  temp = *h;
  m    = koeff->stufenzahl;


  for (step = koeff->maxschritt;/* execute maximally koeff->maxschritt*/
       step != 0;               /* Runge-Kutta steps                  */
       step--,
       firstp = TRUE
      )
  {

    rkerr = rk_schritt(*x, *h, y, n, k, dgl,   /* one integration step*/
                       yhigh, ylow, koeff,
                       newstp, firstp, xzi,
                       yhilf);

    *aufrufe += m;                       /* count calls of dgl()      */
    if (! newstp)
      --*aufrufe;

    if (rkerr || norm_max(yhigh, n) > xzi ||            /* overflow?  */
                 norm_max(ylow,  n) > xzi)
    {
      fehler = -3;
      break;
    }

                                 /* ---- determine new step size  --- */

    for (j = 0; j < n; j++)
      ydiff[j] = yhigh[j] - ylow[j];
    delta = norm_max(ydiff, n);
    if (! hullstp && delta < eps)
      s = TWO;
    else
    {
      epslon = abserr + relerr * norm_max(yhigh, n);
      if (epslon < eps)
        epslon = eps,
        fehler = -1;
      if (hullstp)
        delta /= epslon;
      else
        s = POW(FABS(*h) * epslon / delta, expo);
    }

    if ((! hullstp && s     < ONE)  ||            /* s too small or   */
        (  hullstp && delta > ONE))               /* delta too large? */
    {
      newstp = FALSE;

      if (! hullstp)                            /* find new step size */
        *h *= max(HALF, s);
      else
        delta =  min(POW((REAL)3.6, qg),
                     delta),
        *h    *= (REAL)0.9 / POW(delta, expo);
      if (FABS(*h) < eps)                     /* step size too small? */
      {
        fehler = -4;                             /* report error      */
        break;
      }

      if (laststp)          /* If this would have been the last call, */
        laststp = FALSE;    /* the repeated call would not be so!     */

    }
    else                    /* s large enough or delta small enough?  */
    {
      *x += *h;                                 /* reassign x and y   */
      copy_vector(y, yhigh, n);

      if (laststp)                /* previous call was the final one? */
      {
        *h     = temp;                        /* store new step size  */
        fehler = 0;                           /* report success       */
        break;                                /* stop iterations      */
      }

      else                         /* prepare a new integration step? */
      {
        newstp = TRUE;

        if (! hullstp)            /* compute new step size to be at   */
          *h *= min(TWO ,s);      /* most 2 or four times old one     */
        else
          delta =  max(POW((REAL)0.225, qg), delta),
          *h    *= (REAL)0.9 / POW(delta,  expo);

        xhilf = *x + *h;               /* Next step beyond end point? */
        xend  = beta - (REAL)0.1 * *h;
        if ((*h < ZERO && xhilf < xend) ||
            (*h > ZERO && xhilf > xend)
           )
        {
          laststp = TRUE;             /* set flag for final step      */
          temp    = *h;               /* record new step size in temp */
          *h      = beta - *x;        /* reduce step size as needed   */
        }
      }
    }
  }

  vmfree(vmblock);
  return fehler;
}


/* ------------------------------------------------------------------ */

static int hstart       /* compute starting step size for an IVP .....*/
/*.IX{hstart}*/
                 (
                  dglsysfnk dgl,        /* right hand side of system .*/
                  int       n,          /* number of DEs .............*/
                  REAL      x,          /* Starting x-value ..........*/
                  REAL      beta,       /* desired final x-value .....*/
                  REAL      y[],        /* solution of DE system at x */
                  REAL      relerr,     /* relative error bound ......*/
                  REAL      abserr,     /* absolute error bound ......*/
                  REAL      qg,         /* global error bound for meth*/
                  REAL      *h,         /* step size .................*/
                  long      *aufrufe    /* max. # of calls of dgl() ..*/
                 )                      /* error code ................*/

/***********************************************************************
* Compute the starting step size for an initial value problem.         *
* Here we use the Lipschitz constant and the upper bound of the first  *
* and second derivatives of the differential equation in a neighborhood*
* of the initial x-value x0. The algorithm is derived from the  DEPAC  *
* package (design of a user oriented package of ODE solvers).          *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* dgl      right hand side of DE system                                *
* n        number of DEs                                               *
* x        initial x-value                                             *
* beta     desired final x-value                                       *
* y        [0..n-1] initial y-vector at x                              *
* relerr   relative error bound (relerr >= 0). If relerr = 0 we only   *
*          the absolute accuracy                                       *
* abserr   absolute error bound (abserr >= 0). If abserr = 0 we only   *
*          test for the relative error.                                *
* qg       global error order of methode                               *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* h        starting step size                                          *
* aufrufe  number of calls of dgl()                                    *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
*   =  0: all ok                                                       *
*   = -3: lack of memory                                               *
*                                                                      *
* global names used:                                                   *
* =================                                                    *
* dglsysfnk, REAL, vminit, vmalloc, MATRIX, vmcomplete, vmfree, FABS,  *
* POW, SIGN, max, norm_max, ZERO, ONE, LOG, TEN, HALF, SQRT, MACH_EPS, *
* POSMAX                                                               *
***********************************************************************/

{
  REAL **hilf;            /* [0..4,0..n-1] matrix                     */
  REAL dx;                /* length of interval                       */
  REAL absdx;             /* magnitude of dx                          */
  REAL relper;            /* 0.375 * dsmall                           */
  REAL da;                /* variation in x                           */
  REAL delf;              /* aux variable                             */
  REAL dfdxb;             /* upper bound for second derivative found  */
                          /* from difference quotient                 */
  REAL fbnd;              /* upper bound for first derivative         */
  REAL dely;              /* aux variable                             */
  REAL dfdub;             /* Lipschitz constant                       */
  REAL dy;                /* aux variable                             */
  REAL ydpb;              /* upper bound for second derivative        */
  REAL tolmin;            /* aux variables for tolp                   */
  REAL tolsum;
  REAL tol;
  REAL tolexp;
  REAL tolp;              /* Tolerance                                */
  REAL srydpb;            /* square root of (0.5 * ydpb)              */
  REAL dsmall;            /* Machine constant                         */
  REAL dlarge;            /* largest floating point number            */
  REAL tmp;               /* aux variable                             */
  int  j;                 /* loop counter                             */
  int  k;                 /* loop counter                             */
  int  lk;                /* number of iterations in computing        */
                          /* Lipschitz constant                       */
  void *vmblock;          /* Listof dynamic allocations               */


                      /* ----------- dynamic allocations  ----------- */

  vmblock = vminit();
  hilf = (REAL **)vmalloc(vmblock, MATRIX, 5, n);
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return -5;
  }


  /* --------- upper bound for second derivative (dfdxb)    --------- */
  /* --------- using difference quotient and upper bound of --------- */
  /* --------- first derivative (fbnd)                      --------- */

  dsmall = MACH_EPS;
  dlarge = POSMAX;
  dx     = beta - x;
  absdx  = FABS(dx);
  relper = POW(dsmall, (REAL)0.375);
  tmp    = min(relper * FABS(x), absdx);
  tmp    = max(tmp, (REAL)100.0 * dsmall * FABS(x));
  da     = SIGN(tmp, dx);
  if (FABS(da) < dsmall)
    da = relper * dx;
  dgl(x + da, y, hilf[0]); ++*aufrufe;
  dgl(x,      y, hilf[4]); ++*aufrufe;
  for (j = 0; j < n; j++)
    hilf[1][j] = hilf[0][j] - hilf[4][j];
  delf  = norm_max(hilf[1], n);
  dfdxb = dlarge;
  if (delf < dlarge * FABS(da))
    dfdxb = delf / FABS(da);
  fbnd = norm_max(hilf[0], n);

  /* -------- Find estimate of Lipschitz constant  (dbdub)   -------- */
  /* -------- for DE system and choose upper bound for first -------- */
  /* -------- derivative (fbnd)                              -------- */

  dely = relper * norm_max(y, n);
  if (dely < dsmall)
    dely = relper;
  dely = SIGN(dely, dx);
  delf = norm_max(hilf[4], n);
  fbnd = max(fbnd, delf);
  if (delf < dsmall)
  {
    for (j = 0; j < n; j++)
      hilf[2][j] = ZERO,
      hilf[1][j] = ONE;
    delf = ONE;
  }
  else
    for (j = 0; j < n; j++)
      hilf[2][j] = hilf[4][j],
      hilf[1][j] = hilf[4][j];
  dfdub = ZERO;
  lk    = min(n + 1, 3);
  for (k = 1; k <= lk; k++)
  {
    for (j = 0; j < n; j++)
      hilf[3][j] = y[j] + dely * (hilf[1][j] / delf);
    if (k == 2)
    {
      dgl(x + da, hilf[3], hilf[1]);
      ++*aufrufe;
      for (j = 0; j < n; j++)
        hilf[3][j] = hilf[1][j] - hilf[0][j];
    }
    else
    {
      dgl(x, hilf[3], hilf[1]);
      ++*aufrufe;
      for (j = 0; j < n; j++)
        hilf[3][j] = hilf[1][j] - hilf[4][j];
    }
    fbnd = max(fbnd, norm_max(hilf[1], n));
    delf = norm_max(hilf[3], n);
    if (delf >= dlarge * FABS(dely))
    {
      dfdub = dlarge;
      break;
    }
    dfdub = max(dfdub, delf/ FABS(dely));
    if (k < lk)
    {
      if (delf < dsmall)
        delf = ONE;
      for (j = 0; j < n; j++)
      {
        if (k == 2)
        {
          dy = y[j];
          if (FABS(dy) < dsmall)
            dy = dely / relper;
        }
        else
        {
          dy = FABS(hilf[3][j]);
          if (dy < dsmall)
            dy = delf;
        }
        if (FABS(hilf[2][j]) < dsmall)
          hilf[2][j] = hilf[1][j];
        dy         = SIGN(dy, hilf[2][j]);
        hilf[1][j] = dy;
      }
      delf = norm_max(hilf[1], n);
    }
  }
  ydpb = dfdxb + dfdub * fbnd;

  /* ----------- define tolerance (tolp) for computing  ------------- */
  /* ----------- starting step size                     ------------- */

  tolmin = dlarge;
  tolsum = ZERO;
  for (k = 0; k < n; k++)
  {
    tol = relerr * FABS(y[k]) + abserr;
    if (tol < dsmall)
      tol = FABS(dely) * relerr;
    tolexp =  LOG(tol) / LOG(TEN);
    tolmin =  min(tolmin, tolexp);
    tolsum += tolexp;
  }
  tolp = POW(TEN, (HALF * (tolsum / (REAL)n + tolmin) / (qg + ONE)));

  /* -- Anfangsschrittweite und Richtung der Integration bestimmen -- */

  *h = absdx;
  if (ydpb > dsmall || fbnd > dsmall)
  {
    if (ydpb > dsmall)
    {
      srydpb = SQRT(HALF * ydpb);
      if (tolp < srydpb * absdx)
        *h = tolp / srydpb;
    }
    else if (tolp < fbnd * absdx)
      *h = tolp / fbnd;
  }
  else if (tolp < ONE)
    *h = absdx * tolp;
  if (*h * dfdub > ONE)
    *h = ONE /dfdub;
  *h = max(*h, (REAL)100.0 * dsmall * FABS(x));
  if (*h < dsmall)
    *h = dsmall * FABS(beta);
  *h = SIGN(*h, dx);


  vmfree(vmblock);
  return 0;
}



/* ------------------------------------------------------------------ */
/*.BA*/

int einb_rk      /* solve DE system with one of 15 embedding formulas */
/*.IX{einb\unt rk}*/
           (
            REAL      *x,        /* initial/final x-value ............*/
            REAL      beta,      /* desired final x-value ............*/
            int       n,         /* number od DEs ....................*/
            dglsysfnk dgl,       /* right hand side of system ........*/
            REAL      y[],       /* initial y-value at x .............*/
            REAL      abserr,    /* absolute error bound .............*/
            REAL      relerr,    /* relative error bound .............*/
            int       neinb,     /* Number of embedding formula ......*/
            int       hullstp,   /* step size control according to    */
                                 /* Hull? ............................*/
            int       neu,       /* do not use old data ? ............*/
            int       save,      /* save data for a future call ? ....*/
            long      fmax,      /* maximal # of calls of  dgl() .....*/
            long      *aufrufe   /* actual # of calls of dgl() .......*/
           )                     /* error code .......................*/

/***********************************************************************
* Solve a first order system of ordinary differential equations using  *
* the Runge-Kutta embedding formulas over the interval [x0,beta].      *
* The parameter neinb allows to choose among the embedding formulas    *
* while the parameter hullstp allows to select among the step size     *
* controls.                                                            *
.BE*)
* On first call, the parameter  neu must differ from zero. The function*
* then checks the input and determines several aux parameters for the  *
* integration.                                                         *
* If this function has a return value of -2, the maximum number of     *
* function calls was exceeded before reaching beta. One may then call  *
* up the function again (with neu = 0) for the old parameters.         *
* If the user wants to know several intermediate results in the inter- *
* val, the first call must be started with neu != 0 and subsquent calls*
* can be generated in a loop with neu = 0 and different beta.          *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        Starting x-value                                            *
* beta     final desired x-value                                       *
* n        number of differential equations                            *
* dgl      right hand side of system of DEs                            *
* y        [0..n-1] vector of y-values at x                            *
* abserr   absolute error bound  (abserr >= 0). If abserr = 0  we test *
*          only for the relative error                                 *
* relerr   relative error bound  (relerr >= 0). If relerr = 0  we test *
*          only for the absolute error                                 *
* neinb    the number for the desired embedding formula:               *
*            =  0: RK3(2)                                              *
*            =  1: RKF4(3)   (k1=km)                                   *
*            =  2: RKF5(4)                                             *
*            =  3: RK5(4)6M                                            *
*            =  4: RKE5(4)                                             *
*            =  5: HIHA5     (k1=km)                                   *
*            =  6: RK5(4)7S  (k1=km)                                   *
*            =  7: RK5(4)7M  (k1=km)                                   *
*            =  8: RK5(4)7C  (k1=km)                                   *
*            =  9: RK6(5)8M                                            *
*            = 10: RK6(5)8S                                            *
*            = 11: RK6(5)8C                                            *
*            = 12: RKV6(5)                                             *
*            = 13: RKF6(5)A                                            *
*            = 14: RKF6(5)B                                            *
*            = 15: RKC6(5)    (k1=km)                                  *
*            = 16: RKV6(5)9A  (k1=km)                                  *
*            = 17: RKV6(5)9B  (k1=km)                                  *
*            = 18: RKV7(6)                                             *
*            = 19: RK8(7)13M                                           *
*            = 20: RKF8(7)                                             *
*            = 21: RKV8(7)                                             *
*            = 22: RKV9(8)                                             *
*          In the above list, `k1=km' means that the coefficient       *
*          list of A~ and the last row of the matrix b coincide, so    *
*          that the value for k1 need not be computed after a success- *
*          ful integration step. Instead the value of km from the      *
*          previous step can be used. This is indicated by the flag    *
*          `k1_km' inside the coefficient declarations  `koefftyp'.    *
* hullstp  choses the method of step size control:                     *
*             = 0: standard initial value solver method                *
*            != 0: formula of Hull                                     *
* neu      switch for first or repeated calls:                         *
*             = 0: repeated call                                       *
*            != 0: first call                                          *
* save      = 0: All values saved in the local static variable koeff   *
*                and in h are marked invalid before return.            *
*          != 0: The contents of koeff and von h are saved for a       *
*                subsequent call in order to continue integration.     *
*          If a call ends in an error other than -2, the save flag is  *
*          ignored.                                                    *
* fmax     upper bound for the number of allowed function evaluations  *
*          of the right hand side via dgl()                            *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value of integration (normally x = beta)            *
* y        [0..n-1] vector with y-values at x                          *
* aufrufe  number of calls of dgl()                                    *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
*   =   0:  all ok                                                     *
*                                                                      *
* error in input parameters:                                           *
*   Bit 0 on:  wrong value for neinb or fmax                           *
*   Bit 1 on:  abserr or relerr too small                              *
*   Bit 2 on:  interval of integration too small                       *
*   Bit 3 on:  x or beta cannot be represented in computer             *
*   Bit 4 on:  n <= 0                                                  *
*   Bit 5 on:  initial conditions cannot be represented                *
* Example: error code = 21  =>  improper value for neinb or fmax,      *
*                               interval of integration too small and  *
*                               n <= 0                                 *
*                                                                      *
* errors during application of embedding formula (negative):           *
*   =  -1:  The desired relative accuracy is below 100 * machine       *
*           constant in certain regions of the interval of             *
*           integration. Computations were carried out with            *
*           100 * machine constant as the absolute error bound there.  *
*   =  -2:  maximal number of steps reached.                           *
*   =  -3:  possible overflow                                          *
*   =  -4:  computed step size too small                               *
*   =  -5:  lack of available storage space                            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, dglsysfnk, koefftyp, ZERO, MACH_EPS, POSMAX, FABS, norm_max,   *
* init_koeff, hstart, awpl, NE                                         *
.BA*)
***********************************************************************/
/*.BE*/

{
  static                       /* Structure describing the chosen     */
    koefftyp koeff;            /* embedding formula                   */
  static
    REAL     h = ZERO;         /* current step size; h = 0 indicates  */
                               /* that no previous setting for h is   */
                               /* available as this is the first call */
                               /* or save was not activated           */

  REAL       eps;              /* 100 * machine constant              */
  REAL       xzi;              /* largest representable number / 100  */
  int        fehler;           /* error code: wrong input parameters  */
                               /* or error in hstart() or awpl()      */


  eps      =  (REAL)100.0 * MACH_EPS;
  xzi      =  POSMAX;
  xzi      /= (REAL)100.0;
  *aufrufe =  0l;


  /* check input parameters; in order to track several errors at the  */
  /* same time, each type of error is marked by a different bit.      */

  fehler = 0;

  if (neinb < 0 || neinb >= NE ||    /* invalid embedding number      */
      fmax <= 0)                    /* or number of calls for dgl()?  */
    fehler |= 1;
  if (relerr < eps && abserr < eps)      /* error bounds too small?   */
    fehler |= 2;
  if (FABS(beta - *x) < eps)             /* interval too small?       */
    fehler |= 4;
  if (FABS(*x) > xzi ||         /* values at interval end too large?  */
      FABS(beta) > xzi)
    fehler |= 8;
  if (n <= 0)                /* not even one differential equation?   */
    fehler |= 16;
  if (norm_max(y, n) > xzi)      /* initial value too large?          */
    fehler |= 32;


  if (! fehler)                            /* wrong input parameters? */
  {
    if (neu        ||                /* do not use old data           */
        h == ZERO)                   /* or no old data available?     */
    {
      init_koeff(neinb, fmax, &koeff);         /* determine parameters*/
      fehler = hstart(dgl, n, *x, beta, y,     /* for method and      */
                      relerr, abserr,          /* initial step size   */
                      koeff.qg, &h, aufrufe
                     );
    }

    if (! fehler)                          /* no error initializing? */
      fehler = awpl(x, &h, beta, abserr, relerr,   /* perform        */
                    n, dgl, y, hullstp, eps, xzi,  /*  Integration   */
                    &koeff, aufrufe);
  }


  if (! save ||                   /* do not save old data or error in */
      (fehler != 0 &&             /* integration, but not maximal     */
       fehler != -2))             /* step count reached?              */
    h = ZERO;                     /* record for next call             */

  return fehler;                  /* report error code from awpl()    */
}

/* ------------------------- END einb_rk.c -------------------------- */
