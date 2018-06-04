#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
/* ------------------------- MODULE cuthill.c ----------------------- */

/***********************************************************************
*                                                                      *
* Solving a linear system with a sparse symmetric system matrix using  *
* -------------------------------------------------------------------  *
* the Cuthill-McKee method                                             *
* ------------------------                                             *
*                                                                      *
* exported function:                                                   *
*   - cutgaucho():  Cuthill-McKee method with Gauss or Cholesky        *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Elmar Pohl (FORTRAN)                           *
* Implementation:       Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               FORTRAN code                                   *
* Date:                 3.31.1992 - 2.28.1997                          *
*                                                                      *
***********************************************************************/

#include <basis.h>    /*  for  REAL, abs, TRUE, FALSE, boolean, ZERO, */
                      /*       min                                    */
#include <vmblock.h>  /*  for  vmalloc, vmcomplete, vmfree, vminit,   */
                      /*       VEKTOR, VVEKTOR, MATRIX                */
#include <choband.h>  /*  for  chobnd                                 */
#include <u_proto.h>  /*  for  band                                   */
#include <cuthill.h>  /*  for  cutgaucho                              */



/* ------------------------------------------------------------------ */

#ifdef DEBUG
/* print a REAL variable `var' (and its name)                         */
#define zeigrvar(var)  printf("%-12s%3"LZP"g\n", #var": ", var);

/* print an int variable `var' (and its name)                         */
#define zeigivar(var)  printf("%-12s%4d\n", #var": ", var);

/* print a REAL vector (and its name) of length `n'                   */
#define zeigrvek(vek, n)              \
  {                                   \
    int jjj;                          \
    printf("%-12s", #vek": ");        \
    for (jjj = 0; jjj < n; jjj++)     \
      printf("%3"LZP"g", vek[jjj]);   \
    printf("\n");                     \
  }

/* print an int vector (and its name) of length `n'                   */
#define zeigivek(vek, n)              \
  {                                   \
    int jjj;                          \
    printf("%-12s", #vek": ");        \
    for (jjj = 0; jjj < n; jjj++)     \
      printf("%4d", vek[jjj]);        \
    printf("\n");                     \
  }

#define ZEIGPERMAT(gauss, n, m, ap)   \
        zeigpermat(gauss, n, m, ap)

/* print a packed [0..n-1,0..m]-bandmatrix including heading          */
#define zeigpacmat(gauss, n, m, ap)                         \
  if (! gauss)                                              \
  {                                                         \
    int iii;                                                \
    int jjj;                                                \
    printf("\npacked matrix after Cholesky decomposition:\n"); \
    for (iii = 0; iii < n; iii++)                           \
    {                                                       \
      for (jjj = 0; jjj <= m; jjj++)                        \
        if (iii + jjj < n)                                  \
          printf("%9.5"LZP"f", ap[iii][jjj]);               \
      printf("\n");                                         \
    }                                                       \
  }

#else
#define zeigrvar(var)
#define zeigivar(var)
#define zeigrvek(vek, n)
#define zeigivek(vek, n)
#define ZEIGPERMAT(gauss, n, m, ap)
#define zeigpacmat(gauss, n, m, ap)
#endif



/* ------------------------------------------------------------------ */

static
void srtdeg         /* sorts a subset of nodes by increasing degree...*/
/*.IX{srtdeg}*/
    (
     int     node[],      /* vector with node numbers.................*/
     int     ideg[],      /* vector with node degrees.................*/
     int     ibeg,        /* starting index...........................*/
     int     iend         /* final index..............................*/
    )

/***********************************************************************
* sorts a subset of nodes by increasing degree                         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* node  [0..] vector with node numbers. Every element in node must be  *
*       a valid index in ideg.                                         *
* ideg  [0..] vector with node degrees. node node[i] has degree        *
*       ideg[node[i].                                                  *
* ibeg  starting index of the nodes to be sorted                       *
* iend  final index of the nodes to be sorted                          *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* node  sorted node vector                                             *
***********************************************************************/

{
  int node0;
  int ideg0;
  int i, j;                                /* Loop variables for node */
  int j0;


  for (i = ibeg + 1; i <= iend; i++)
  {
    node0 = node[i];
    ideg0 = ideg[node0];

    for (j = j0 = i; j-- != ibeg; j0 = j)
    {
      if (ideg0 >= ideg[node[j]])
        break;
      node[j + 1] = node[j];
    }
    node[j0] = node0;
  }
}



/* ------------------------------------------------------------------ */

static
void cuth1k         /* Determines the component of a graph............*/
/*.IX{cuth1k}*/
    (
     int     iroot,       /* Number of starting node..................*/
     int     istart,      /* starting number..........................*/
     int     neighb[],    /* vector with lists of adjacent nodes......*/
     int     inb[],       /* vector with indices......................*/
     int     ideg[],      /* vector with the degree of each node......*/
     boolean mark[],      /* vector with node markings................*/
     int     icm[]        /* vector with previous components..........*/
    )

/***********************************************************************
* Determines the component of a graph induced by node iroot and its    *
* Cuthill-McKee numbering.                                             *
*                                                                      *
* To determine the Cuthill-McKee numbering of a graph, do not call     *
* this routine, but rather cuthill(), which in turn calls this one     *
* repeatedly until all components of the graph are covered. Note that  *
* cuthill() works for unconnected graphs as well.                      *
*                                                                      *
* By successive calls cuth1k for increasing values of istart, all      *
* entries of icm and mark are asigned.                                 *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* iroot   Number of starting node of the new component of the graph    *
* istart  starting number for the Cuthill-McKee numbering of this      *
*         component                                                    *
* neighb  [0..nv-1] vector with lists of adjacent nodes.               *
*         For i=0,...,n-1 neighb contains the numbers of the           *
*         adjacent nodes to node i in neighb[k], where                 *
*         k=inb[i], inb[i]+1, ..., inb[i+1]-1.                         *
* inb     [0..n] vector with indices for the individual lists in       *
*         neighb. inb[n] must equal the number of elements in neighb.  *
* ideg    [0..n-1] vector, which contains the degree of each node,     *
*         i.e., the number of adjacent nodes                           *
* mark    [0..n-1] vector with node markings.                          *
*         If mark[i] is FALSE, the node i can be used to start         *
*         a new component, otherwise it belongs to a previously        *
*         formed component.                                            *
* icm     [0..n-1] vector with previous components of the graph.       *
*         The Cuthill-McKee numbering of these components are          *
*         at the positions  i=0, ..., istart-1.                        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* mark    as on input, but the nodes of newly formed components are    *
*         marked with TRUE.                                            *
* icm     [0..n-1] vector containing the permutations of the node      *
*         numbering after Cuthill-McKee. For i=0, 1, ..., n-1, icm[i]  *
*         denotes the original number of the node and i denotes its    *
*         Cuthill-McKee number. When calling this function, the        *
*         elements i=istart, isart+1, ..., istart+nnew-1 are           *
*         reassigned. Here nnew stands for the number of nodes of the  *
*         new component.                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* srtdeg, TRUE, boolean                                                *
***********************************************************************/

{
  int newbeg;                    /* starting index of the unmarked    */
                                 /* neighbors of a node               */
  int newend;                    /* end index of the new component,   */
                                 /* as it is formed                   */
  int levbeg;                    /* start of the last level in icm    */
  int levend;                    /* end of the last level in icm      */
  int i;                         /* Loop variable for icm             */
  int j;                         /* Loop variable for neighb          */

  /* -------------------------------------------------------------------
  * This algorithm is similar to the one that constructs the level
  * structure, which here is formed inside the vector icm. Additionally
  * we order the node lists which are adjoined to a level by increasing
  * degree.
  * ----------------------------------------------------------------- */

  icm[istart] = iroot;
  newend      = istart;
  levend      = istart - 1;
  mark[iroot] = TRUE;

  do                      /* form the level structure of iroot in icm */
  {
    levbeg = levend + 1;
    levend = newend;
    /* levbeg points to the start, levend to the end of the most  */
    /* recently found level in icm. (The first level consists of  */
    /* iroot).                                                    */

    /* -----------------------------------------------------------------
    * determine the nodes for the next level :
    * find all unmarked neighbors of nodes of the last level and
    * store in  icm
    * --------------------------------------------------------------- */
    for (i = levbeg; i <= levend; i++)
    {
      /* ---------------------------------------------------------------
      * Find the list of unmarked neighbors of the node with
      * original number icm[i]
      *  ------------------------------------------------------------ */

      newbeg = newend + 1;          /* starting index for list in icm */
      /* newbeg points to the beginning of the new level in icm       */
      /* newend always points to the most recently found node in      */
      /* icm.                                                         */

      for (j = inb[icm[i]]; j < inb[icm[i] + 1]; j++)
        if (! mark[neighb[j]])                      /* node unmarked? */
          newend++,
          icm[newend]     = neighb[j],
          mark[neighb[j]] = TRUE;

      /* sort the elements in icm[newbeg],..., icm[newend] */
      /* by increasing degree                              */
      srtdeg(icm, ideg, newbeg, newend);
    }
  }
  while (newend > levend);         /* as long as new nodes are found  */
}



/* ------------------------------------------------------------------ */

static
void lvstru         /* Construct the level structure..................*/
/*.IX{lvstru}*/
    (
     int     iroot,       /* starting node number.....................*/
     int     neighb[],    /* vector with the list if adjacent nodes...*/
     int     inb[],       /* vector with indices......................*/
     boolean mark[],      /* vector with node markings................*/
     int     *nlv,        /* number of levels.........................*/
     int     level[],     /* vector with nodes of the same level......*/
     int     ilv[],       /* vector with indices for the level........*/
     int     *lvnodes     /* number of nodes in the component.........*/
    )

/***********************************************************************
* Construct the level structure for the component created by iroot     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* iroot   starting node number                                         *
* neighb  [0..nv-1] vector with the list if adjacent nodes             *
*         For i=0, 1, ..., n-1 neighb contains the indices             *
*         of the neighbors of node i in  neighb[k] for                 *
*         k=inb[i], inb[i]+1, ..., inb[i+1]-1.                         *
* inb     [0..n] vector with indices for the individual lists in       *
*         neighb. inb[n] equals the number of elements in neighb.      *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* mark     [0..n-1] vector with node markings                          *
* nlv      number of levels                                            *
* level    [0..n-1] vector with nodes of the same level.               *
*          For i=0, 1,..., nlv-1 level contains the indices            *
*          of the nodes of level i in level[k] for                     *
*          k = ilv[i], ilv[i]+1, ..., ilv[i+1]-1.                      *
* ilv      [0..nlv] vector with indices for the level lists in level   *
*          ilv[nlv] equals the number of elements in level.            *
*          If the graph is connected, this number is n.                *
* lvnodes  number of nodes in the component                            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* TRUE, FALSE, boolean                                                 *
***********************************************************************/

{
  int my;                  /* number of nodes found thus far - 1      */
  int levbeg;              /* points to the start of the most recently*/
                           /* found level (nlv-1) in level            */
  int levend;              /* points to the end of the most recently  */
                           /* found level (nlv-1) in level            */
  int i;                   /* indexing variable for level             */
  int j;                   /* indexing variable for neighb            */


  *nlv        = 0;
  level[0]    = iroot;
  mark[iroot] = TRUE;
  my          = 0;
  levend      = -1;
  do
  {
    levbeg    = levend + 1;
    levend    = my;
    ilv[*nlv] = levbeg;
    (*nlv)++;                     /* nlv denotes the number of levels */
                                  /* thus far found. (The first level */
                                  /* consists entirely of iroot.)     */

    /* -----------------------------------------------------------------
    * Find the nodes of the next level :
    * search all unmarked neighbors of nodes of level nlv-1
    * ans store in level
    * --------------------------------------------------------------- */
    for (i = levbeg; i <= levend; i++)
      for (j = inb[level[i]]; j < inb[level[i] + 1]; j++)
        if (! mark[neighb[j]])
          my++,
          level[my]       = neighb[j],
          mark[neighb[j]] = TRUE;
  }                            /* as long as unmarked nodes are found */
  while (my > levend);

  /* The level structure of the component of iroot is now established */
  *lvnodes  = levend + 1;        /* number of nodes in this component */
  ilv[*nlv] = *lvnodes;

  for (i = 0; i < *lvnodes; i++) /* erase all node markings of the run*/
    mark[level[i]] = FALSE;
}



/* ------------------------------------------------------------------ */

static
void fndroo         /* construct optimum structure and start node.....*/
/*.IX{fndroo}*/
    (
     int     *iroot,   /* Number of a node that lies in this component*/
     int     neighb[], /* vector with the list........................*/
     int     inb[],    /* vector with indices.........................*/
     int     ideg[],   /* vector with the degree of each node.........*/
     boolean mark[],   /* vector for marking nodes....................*/
     int     *nlv,     /* number of levels............................*/
     int     level[],  /* vector with lists of nodes of the same level*/
     int     ilv[],    /* vector with indices for the level...........*/
     int     *lvnodes  /* number of nodes in this component...........*/
    )

/***********************************************************************
* Construct the level structure of the component of iroot of a graph;  *
* then try to select the starting node so that a structure with as many*
* as possible levels is generated.                                     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* iroot   Number of a node that lies in this component                 *
* neighb  [0..nv-1] vector with the list of adjacent nodes             *
*         For i=0, 1, ..., n-1, neighb contains the indices            *
*         of neighbors of node i in neighb[k] where                    *
*         k=inb[i], inb[i]+1, ...,  inb[i+1]-1.                        *
* inb     [0..n] vector with indices for the individual lists in       *
*         neighb. inb[n] equal the number of elements in neighb        *
* ideg    [0..n-1] vector with the degree of each node                 *
* mark    [0..n-1] vector for marking nodes                            *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* iroot    new starting node                                           *
* mark     [0..n-1] vector for marking nodes                           *
* nlv      number of levels                                            *
* level    [0..n-1] vector with lists of nodes of the same level.      *
*          For i=0, 1, ..., nlv-1, level contains the indices          *
*          of the nodes of level i in level[k] for                     *
*          k = ilv[i], ilv[i]+1, ..., ilv[i+1]-1.                      *
* ilv      [0..nlv] vector with indices for the level lists in level   *
*          ilv[nlv] always equals the number of elements in            *
*          level. If the graph is connected, then ilv[nlv] =  n        *
* lvnodes  number of nodes in this component                           *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* lvstru                                                               *
***********************************************************************/

{
  int nlvold;   /* maximal level length found thus far                */
  int imin;     /* Index for level, so that level[imin] is the index  */
                /* of a node of minimal degree at the maximal level   */
  int idegmin;  /* degree of the node indexed imin                    */
  int i;        /* indexing variable for level                        */


  for (nlvold = 0; ; )
  {
    lvstru(*iroot, neighb, inb, mark, nlv, level, ilv, lvnodes);

    if (*nlv <= nlvold)         /* level length increased during last */
                                /* step ?                             */
      return;                   /* stop the algorithm                 */

    nlvold = *nlv;

    imin    = ilv[*nlv - 1];          /* search for a node of minimal */
    idegmin = ideg[level[imin]];      /* degree at the last level     */
    for (i = ilv[*nlv - 1] + 1; i < ilv[*nlv]; i++)
      if (ideg[level[i]] < idegmin)
        imin    = i,
        idegmin = ideg[level[i]];

    *iroot = level[imin];       /* use this node to start a new level */
  }                             /* structure                          */
}



/* ------------------------------------------------------------------ */

static
int cuthill         /* Determines the CM-permutation of a graph.......*/
/*.IX{cuthill}*/
    (
     int n,        /* number of nodes of the graph....................*/
     int neighb[], /* vector with lists of adjacent nodes.............*/
     int inb[],    /* vector with indices.............................*/
     int ideg[],   /* vector with the degree of every node............*/
     int icm[],    /* vector with the permutation of the node indices.*/
     int icmrev[]  /* vector with the inverse permutation of icm......*/
    )

/***********************************************************************
* Determines the Cuthill-McKee permutation of a graph.                 *
* This numbering is used to solve sparse symmetric linear systems      *
* in order to economize on memory requirements and computing time.     *
* By using the Cuthill-McKee permutation on the incidence graph of the *
* symmetric system matrix, this matrix is transformed to a symmetric   *
* band matrix, whose bandwidth usually is smaller than that of the     *
* original matrix.                                                     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n       number of nodes of the graph. The graph is given in the      *
*         following two vectors:                                       *
* neighb  [0..nv-1] vector with lists of adjacent nodes.               *
*         For i=0,1,...,n-1, neighb contains the indices of            *
*         the neighbors of node i in neighb[k] where                   *
*         k=inb[i], inb[i]+1, ..., inb[i+1]-1.                         *
* inb     [0..n] vector with indices for the individual lists in       *
*         neighb. inb[n] equals the number of elements in neighb.      *
* ideg    [0..n-1] vector with the degree of every node                *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* icm     [0..n-1] vector with the permutation of the node indices     *
*         according to Cuthill-McKee. For i=0,1,...,n-1, icm[i]        *
*         denotes the original number of a node, while i is its        *
*         Cuthill-McKee number                                         *
* icmrev  [0..n-1] vector with the inverse permutation of icm.         *
*         For i=0,1,...,n-1, icmrev[i] denotes the Cuthill-            *
*         McKee number of the node originally indexed by i.            *
*                                                                      *
* Remark : One of the vectors icm or icmrev is of course redundant.    *
*          Both vectors are kept in order to facilitate condensing the *
*          matrix without laborious look up procedures in only one of  *
*          them.                                                       *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 1: lack of storage space                                           *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* fndroo, cuth1k, FALSE, boolean, vminit, vmalloc, vmcomplete, vmfree, *
* VEKTOR                                                               *
***********************************************************************/

{
  boolean *mark;         /* [0..n-1] vector for node markings         */
  int     *level;        /* [0..n-1] vector for level structure of one*/
                         /* component                                 */
  int     *ilv;          /* [0..n] vector for starting index of levels*/
                         /* in vector level                           */
  int     lvnodes;       /* number of nodes in the level structure    */
                         /* from fndroo()                             */
  int     nfound;        /* total number of nodes of the components,  */
                         /* which have thus far been found by fndroo()*/
  int     iroot;         /* Starting node for one level               */
  int     nlv;           /* number of levels constructed in fndroo()  */
  int     i;             /* indexing variables for mark, icm, icmrev  */
  void    *vmblock;      /* List of dynamically allocated vectors     */


  vmblock = vminit();
  mark  = (boolean *)vmalloc(vmblock, VVEKTOR, n,     sizeof(boolean));
  level = (int *)    vmalloc(vmblock, VVEKTOR, n,     sizeof(int));
  ilv   = (int *)    vmalloc(vmblock, VVEKTOR, n + 1, sizeof(int));
  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return 1;
  }

  for (i = 0; i < n; i++)
    mark[i]   = FALSE,
    icm[i]    = 0,
    icmrev[i] = 0;

  for (i = 0, nfound = 0; i < n; i++)
    if (! mark[i])                              /* node i unlabelled? */
    {                    /* start a new component of the graph  by    */
      /* searching for a starting node iroot in the component         */
      /* of i which has a level structure of maxo=imal length         */
      iroot = i;
      fndroo(&iroot, neighb, inb, ideg, mark, &nlv, level, ilv,
             &lvnodes);

      /* compute the Cuthill-McKee numbering of this component   */
      /* for the starting node and starting index nfound         */
      cuth1k(iroot, nfound, neighb, inb, ideg, mark, icm);

      nfound += lvnodes;
    }

                                 /* now form the inverse permutation */
  for (i = 0; i < n; i++)
    icmrev[icm[i]] = i;
  zeigivek(icm, n);
  zeigivek(icmrev, n);


  vmfree(vmblock);
  return 0;                                         /* successful run */
}



/* ------------------------------------------------------------------ */

static
int symmetrisch     /* Determine whether the matrix is symmetric......*/
/*.IX{symmetrisch}*/
    (
     int     i,      /* ic[i] is the column index of A................*/
     int     j,      /* j is the row index of A.......................*/
     REAL    v[],    /* vector containing the nonzero matrix elements.*/
     int     ic[],   /* vector with the column indices................*/
     int     ir[],   /* vector with indices of beginning rows in ic...*/
     boolean gauss   /* flag, which indicates the method of solution..*/
    )

/***********************************************************************
* Determine whether the matrix described by ic, ir and v is symmetric  *
* at (ic[i],j), i.e. check whether the matrix element A[ic[i]][j]      *
* appears in v and is equal to  A[j][ic[i]].                           *
* If  Gaussian elimination will be used, it is sufficient to check     *
* whether A[ic[i]][j] appears in v. (called weakly symmetric)          *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* i      ic[i] is the column index of A                                *
* j      j is the row index of A                                       *
* v      [0..nv-1] vector containing the nonzero matrix elements.      *
*        (nv denotes the length of v)                                  *
* ic     [0..nv-1] vector with the column indices of the nonzero       *
*        matrix elements.                                              *
* ir     [0..n] vector with indices of beginning rows in ic.           *
*        ir[n] equals the length of ic.                                *
*        (n is the size of the matrix.)                                *
* gauss  Flag, which indicates the method of solution of the linear    *
*        system : Gauss or Cholesky.                                   *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = FALSE: The matrix is not symmetric at (ic[i],j) in case of         *
*          Cholesky or not wekly symmetric for Gauss.                  *
* = TRUE:  Matrix ok at  (ic[i],j).                                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, boolean                                                        *
***********************************************************************/

{
  int k;                                  /* index variable for  ic   */
  int endwert;                            /* stopping value of k loop */


  for (k = ir[ic[i]],                   /* loop over all nonzero      */
       endwert = ir[ic[i] + 1];         /* elemente of A in row ic[i] */
       k < endwert; k++)
    if (ic[k] == j)      /* does A[ic[i]]][j] appear at index k in v? */
      if (gauss)                   /* using Gaussian limination?      */
        return TRUE;               /* occurence in v sufficient       */
      else                         /* Cholesky method used?           */
        return v[i] == v[k];       /* true symmetry needed            */


  return FALSE;/* element does not occur in v => not weakly symmetric */
}



/* ------------------------------------------------------------------ */

static
int graphen_bauen   /* Construct the graph of a symmetric matrix......*/
/*.IX{graphen\unt bauen}*/
    (
     int     n,        /* size of the sparse matrix...................*/
     REAL    v[],      /* vector with nonzero elements of matrix......*/
     int     ic[],     /* vector with column indices..................*/
     int     ir[],     /* vector with the indic. of beg. rows in ic...*/
     boolean gauss,    /* Flag, for system solver to be used..........*/
     int     neighb[], /* vector with indices of adjacent nodes.......*/
     int     inb[],    /* vector with indices for neighb..............*/
     int     ideg[]    /* vector with the degree of each node.........*/
    )

/***********************************************************************
* Construct the graph of a symmetric matrix; prepare for Cuthill-McKee *
* method.                                                              *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n      size of the system matrix                                     *
* v      [0..nv-1] vector with nonzero elements of matrix.             *
* ic     [0..nv-1] vector with column indices of the nonzero matrix    *
*        elements                                                      *
* ir     [0..n] vector with the indices of beginning rows in ic.       *
*        ir[n] equals the length of ic.                                *
* gauss  Flag, for system solver to be used: Gauss or Cholesky.        *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* neighb  [0..nv-1] vector with indices of adjacent nodes. Each row i  *
*         corresponds to a node i. If i != k, node i is                *
*         adjacent to node k if the matrix element a[i][k] != 0        *
*         By symmetry if node i is adjacent to node k, so is           *
*         k to i.                                                      *
*         For i=0,1,...,n-1, neighb contains the indices of            *
*         adjacent nodes from position inb[i] to inb[i+1]-1).          *
* inb     [0..n] vector with indices for neighb.                       *
*         For i=0,1,..,n-1  inb[i] is the starting index of the        *
*         neighbors of node i in neighb. This list continues until     *
*         index inb[i+1]-1. inb[n] equals the length of neighb, which  *
*         is usually less than nv.                                     *
* ideg    [0..n-1] vector with the degree of each node                 *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all is ok                                                       *
* = 1: Matrix not symmetric (Cholesky) or not weakly symmetric (Gauss).*
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, symmetrisch, boolean                                           *
***********************************************************************/

{
  int i;    /* row index of A                                         */
  int j;    /*     index of ic                                        */
  int my;   /* counter of nonzero off-diagonal matrix elements,       */
            /* serves as index for neighb, and supplies values ro inb */


  for (i = my = 0; i < n; i++)          /* for all nodes of the graph */
  {
    inb[i] = my;                           /* starting index of the   */
                                           /* neighbor list of node i */
    for (j = ir[i]; j < ir[i + 1]; j++)    /* for all nonzero entries */
                                           /* in row i                */
      if (ic[j] != i)                      /* not on diagonal?        */
      {
        neighb[my] = ic[j],           /* label the node ic[j] as      */
        my++;                         /* adjacent to node i in neighb */
        if (! symmetrisch(j, i, v, ic, ir, gauss))
          return 1;
      }
    ideg[i] = my - inb[i];           /* number of neighbors of node i */
  }

  inb[n] = my;             /* number of nonzero off-diagonal elements */
  zeigivek(neighb, my);
  zeigivek(inb, n + 1);
  zeigivek(ideg, n);


  return 0;                                      /* run finished      */
}



/* ------------------------------------------------------------------ */

static
int mach_zeilenindizes /* prepare graph construction .................*/
/*.IX{mach\unt zeilenindizes}*/
    (
     int     n,           /* size of matrix...........................*/
     int     *nv,         /* original length of ic and v..............*/
     int     ic[],        /* vector of column indices.................*/
     REAL    v[],         /* vector with nonzero matrix elements......*/
     int     ir[]         /* vector with indices for v................*/
    )

/***********************************************************************
* store the starting indices in v of the matrix rows in  ir.           *
* v now only contains the nonzero elements of A. ic is indexed from 0  *
* on. This is a preparation for forming the graph.                     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n   size of matrix                                                   *
* nv  original length of ic and v                                      *
* ic  [0..nv-1] vector of column indices, ordered row wise, of the     *
*     elements in v. In addition, to mark the rows, at the end of      *
*     each row, ic has one entry of 0.                                 *
* v   [0..nv-1] vector with nonzero matrix elements whereever ic       *
*     has a nonzero index.                                             *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* nv  new length of ic and v                                           *
* ic  [0..nv-1] vector of column indices, ordered row wise, of the     *
*     nonzero matrix elements in v                                     *
* v   [0..nv-1] vector of nonzero matrix elements                      *
* ir  [0..n] vector with indices for v, that start a new matrix row    *
*     ir[n] contains the value of nv.                                  *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 1: ic does not have as many zero indices as the value n.           *
* = 2: There is no matrix row, for which the column index in ic        *
*      is not monotonically increasing or in which a column index      *
*      larger than n-1 was used.                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  int i;           /* Index for ir (row counter for matrix)           */
  int j;           /* index for ic and iv; after the first sweep this */
                   /* equals the number of nonzero matrix elements    */
                   /* or the new length of ic and v.                  */
  int k;           /* read index for ic and iv                        */
  int maxic;       /* maximal column index of a row, to check on      */
                   /* monotonicity and range of the column indices in */
                   /* ic                                              */


  if (ic[*nv - 1] != 0)                  /* no zero at the end of ic? */
    return 1;

  for (i = j = k = 0; k < *nv; k++, i++)           /* over all of  ic */
  {
    ir[i] = j;           /* store the index needed for v of the first */
                         /* nonzero entry row i in ir[i] eintragen    */
    for (maxic = 0; ic[k] != 0; k++, j++)      /* for all nonzero     */
    {                                          /* elements in row i   */
      if (ic[k] <= maxic)                      /* not monotone?       */
        return 2;                              /* return error        */
      else                              /* new index exceeds maximum? */
        maxic = ic[k];                  /* use as new maximum         */
      ic[j] = ic[k] - 1;     /* change book keeping info in ic and v  */
      v[j]  = v[k];
    }
    if (maxic > n)                      /* column indes too large?    */
      return 2;                         /* report error               */
  }

  if (i != n)                  /* computed number of rows i incorrect */
    return 1;

  ir[n] = j;
  *nv   = j;


  return 0;                                         /* run successful */
}



/* ------------------------------------------------------------------ */

static
void permut         /* permute REAL vector xold to xnew via perm......*/
/*.IX{permut}*/
    (
     int     n,           /* size of the matrix.......................*/
     int     perm[],      /* permutation vector.......................*/
     REAL    xold[],      /* vector to be permuted....................*/
     REAL    xnew[]       /* permuted vector..........................*/
    )

/***********************************************************************
* permute vector vector xold to xnew via perm.                         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n     size of the matrix                                             *
* perm  [0..n-1] permutation vector                                    *
* xold  [0..n-1] vector to be permuted                                 *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* xnew  [0..n-1] permuted vector                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL                                                                 *
***********************************************************************/

{
  while (n-- != 0)
    xnew[perm[n]] = xold[n];
}



/* ------------------------------------------------------------------ */

static
int ibdwid          /* Compute the bandwidth of a matrix..............*/
/*.IX{ibdwid}*/
    (
     int     n,           /* size of the matrix.......................*/
     int     neighb[],    /* vector with indices of adjacent nodes....*/
     int     inb[],       /* vector with indices for neighb...........*/
     int     nold[],      /* Cuthill-Mc-Kee-permutation vector........*/
     int     nnew[]       /* inverse permutation of nold..............*/
    )

/***********************************************************************
* Compute the bandwidth of a matrix for the given permutation of nodes *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n      \   Description of the graph                                  *
* neighb  >  (see graphen_bauen())                                     *
* inb    /                                                             *
* nold       [0..n-1] permutation vector. For i=0, ..., n-1,           *
*            nold[i] contains the original index of the                *
*            node now numbered i                                       *
* nnew       [0..n-1] inverse permutation of nold.                     *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* bandwidth, given by the maximal distance of two neighbors            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* abs                                                                  *
***********************************************************************/

{
  int maxbw;
  int diff;
  int i;
  int j;


  for (i = 0, maxbw = 0; i < n - 1; i++)
    for (j = inb[nold[i]]; j < inb[nold[i] + 1]; j++)
      if ((diff = abs(i - nnew[neighb[j]])) > maxbw)
        maxbw = diff;

  return maxbw;
}



/* ------------------------------------------------------------------ */

static
int checkbdwid      /* prevention of increasing of the band width.....*/
/*.IX{checkbdwid}*/
    (
     int     n,        /* size of the sparse matrix...................*/
     int     neighb[], /* vector with indices of adjacent nodes.......*/
     int     inb[],    /* vector with indices for neighb..............*/
     int     m,        /* halfbandwidth of the CM-numerated matrix....*/
     int     icm[],    /* Cuthill-Mc-Kee-permutation vector...........*/
     int     icmrev[]  /* invers Cuthill-Mc-Kee-permutation vector....*/
    )

/***********************************************************************
* test whether the Cuthill-McKee-permutation increases the bandwidth   *
* of the matrix instead of decreasing it. In this unfortunate case the *
* matrix shall not be permutated.                                      *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n      \   description of the graph                                  *
* neighb  >  (see graphen_bauen())                                     *
* inb    /                                                             *
* m          halfbandwidth of the Cuthill-McKee numerated matrix       *
* icm        [0..n-1]-Cuthill-McKee-permutation vector.                *
*            icm[i] (i=0(1)n-1) is the original number of the node     *
*            with the new number i.                                    *
* icmrev     [0..n-1]-invers CM-permutation vector of icm.             *
*            icmrev[i] (i=0(1)n-1) is the new number of the node with  *
*            the old number i.                                         *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* icm        [0..n-1]-vector with identical permutation, if m had to   *
*            be corrected, otherwise the CM-Permutation                *
* icmrev     [0..n-1]-vector containing the inverse permutation of icm *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* corrected band width, i.e. maximum distance between two neighbours:  *
* either m, if m <= original band width, or the original band width.   *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* abs                                                                  *
***********************************************************************/

{
  int morig;              /* halfbandwidth of the CM-numerated matrix */
  int diff;
  int i;
  int j;


  for (i = 0, morig = 0; i < n - 1; i++)
    for (j = inb[i]; j < inb[i + 1]; j++)
      if ((diff = abs(i - neighb[j])) > morig)
        morig = diff;
  zeigivar(morig);

  if (m > morig)
  {
    m = morig;
    for (i = 0; i < n; i++)
      icm[i] = icmrev[i] = i;
  }

  return m;
}



/* ------------------------------------------------------------------ */

static
void cutpak         /* build full perm. matrix for Cholesky...........*/
/*.IX{cutpak}*/
    (
     int     n,           /* size of matrix...........................*/
     int     m,           /* half band width of the matrix............*/
     REAL    v[],         /* nonzero matrix element vector............*/
     int     ic[],        /* column index vector......................*/
     int     ir[],        /* starting indices for each matrix row.....*/
     int     icmrev[],    /* inverse permutation vector...............*/
     REAL    *ap[]        /* full matrix in condensed form............*/
    )

/***********************************************************************
* perform the Cuthill-McKee permutation given in icm on v, ir and ic.  *
* Store the upper band of the resulting matrix in ap.                  *
* ap is then used to solve the linear system via chobnd().             *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n       size of matrix                                               *
* m       number of upper codiagonals after Cuthill-McKee permutation  *
* v       Vector, containing the nonzero matrix elements rowwise       *
* ic      Vector, with the column indices for every element in v       *
* ir      [0..n] vector, of starting indices for each matrix row in    *
*         v and ic. ir[n] equals the number of nonzero matrix          *
*         elements.                                                    *
* icmrev  [0..n-1] vector of the inverse Cuthill-McKee permutation.    *
*         For i=0,...,n-1,  icmrev[i] denotes the Cuthill-             *
*         McKee number of the node with original index  i.             *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* ap      [0..n-1,0..m] matrix in condensed form                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO                                                           *
***********************************************************************/

{
  int irev;
  int jrev;
  int i;
  int j;


  for (i = 0; i < n; i++)
    for (j = 0; j <= m; j++)
      ap[i][j] = ZERO;

  for (i = 0; i < n; i++)
    for (j = ir[i], irev = icmrev[i]; j < ir[i + 1]; j++)
      if ((jrev = icmrev[ic[j]]) >= irev)
        ap[irev][jrev - irev] = v[j];
}



/* ------------------------------------------------------------------ */

static
void cutpk2         /* build full perm. matrix for Gauss..............*/
/*.IX{cutpk2}*/
    (
     int     n,           /* size of matrix...........................*/
     int     m,           /* half band width of the matrix............*/
     REAL    v[],         /* nonzero matrix element vector............*/
     int     ic[],        /* column index vector......................*/
     int     ir[],        /* starting indices for each matrix row.....*/
     int     icmrev[],    /* inverse permutation vector...............*/
     REAL    *ap[]        /* full matrix in condensed form............*/
    )

/***********************************************************************
* perform the Cuthill-McKee permutation given in icm on v, ir and ic.  *
* Store the upper band of the resulting matrix in ap in condensed form.*
* ap is then used to solve the linear system via band().               *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n       size of matrix                                               *
* m       number of upper codiagonals after Cuthill-McKee permutation  *
* v       Vector, containing the nonzero matrix elements rowwise       *
* ic      Vector, with the column indices for every element in v       *
* ir      [0..n] vector, of starting indices for each matrix row in    *
*         v and ic. ir[n] equals the number of nonzero matrix          *
*         elements.                                                    *
* icmrev  [0..n-1] vector of the inverse Cuthill-McKee permutation.    *
*         For i=0,...,n-1,  icmrev[i] denotes the Cuthill-             *
*         McKee number of the node with original index  i.             *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* ap      [0..n-1,0..3*m] matrix in condensed form                     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, ZERO                                                           *
***********************************************************************/

{
  int irev;
  int i;
  int j;


  for (i = 0; i < n; i++)
    for (j = 0; j <= 2 * m; j++)
      ap[i][j] = ZERO;

  for (i = 0; i < n; i++)
    for (j = ir[i], irev = icmrev[i]; j < ir[i + 1]; j++)
      ap[irev][m + icmrev[ic[j]] - irev] = v[j];
}
#ifdef DEBUG



/* ------------------------------------------------------------------ */

static
void zeigpermat     /* print full matrix in condensed form............*/
/*.IX{zeigpermat}*/
    (
     boolean gauss,       /* flag: Gauss oder Cholesky? ..............*/
     int     n,           /* size of the sparse matrix................*/
     int     m,           /* half band width of matrix ...............*/
     REAL    *ap[]        /* full matrix in condensed form. ..........*/
    )

/***********************************************************************
* print full matrix in condensed form                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* gauss   flag: Gauss oder Cholesky?                                   *
* n       size of the sparse matrix                                    *
* m       half band width of matrix                                    *
* ap      `gauss' set:     [0..n-1,0..3*m]-matrix in condensed form    *
*         `gauss' cleared: [0..n-1,0..m]-matrix in condensed form      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* boolean, REAL, printf, min, abs, ZERO                                *
***********************************************************************/

{
  int  abstand;
  int  i;
  int  j;
  int  i2;
  int  j2;


  abstand = (gauss) ? m: 0;

  printf("\npermutated matrix:\n");
  for (i = 0; i < n; i++)            /* full condensed matrix A[i][j] */
  {
    for (j = 0; j < n; j++)
    {
      i2 = min(i, j);
      j2 = abs(j - i);
      printf("%3.2"LZP"g", (j2 < m + 1)            ?
                           (ap[i2][abstand + j2])  :
                           (ZERO)
            );
    }
    printf("\n");
  }
}
#endif



/* ------------------------------------------------------------------ */
/*.BA*/

int cutgaucho       /* sparse matrices via  Cuthill-McKee ............*/
/*.IX{cutgaucho}*/
    (
     boolean gauss,       /* Flag: Gauss or Cholesky .................*/
     int     n,           /* size of the sparse matrix ...............*/
     int     nv,          /* number of nonzero matrix elements + n ...*/
     int     ic[],        /* vector with column indices of v-elements */
     REAL    v[],         /* vector with nonzero matrix elements .....*/
     int     nrs,         /* number of right hand sides ..............*/
     REAL    rs[],        /* vector with all right hand sides ........*/
     REAL    x[],         /* vector with all solutions ...............*/
     int     *m           /* half band width of matrix ...............*/
    )                     /* error code ..............................*/

/***********************************************************************
* Solve a linear sparse and symmetric system of equations using the    *
* Cuthill-McKee method followed by Gaussian elimination with column    *
* pivot search, or if specified, by the Cholesky method.               *
.BE*)
*                                                                      *
* The nonzero matrix elements are stored in v, together with some      *
* bookkeeping zeros. The Cuthill-McKee method is applied to v to reduce*
* the bandwidth of the sparse system matrix. The compressed matrix or  *
* its upper half band is stored in condensed form and the linear       *
* systems are solved via Gaussian elimination with pivot search for    *
* condensed banded matrices using the function band(), or via the      *
* Cholesky method for condensed banded matrices in chobnd(). If the    *
* attempted factorization can be performed, we solve for the given     *
* right hand sides rs using band() or chobnd(). The solutions are      *
* stored in x.                                                         *
*                                                                      *
* If the system matrix is not positive definite, the Cholesky method   *
* will fail. In this case Gaussian elimination must be used. This,     *
* however, requires more computation time and triple the storage of    *
* the condensed matrix.                                                *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* gauss  Flag, TRUE for Gauss, FALSE for Cholesky                      *
* n      size of the sparse [0..n-1,0..n-1] matrix                     *
* nv     number of column indices in ic                                *
*        (= number of nonzero matrix elements + n)                     *
* ic     [0..nv-1] vector with column indices stored row wise followed *
*        by a bookkeeping zero.                                        *
* v      [0..nv-1] vector with the nonzero elements. Each nonzero index*
*        in ic corresponds to a nonzero matrix entry of A in v         *
* nrs    number of right hand sides                                    *
* rs     [0..n*nrs-1] vector with all right hand sides                 *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x  [0..n*nrs-1] vector of all solutions                              *
* m  number of upper codiagonals of the condensed matrix               *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 1: n <= 0  or  nv <= 0  or  nrs <= 0  or    m < 0                  *
* = 2: lack of memory                                                  *
* = 3: the matrix is numerically singular in Gauss, or not strongly    *
*      nonsingular in Cholesky.                                        *
* = 5: ic and n do not match.                                          *
* = 6: There is a matrix row for which the column index in ic is not   *
*      monotonically increasing or in which a column index larger      *
*      n-1 is used.                                                    *
* = 7: matrix not symmetric in Cholesky or not weakly symmetric in     *
*      Gauss.                                                          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* mach_zeilenindizes, graphen_bauen, cuthill, ibdwid, cutpk2,          *
* cutpak, permut, REAL, vminit, vmalloc, vmcomplete, vmfree, VEKTOR,   *
* VVEKTOR, MATRIX, band, chobnd, boolean, checkbdwid                   *
.BA*)
***********************************************************************/
/*.BE*/

{
  void *vmblock;        /* List of dynamically allocated vectors and  */
                        /* matrizes                                   */
  int  *ir;             /* [0..n] vector with starting indices of     */
                        /* rows in v                                  */
  int  *neighb;         /* [0..nv-1] vector with the graph lists      */
  int  *inb;            /* [0..n] vector with pointers to sublists    */
                        /* of neighb                                  */
  int  *ideg;           /* [0..n-1] vector of node degrees            */
  int  *icm;            /* [0..n-1] Cuthill-McKee numbering           */
  int  *icmrev;         /* [0..n-1] inverse permutation to icm        */
  REAL *rshilf;         /* [0..n-1] aux vector used for r.h. sides    */
  int  *perm = NULL;    /* [0..n-1] aux vector for band()             */
                        /* (for row permutation when pivoting)        */
  REAL **ap;            /* [0..n-1,0..3*m] matrix with the upper band */
                        /* of the compressed matrix in condensed form.*/
                        /* A [0..n-1,0..2*m] matrix will not suffice  */
                        /* since the factorization in band() requires */
                        /* additional columns in Gauss.               */
                        /* in Cholesky: [0..n-1,0..m] matrix with the */
                        /* upper half band of the compressed matrix   */
                        /* in condensed form                          */
  int  fehler;          /* error code of cuthill() and chobnd()       */
  int  irs;             /* loop counter for all r.h. sides            */
  int  hilf;            /* number of columns of ap, later replaced by */
                        /* sign of the system matrix determinant, if  */
                        /* band() is used.                            */


  /* ------------------- check input parameters --------------------- */

  if (n <= 0 || nv <= 0 || nrs <= 0)
    return 1;


  /* -------------------- create dynamic vectors -------------------- */

  vmblock = vminit();
  ir     = (int *) vmalloc(vmblock, VVEKTOR, n + 1, sizeof(int));
  neighb = (int *) vmalloc(vmblock, VVEKTOR, nv,    sizeof(int));
  inb    = (int *) vmalloc(vmblock, VVEKTOR, n + 1, sizeof(int));
  ideg   = (int *) vmalloc(vmblock, VVEKTOR, n,     sizeof(int));
  icm    = (int *) vmalloc(vmblock, VVEKTOR, n,     sizeof(int));
  icmrev = (int *) vmalloc(vmblock, VVEKTOR, n,     sizeof(int));
  rshilf = (REAL *)vmalloc(vmblock, VEKTOR,  n,     0);
  if (gauss)
    perm = (int *)   vmalloc(vmblock, VVEKTOR, n,     sizeof(int));
  if (! vmcomplete(vmblock))                       /* lack of memory? */
  {
    vmfree(vmblock);
    return 2;
  }


  /* ----------- prepare and execute Cuthill-McKee method ----------- */

  fehler = mach_zeilenindizes(n, &nv, ic, v, ir);
                                 /* determine ir and reduce ic and v  */
  if (fehler)                    /* inconsistent values in ic and n?  */
  {
    vmfree(vmblock);
    return fehler + 4;
  }

  if (graphen_bauen(n, v, ic, ir, gauss,  neighb, inb, ideg))
                                                        /* form graph */
  {
    vmfree(vmblock);
    return 7;
  }


  fehler = cuthill(n, neighb, inb, ideg, icm, icmrev);
                       /* compute Cuthill-McKee permutation berechnen */
  if (fehler)                                /* lack of memory?       */
  {
    vmfree(vmblock);
    return 2;
  }


  /* ----- apply Cuthill-McKee numbering to the matrix -------------- */

  *m = ibdwid(n, neighb, inb, icm, icmrev);
                  /* Compute half band width of the compressed matrix */

  *m = checkbdwid(n, neighb,         /* einschreiten, falls die neue  */
                  inb, *m,           /* Numerierung die Bandbreite zu */
                  icm, icmrev);      /* vergroessern droht            */

  if (gauss)                           /* assign number of rows in ap */
    hilf = 3 * *m + 1;                 /* according to method used    */
  else
    hilf = *m + 1;
  ap = (REAL **)vmalloc(vmblock, MATRIX, n, hilf);
                                      /* assign matrix ap dynamically */
  if (! vmcomplete(vmblock))          /* lack of memory?              */
  {
    vmfree(vmblock);
    return 2;
  }

  if (gauss)
    cutpk2(n, *m, v, ic, ir, icmrev, ap);
                                         /* condense matrix for Gauss */
  else
    cutpak(n, *m, v, ic, ir, icmrev, ap);
                                      /* condense matrix for Cholesky */
  ZEIGPERMAT(gauss, n, *m, ap);


  /* ----------------- solve all linear systems --------------------- */

  if (gauss)
    fehler = band(1, n, *m, *m, ap, rshilf, perm, &hilf);
                                             /* Compute Gauss factors */
  else
    fehler = chobnd(1, n, *m, ap, rshilf);
                                          /* Compute Cholesky factors */
  if (fehler)                             /* unsuccessful?            */
  {
    vmfree(vmblock);
    return fehler;
  }
  zeigpacmat(gauss, n, *m, ap);

  for (irs = 0; irs < nrs; irs++)               /* for all r.h. sides */
  {
    permut(n, icmrev, rs, rshilf);
                         /* copy icmrev permutation from rs to rshilf */
    if (gauss)
      (void)band(2, n, *m, *m, ap, rshilf, perm, &hilf);
                        /* solve system and store solutionm in rshilf */
    else
      (void)chobnd(2, n, *m, ap, rshilf);
                         /* solve system and store solution in rshilf */
    permut(n, icm, rshilf, x);
          /* reverse the permutation in rshilf and copy solution to x */

    rs += n;                             /* proceed to next r.h. side */
    x  += n;                             /* proceed to next solution  */
  }


  vmfree(vmblock);
  return 0;                                         /* return message */
}

/* ------------------------- END cuthill.c -------------------------- */
