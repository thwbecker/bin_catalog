#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#define DOUBLE_PREC

#ifdef DOUBLE_PREC

#define COMP_PRECISION double
#define FLT_FORMAT "%lf"
#define GROUT rg_
#define SROUT rs_
#define EPS_COMP_PREC 5e-15
#else
#define COMP_PRECISION float
#define FLT_FORMAT "%f"
#define GROUT srg_
#define SROUT srs_
#define EPS_COMP_PREC 5e-7
#endif

#define RR 0 // first six are used to initialize symmetric matrices
#define RT 1
#define RP 2
#define TT 4
#define TP 5
#define PP 8
#define TR 3 // these are the (possibly) symmetric parts
#define PR 6
#define PT 7
/* those are the same, really */
#define XX 0 // first six are used to initialize symmetric matrices
#define XY 1
#define XZ 2
#define YY 4
#define YZ 5
#define ZZ 8
#define YX 3 // these are the (possibly) symmetric parts
#define ZX 6
#define ZY 7


// normalize eigenvectors to unity length
#define NORMALIZE

// general real matrix eigenvectors and values
extern void GROUT(int *, int *, COMP_PRECISION *,COMP_PRECISION *, 
		  COMP_PRECISION *, int *, COMP_PRECISION *, int *, 
		  COMP_PRECISION *, int *);
// symmetric real
extern void SROUT(int *, int *, COMP_PRECISION *,COMP_PRECISION *, int *, 
		  COMP_PRECISION *, COMP_PRECISION *,COMP_PRECISION *, int *);

void indexx(int ,COMP_PRECISION *,int *);

void calc_eigensystem_sym_9(COMP_PRECISION *,COMP_PRECISION *,
			    COMP_PRECISION *,
			    short int );
void calc_eigensystem_sym_3x3(COMP_PRECISION [3][3],COMP_PRECISION *,
			      COMP_PRECISION *,
			      short int );
