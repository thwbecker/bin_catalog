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

/* matrix */
#define RR 0 // first six are used to initialize symmetric matrices
#define RT 1
#define RP 2
#define TT 4
#define TP 5
#define PP 8

#define TR 3 // these are the (possibly) symmetric parts
#define PR 6
#define PT 7

#ifndef BC_BOOLEAN 
#define BC_BOOLEAN unsigned short
#endif



// normalize eigenvectors to unity length
#define NORMALIZE

// general real matrix eigenvectors and values
extern void GROUT(int *, int *, COMP_PRECISION *,COMP_PRECISION *, 
		  COMP_PRECISION *, int *, COMP_PRECISION *, int *, 
		  COMP_PRECISION *, int *);
// symmetric real
extern void SROUT(int *, int *, COMP_PRECISION *,COMP_PRECISION *, int *, 
		  COMP_PRECISION *, COMP_PRECISION *,COMP_PRECISION *, int *);



void calc_eigensystem_vec6(double *, double *, double *, unsigned short, unsigned short);
void calc_eigensystem_sym_3x3(double [3][3], double *, double *, unsigned short, unsigned short);
void calc_eigensystem_sym_9(double *, double *, double *, unsigned short, unsigned short);
void indexx(int, double *, int *);
void swap(double *, double *);
