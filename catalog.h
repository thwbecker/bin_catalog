#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


/* modes */
#define BC_AKI 0
#define BC_CMT 1 			/* mdat format with time as code */
#define BC_ENG 2			/* lon lat dep mw time */
#define BC_CMT_FP 3 		/* Harvard FP mode */

/* catalog and our stress tensors are is in spherical coordinate system , r theta phi */
#define BC_RR 0			/* vector */
#define BC_RT 1
#define BC_RP 2
#define BC_TT 3
#define BC_TP 4
#define BC_PP 5

#define BC_R 0
#define BC_THETA 1
#define BC_PHI 2

#define BC_NDIM 3 			/* three dimensions */
#define BC_MICHAEL_NPAR  5	/* five parameters for michael inversion */
#ifndef BC_BOOLEAN 
#define BC_BOOLEAN unsigned short
#endif
#define BC_SWITCH unsigned short int
#define BC_TRUE 1
#define BC_FALSE 0

#define BC_CPREC double
#define BC_EPS 5e-15
#define BC_PREC_FMT "%lf"
#define BC_PREC2_FMT "%lf %lf"

#define BC_RGEN ran2

/* from GMT 4.5.18 */
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef HALF_PI
#define HALF_PI 1.570157963267948966192313216916 /* 90 dege */
#endif
#ifndef TWO_PI
#define TWO_PI 6.283185307179586476925286766559
#endif

#ifndef BC_D2R_FAC
#define BC_D2R_FAC 0.0174532925199432957692369076849 // (M_PI / 180.0)
#endif
#ifndef BC_R2D_FAC
#define BC_R2D_FAC 57.295779513082320876798154814105 // (180.0 / M_PI)
#endif

#define BC_R2D(x) ((x)* BC_R2D_FAC)
#define BC_D2R(x) ((x)* BC_D2R_FAC)
#ifndef BC_MEMERROR
#define BC_MEMERROR(x) {fprintf(stderr,"memory allocation error: %s\n",x);exit(-1);}
#endif

/* 


   constants

   
*/
/* from meca */
#define BC_EPSIL 0.0001

/* angle in deg for random fluctuations */
#define BC_EPS_ANGLE_FOR_RANDOM_DEG 30
/*  */
#define BC_FRIC_DEF 0.6
#define BC_FRIC_SCAN_INC 0.025	/* increment for friction scan */
#define BC_MAX_SWEEP_FAC 2	/* multiple of number of bins */
/* 

   catalog closeness
*/

#define BC_DEP_CLOSE 0.4 /* vertical  distance, km */
#define BC_TIME_CLOSE 0.0003 /* 0.00025 is <~1 min */

#define BC_RADIUS 6371.
#define BC_NQUAKE_LIM_FOR_STRESS 2	/* min number of entries, has
					   to be greater */
/*  */
#define BC_MICHAEL_NMC 5000	/* monte carlo for michael inversion, 2000 good number for 95% confidence? */

/* two bail crit */
#define BC_MICHAEL_RSWEEP_MAX 50000	/* max sweep number for random sampling */
#define BC_MICHAEL_RACC 1e-4	/* random sweep accuracy */


/* structures */

/* from GMT 4.5.18 */
struct AXIS {
    double str;
    double dip;
    double val;
    long int e;
};
struct M_TENSOR {
    long int expo;
    double f[6];
};
struct nodal_plane {
    double str;
    double dip;
    double rake;
}; 


/* 

 */

//#define USE_POT

/* quake structure */
struct qke{
  BC_CPREC strike,dip,rake;
  BC_CPREC strike2,dip2,rake2;	/* alternate fault plane */
  BC_CPREC dlon,dlat,ddlon,ddlat,depth,coslat,lkm;
  BC_CPREC lon,lat;		/* radian versions, for now, keep both */
  BC_CPREC mag,m0,m[6],tsec;
  BC_BOOLEAN deleted;
  int exp;
  BC_CPREC plon,plat;
};
/* Kostrov bin structure */
struct bn{
  BC_CPREC lon,lat,dlon,dlat,area,coslat;
  BC_CPREC sumw;			/* sum of weights, number of entries
					   if no weighting*/
  unsigned int n;		  /* number of entries */
  unsigned int *quake ;	  /* quake list */
  BC_CPREC s[6],ds[6];		  /* stress tensor, from randomized scan*/
  BC_CPREC mn[6],m[6];		/* normalized and scaled moment */
  BC_CPREC me,mens,men;		/* mean horizontal strain */
  BC_CPREC smn[6],mnloc[6];		/* std for monte carlo, local
					   realization */
  BC_CPREC inst[3];		/* instability for three different inversions */
  BC_CPREC best_fric;
  BC_CPREC def_s[6],best_s[6];		/* default friction and best friction detemined stress */
  BC_CPREC *weight;
};
/* kostrov summation structure */
struct kostrov_sum{
  /* for binning */
  BC_CPREC dlonmin,dlonmax,dlatmin,dlatmax,ddlon,ddlat;
  BC_CPREC minmag,maxmag,mindepth,maxdepth,mtot;
  BC_BOOLEAN init;
  int weighting_method;
  struct bn *bin;
  int nx,ny,nxny;
};
/* catalog structure */
  
struct cat{
  struct kostrov_sum sum[1];
  struct qke *quake;
  int n;
  BC_CPREC minlond,maxlond,minlatd,maxlatd,lkm_min,lkm_max;
  BC_CPREC minmag,maxmag,mindepth,maxdepth;
  /* those are auto-determined */
  BC_CPREC tcenter,tmin,tmax,trange;
  /* 0: lon-lat coordinates 1: x(e.g. distance along profile)-y(e.g. depth) coordinates */
  BC_BOOLEAN is_xy,use_friction_solve;
  /* for random numbers */
  long int seed;
};


extern void GROUT(int *, int *, double *,double *, 
		  double *, int *, double *, int *, 
		  double *, int *);
// symmetric real
extern void SROUT(int *, int *, double *,double *, int *, 
		  double *, double *,double *, int *);
/* 
 cproto -f2 *.c | grep -v main\( > proto.h
*/
#include "proto.h"
