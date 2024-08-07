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

/* catalog is in spherical coordinate system */
#define BC_RR 0
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

#define BC_BOOLEAN unsigned short
#define BC_TRUE 1
#define BC_FALSE 0

#define BC_CPREC double
#define BC_PREC_FMT "%lf"

#define BC_RGEN ran2

/* from GMT 4.5.18 */
#ifndef TWO_PI
#define TWO_PI        6.28318530717958647692
#endif
#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif
#ifndef BC_D2R
#define BC_D2R (M_PI / 180.0)
#endif
#ifndef BC_R2D
#define BC_R2D (180.0 / M_PI)
#endif
#define atan2d(y,x) (atan2(y,x) * BC_R2D)
#define d_atan2d(y,x) ((x) == 0.0 && (y) == 0.0 ? 0.0 : atan2d(y,x))
/* from meca */
#define BC_EPSIL 0.0001


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
#define BC_NQUAKE_LIM_FOR_STRESS 2	/* min number of entries, has
					   to be greater */
#define BC_MICHAEL_NMC 5000	/* monte carlo for michael inversion, 2000 good number for 95% confidence? */
//#define USE_POT

/* quake structure */
struct qke{
  BC_CPREC strike,dip,rake;
  BC_CPREC strike2,dip2,rake2;	/* alternate fault plane */
  BC_CPREC lon,lat,depth,coslat,lkm;
  BC_CPREC mag,m0,m[6],tsec;
  BC_BOOLEAN deleted;
  int exp;
  BC_CPREC plon,plat;
};
/* Kostrov bin structure */
struct bn{
  float lon,lat,area,coslat;
  BC_CPREC sumw;			/* sum of weights, number of entries
				   if no weighting*/
  unsigned int n;		  /* number of entries */
  unsigned int *quake ;	  /* quake list */
  BC_CPREC s[6],ds[6];		  /* stress tensor */
  BC_CPREC mn[6],m[6];		/* normalized and scaled moment */
  BC_CPREC me,mens,men;		/* mean horizontal strain */
  BC_CPREC smn[6],mnloc[6];		/* std for monte carlo, local
					   realization */
  BC_CPREC *weight;
};
/* kostrov summation structure */
struct kostrov_sum{
  /* for binning */
  BC_CPREC lonmin,lonmax,latmin,latmax,dlon,dlat;
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
  BC_CPREC minlon,maxlon,minlat,maxlat,lkm_min,lkm_max;
  BC_CPREC minmag,maxmag,mindepth,maxdepth;
  /* those are auto-determined */
  BC_CPREC tcenter,tmin,tmax,trange;
  /* 0: lon-lat coordinates 1: x(e.g. distance along profile)-y(e.g. depth) coordinates */
  BC_BOOLEAN is_xy;
  /* for random numbers */
  long int seed;
};




#ifndef BC_MEMERROR
#define BC_MEMERROR(x) {fprintf(stderr,"memory allocation error: %s\n",x);exit(-1);}
#endif

#define BC_PIF 57.295779513082320876798154814105
#define BC_RADIUS 6371.

#define BC_DEG2RAD(x) ((x)/BC_PIF)


/* handle_catalog.c */
void relocate_catalog(struct cat *, struct cat *, struct cat *);
void sum_kostrov_bins(struct cat *, BC_BOOLEAN , BC_BOOLEAN, BC_BOOLEAN);
void setup_kostrov(struct cat *, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC,
		   int);
void clear_bins(struct cat *);
void print_kostrov_bins(struct cat *, char *, BC_BOOLEAN);
void aki2mom(BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC *, BC_CPREC *);
void merge_catalog(struct cat *, struct cat *, struct cat *, BC_CPREC, BC_CPREC);
int print_catalog(char *, struct cat *, int);
int read_catalog(char *, struct cat *, int);
void copy_quake(struct qke *, struct qke *);
int read_quake(FILE *, struct qke *, int);
char *mode_name(int);
void print_quake(FILE *, struct qke, int);
BC_CPREC max_x_from_int_vector(BC_CPREC *, int *, int );
int read_quake_aki(FILE *, struct qke *);
int read_quake_cmt(FILE *, struct qke *);
int read_quake_eng(FILE *, struct qke *);
void rotate_vec6(BC_CPREC *, BC_CPREC *, BC_CPREC, BC_CPREC, BC_CPREC);
void sixsymtomat(BC_CPREC *, BC_CPREC [3][3]);
void mattosixsym(BC_CPREC [3][3], BC_CPREC *);
void get_gen_rot(BC_CPREC [3][3], BC_CPREC, BC_CPREC, BC_CPREC);
void rotate_mat(BC_CPREC [3][3], BC_CPREC [3][3], BC_CPREC [3][3]);
void print_quake_aki(FILE *, struct qke);
void print_quake_cmt(FILE *, struct qke);
void create_catalog(struct cat *,long int);
void make_room_for_quake(struct cat *);
FILE *myopen(char *, char *, char *);
BC_CPREC distance_geo(BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC);
BC_CPREC distance_cart(BC_CPREC, BC_CPREC, BC_CPREC, BC_CPREC);
BC_CPREC distance(struct cat *,struct cat *,int, int);
void remove_trace(BC_CPREC *);
BC_CPREC mean_hor_strain(BC_CPREC *);
void get_index_vector(int **, int, int, long *);
BC_CPREC gauss_ran(long int *, BC_CPREC);
BC_CPREC gasdev(long int *);
BC_CPREC ran2(long int *);
BC_CPREC mag2mom(BC_CPREC);
BC_CPREC mom2mag(BC_CPREC);
BC_CPREC mag2pot(BC_CPREC);
BC_CPREC scalar_mom(BC_CPREC *);
BC_CPREC tensor6_norm(BC_CPREC *);
void tens6to3by3(BC_CPREC *, BC_CPREC [3][3]);
BC_CPREC quake_weight(BC_CPREC , BC_CPREC , BC_CPREC ,int );
void sum_smoothed_seismicity(struct cat *, int );
BC_CPREC quake_scale_lkm(BC_CPREC);
void print_summed_moment(struct cat *, char *);

void calc_b_value_thomas(BC_CPREC *, int , BC_CPREC , BC_CPREC ,BC_CPREC *, BC_CPREC *);
void calc_b_value_marzocci(BC_CPREC *, int , BC_CPREC , BC_CPREC ,BC_CPREC *, BC_CPREC *);
void calc_b_value_ml(BC_CPREC *, int, BC_CPREC,BC_CPREC,BC_CPREC *,BC_CPREC *);

void normalize_tens6(BC_CPREC *);
BC_CPREC std_quick(int , BC_CPREC , BC_CPREC );
void stridip(BC_CPREC ,BC_CPREC ,BC_CPREC ,BC_CPREC *,BC_CPREC *);
void ranger(BC_CPREC *);
void find_alt_plane(BC_CPREC ,BC_CPREC ,BC_CPREC ,BC_CPREC *,BC_CPREC *,BC_CPREC *);
void tensor2fpangle(BC_CPREC *, BC_CPREC *, BC_CPREC *, BC_CPREC *,
		    BC_CPREC *, BC_CPREC *, BC_CPREC *);
void print_quake_cmt_fp(FILE *, struct qke );
void print_stress_tensors(struct cat *, char *);
void make_histogram(BC_CPREC *, int , BC_CPREC , BC_CPREC *, BC_CPREC *,int **,BC_CPREC **, int *);
void print_histogram(int *, BC_CPREC *, int, FILE *);

BC_CPREC kostrov_blon(int , struct kostrov_sum *);
BC_CPREC kostrov_blat(int , struct kostrov_sum *);

void add_quake_to_bin_list(unsigned int , struct bn *,BC_CPREC);
void calc_stress_tensor_for_kbins(struct cat *);
void solve_stress_michael(int , BC_CPREC *,BC_CPREC *,BC_CPREC *, BC_CPREC *, long int *);
void sincos(double , double *, double *);
/* Andy Michael routines */
void michael_assign_to_matrix(BC_CPREC ,BC_CPREC ,BC_CPREC ,int *,BC_CPREC **,BC_CPREC **);
void michael_leasq(BC_CPREC *,int ,int ,BC_CPREC *,BC_CPREC *,BC_CPREC *,BC_CPREC *,BC_CPREC *) ;
void michael_solve_lsq(int ,int , int , BC_CPREC *, BC_CPREC *, BC_CPREC *,BC_CPREC *);

void sincosd(double , double *, double *);

/* GMT 4.5.18 and psmeca routine */
/* for fault planes */


double computed_rake2(double ,double ,double ,double ,double );
int GMT_jacobi (double *, int *, int *, double *, double *, double *, double *,  int *);

/* from meca stuff */
void GMT_momten2axe(struct M_TENSOR ,struct AXIS *,struct AXIS *,struct AXIS *);
void axe2dc(struct AXIS ,struct AXIS ,struct nodal_plane *,struct nodal_plane *);
