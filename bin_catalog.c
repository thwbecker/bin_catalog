#include "catalog.h"
/* 

   perform a kostrov summation and / or Michael (1984) type stress
   inversion based on AKI or gCMT style earthquake focal
   mechanism/moment tensor catalogs for simple binning

   there is also nsample_catalog which tries to find the closest
   events and has more flexibility in terms of sampling

   (c) 2020 - 2025, Thorsten Becker, thwbecker@post.harvard.edu, see README
*/

int main(int argc, char **argv)
{
  struct cat *catalog;
  struct kostrov_sum *kostrov;
  int itmp,ic;
  char out_filename[500],out_filename2[500],out_istring[500];
  int monte_carlo = 0;
  BC_BOOLEAN use_aki     =   BC_TRUE;
  BC_BOOLEAN remove_trace =  BC_TRUE;	/* remove trace from summations */
  BC_BOOLEAN calc_stress =   BC_TRUE;

  BC_BOOLEAN compute_dtree = BC_FALSE; /* not needed */
  /* default */
  int min_events_for_stress = 5; /* at least so many events to attempt a stress inversion */

  int  weighting_method = 0; /* 0: none 
				1: distance from center of time 
				2: normalized by numbers over time 
			     */
  catalog=(struct cat *)calloc(1,sizeof(struct cat));
  
  catalog->use_friction_solve = 2; /* 1: additional stress inversion 2: optimize friction*/
  //catalog->use_friction_solve = 1; /* 1: additional stress inversion 2: optimize friction*/
  sprintf(out_istring,"kostrov");
  
  if(!catalog)
    BC_MEMERROR(argv[0]);
  
  kostrov = catalog->sum;
  kostrov_set_defaults(kostrov); /* set binning defaults */
  
  
  ic=2;
  if(argc < ic){
    fprintf(stderr,"%s catalog.aki [dx, %g] [min_mag, %g] [max_mag, %g] [monte_carlo, %i] [min_lon, %g] [max_lon, %g] [min_lat, %g] [max_lat, %g] [max_depth, %g] [use_aki, %i] [weighting_method (0/1/2), %i] [min_depth, %g] [is_xy, %i] [out_istring, %s] [dy, dx] [nmin_for_stress, %i]\n",
	    argv[0],kostrov->dx,kostrov->minmag,kostrov->maxmag,monte_carlo,
	    kostrov->dlonmin, kostrov->dlonmax, 
	    kostrov->dlatmin, kostrov->dlatmax,
	    kostrov->maxdepth,(int)use_aki,weighting_method,
	    kostrov->mindepth,(int)catalog->is_xy,out_istring,min_events_for_stress);
    exit(-1);
  }
  if(argc>ic){
    sscanf(argv[ic],BC_PREC_FMT,&kostrov->dx);
    kostrov->dy = kostrov->dx;
  }
  ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->minmag);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->maxmag);ic++;
  if(argc>ic)sscanf(argv[ic],"%i",       &monte_carlo);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->dlonmin);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->dlonmax);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->dlatmin);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->dlatmax);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->maxdepth);ic++;
  if(argc>ic){
    sscanf(argv[ic],"%i",&itmp);use_aki = (BC_BOOLEAN)itmp;
  }
  ic++;
  if(argc>ic)sscanf(argv[ic],"%i",&weighting_method);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->mindepth);ic++;
  if(argc>ic){
    sscanf(argv[ic],"%i",&itmp);
    catalog->is_xy = (BC_BOOLEAN)itmp;
  }
  ic++;
  if(argc>ic)sprintf(out_istring,"%s",argv[ic]);ic++;
  if(argc>ic)sscanf(argv[ic],BC_PREC_FMT,&kostrov->dy);ic++;
  if(argc>ic)sscanf(argv[ic],"%i",&min_events_for_stress);ic++;
  /*  */
  kostrov->nmin = min_events_for_stress;
  /* output filename for Kostrov summations */
  snprintf(out_filename,sizeof(out_filename),"%s.%g.%g.%i",out_istring,kostrov->dx,kostrov->dy,monte_carlo);
  
  if(use_aki){			/* aki with last column time */
    fprintf(stderr,"%s: assuming AKI format (last column is UNIX time)\n",argv[0]);
    read_catalog(argv[1],catalog,BC_AKI,compute_dtree);
  }else{
    fprintf(stderr,"%s: assuming CMT format\n",argv[0]);
    read_catalog(argv[1],catalog,BC_CMT,compute_dtree);
  }
  if(!catalog->n){
    fprintf(stderr,"%s: error: zero events read\n",argv[0]);
    exit(-1);
  }
  switch(weighting_method){
  case 0:
    break;
  case 1:
    fprintf(stderr,"%s: weighing sum by temporal distance from tsec %g\n",
	    argv[0],catalog->tcenter);
    break;
  case 2:
    fprintf(stderr,"%s: frequency normalized weight\n",argv[0]);
    break;
  default:
    fprintf(stderr,"%s: normalizing method %i undefined\n",argv[0],weighting_method);
    exit(-1);
    break;
  }
  
  /* 
     
     setup bins
  */
  /* parameters */

  setup_kostrov(catalog,weighting_method);
  
  /* 
     sum 
  */
  sum_kostrov_bins(catalog,remove_trace,monte_carlo,BC_TRUE);
  /* 
     print non-zero 
  */
  print_kostrov_bins(catalog,out_filename,monte_carlo);
  if(calc_stress){
    /* compute Andy Michael style stress tensors */
    snprintf(out_filename2,sizeof(out_filename2),"%s.%g.%g",out_istring,kostrov->dx,kostrov->dy);
    calc_stress_tensor_for_kbins(catalog);
    print_stress_tensors(catalog,out_filename2);
  }
  return 0;
}

