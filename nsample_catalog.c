#include "catalog.h"
/* 

   perform a kostrov summation and / or Michael (1984) type stress
   inversion based on AKI or gCMT style earthquake focal
   mechanism/moment tensor catalogs for spatial selection criteria 

   there is also bin_catalog which uses simple binning

   nmin, dist_max > 0, use bins with at least nmin entries within dist_max[km]
   nmin<0, dist_max> 0, use bins with at -nmin entries as long as they are within dist_max[km]

   (c) Thorsten Becker, thwbecker@post.harvard.edu, see README

*/

int main(int argc, char **argv)
{
  
  /*  */
  
  struct cat *catalog;
  struct kostrov_sum *kostrov;
  
  int itmp,ic;
  char out_filename[500],out_filename2[500],out_istring[500];
  BC_BOOLEAN use_aki = BC_TRUE;
  BC_BOOLEAN remove_trace =  BC_TRUE;	/* remove trace from summations */
  BC_BOOLEAN calc_stress =   BC_TRUE;
  BC_BOOLEAN compute_dtree = BC_TRUE;
  int  use_weights = 0; /* 0: none 1: by distance from center of bin*/
  const int monte_carlo = 0;

  catalog=(struct cat *)calloc(1,sizeof(struct cat));
  /* default parameters */
  kostrov = catalog->sum;
  kostrov_set_defaults(kostrov); /* set binning defaults */
  
  /*  */
  catalog->use_friction_solve = 2; /* 1: additional stress inversion 2: optimize friction*/
  //catalog->use_friction_solve = 1; /* 1: additional stress inversion 2: optimize friction*/
  
  snprintf(out_istring,sizeof(out_istring),"nsample");
  
  if(!catalog)
    BC_MEMERROR(argv[0]);

  /* 
     range of summation 
  */
 
  ic=2;
  if(argc < ic){
    fprintf(stderr,"%s catalog.aki [dx, %g] [min_mag, %g] [max_mag, %g] [min_lon, %g] [max_lon, %g] [min_lat, %g] [max_lat, %g] [max_depth, %g] [use_aki, %i] [weights, %i] [nmin, %i] [dist_max, %g] [nmin, %g] [is_xy, %i] [out_istring, %s] [dy, dx]\n",
	    argv[0],kostrov->dx,kostrov->minmag,kostrov->maxmag,kostrov->dlonmin, kostrov->dlonmax, kostrov->dlatmin, kostrov->dlatmax,
	    kostrov->maxdepth,(int)use_aki,use_weights,
	    kostrov->nmin,kostrov->dist_max,kostrov->mindepth,(int)catalog->is_xy,out_istring);
    exit(-1);
  }
  if(argc > ic){
    sscanf(argv[ic],"%lf",&kostrov->dx);
    kostrov->dy = kostrov->dx;
  }
  ic++;
  if(argc>ic)
    sscanf(argv[ic],"%lf",&kostrov->minmag);
  ic++;
  if(argc>ic)
    sscanf(argv[ic],"%lf",&kostrov->maxmag);
  ic++;
  if(argc>ic)sscanf(argv[ic],"%lf",&kostrov->dlonmin);
  ic++;
  if(argc>ic)sscanf(argv[ic],"%lf",&kostrov->dlonmax);
  ic++;
  if(argc>ic)sscanf(argv[ic],"%lf",&kostrov->dlatmin);
  ic++;
  if(argc>ic)sscanf(argv[ic],"%lf",&kostrov->dlatmax);
  ic++;
  if(argc>ic)sscanf(argv[ic],"%lf",&kostrov->maxdepth);
  ic++;
  if(argc>ic){
    sscanf(argv[ic],"%i",&itmp);
    use_aki = (BC_BOOLEAN)itmp;
  }
  ic++;
  if(argc>ic)
    sscanf(argv[ic],"%i",&use_weights);
  ic++;
  if(argc>ic)
    sscanf(argv[ic],"%i",&kostrov->nmin);
  ic++;
  if(argc>ic)
    sscanf(argv[ic],"%lf",&kostrov->dist_max);
  ic++;
  if(argc>ic)
    sscanf(argv[ic],"%lf",&kostrov->mindepth);
  ic++;
  if(argc>ic){
    sscanf(argv[ic],"%i",&itmp);
    catalog->is_xy = (BC_BOOLEAN)itmp;
  }
  ic++;
  if(argc>ic){
    snprintf(out_istring,sizeof(out_istring),"%s",argv[ic]);
  }
  ic++;
  if(argc>ic){
    sscanf(argv[ic],"%lf",&kostrov->dy);
  }
  ic++;
  snprintf(out_filename,sizeof(out_filename),"%s.%g.%g.%i.%g",out_istring,
	   kostrov->dx,kostrov->dy,kostrov->nmin,kostrov->dist_max);
  
  fprintf(stderr,"%s: catalog: %s dx: %g dy: %g min_mag: %g max_mag: %g min_lon: %g max_lon: %g min_lat: %g: max_lat: %g\n", 
	  argv[0],argv[1],kostrov->dx,kostrov->dy,kostrov->minmag,kostrov->maxmag,kostrov->dlonmin, kostrov->dlonmax, kostrov->dlatmin, kostrov->dlatmax);


 fprintf(stderr,"%s: min_depth: %g max_depth: %g use_aki: %i use_weights: %i nmin: %i dist_max: %g is_xy: %i out_name: %s\n", 
	 argv[0],kostrov->mindepth,kostrov->maxdepth,(int)use_aki,use_weights,
	 kostrov->nmin,kostrov->dist_max,(int)catalog->is_xy,out_filename);

  
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
 
  /* 
     
     setup bins

  */
  /* parameters */
  /*  */
  setup_kostrov(catalog,use_weights);
  /* 
     assemble based on nmin and dist_max criteria
  */
  assemble_bins_based_on_distance(catalog,remove_trace,monte_carlo,BC_TRUE);
  /* 
     print non-zero 
  */
  print_kostrov_bins(catalog,out_filename,monte_carlo);

  if(calc_stress){
    /* compute Andy Michael style stress tensors */
    snprintf(out_filename2,sizeof(out_filename2),"%s.%g.%g.%i.%g",out_istring,kostrov->dx,kostrov->dy,
	     kostrov->nmin,kostrov->dist_max);
    calc_stress_tensor_for_kbins(catalog);
    print_stress_tensors(catalog,out_filename2);
  }
  geo_tree_destroy(catalog->tree);
  return 0;
}

