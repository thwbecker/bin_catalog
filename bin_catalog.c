#include "catalog.h"
/* 

   perform a kostrov summation 

*/

int main(int argc, char **argv)
{
  double dx = 0.2,min_mag=2.49,max_mag=8.01;
  double lonmin,lonmax,latmin,latmax;
  double mindepth,maxdepth;
  double dy;
  
  struct cat *catalog;
  int itmp;
  char out_filename[500],out_filename2[500],out_istring[500];
  BC_BOOLEAN  monte_carlo = BC_FALSE, use_aki = BC_TRUE;
  BC_BOOLEAN remove_trace = BC_TRUE;	/* remove trace from summations */
  BC_BOOLEAN  calc_stress = BC_TRUE;
  int  weighting_method = 0; /* 0: none 
				1: distance from center of time 
				2: normalized by numbers over time 
			     */
  dy = dx;
  catalog=(struct cat *)malloc(sizeof(struct cat)); 
  sprintf(out_istring,"kostrov");
  
  if(!catalog)
    BC_MEMERROR(argv[0]);

  /* 
     range of summation 
  */
  lonmin = 232;lonmax=250;
  latmin=30;latmax =45;
  mindepth = -10;maxdepth = 15;	/* in km */
  
  if(argc < 2){
    fprintf(stderr,"%s catalog.aki [dx, %g] [min_mag, %g] [max_mag, %g] [monte_carlo, %i] [min_lon, %g] [max_lon, %g] [min_lat, %g] [max_lat, %g] [max_depth, %g] [use_aki, %i] [weighting_method (0/1/2), %i] [dmin, %g] [is_xy, %i] [out_istring, %s] [dy, dx]\n",
	    argv[0],dx,min_mag,max_mag,(int)monte_carlo,
	    lonmin, lonmax, 
	    latmin, latmax,
	    maxdepth,(int)use_aki,weighting_method,
	    mindepth,(int)catalog->is_xy,out_istring);
    exit(-1);
  }
  if(argc>2){
    sscanf(argv[2],"%lf",&dx);
    dy = dx;
  }
  if(argc>3)
    sscanf(argv[3],"%lf",&min_mag);
  if(argc>4)
    sscanf(argv[4],"%lf",&max_mag);
  if(argc>5){
    sscanf(argv[5],"%i",&itmp);
    monte_carlo = (BC_BOOLEAN)itmp;
  }
  if(argc>6)sscanf(argv[6],"%lf",&lonmin);
  if(argc>7)sscanf(argv[7],"%lf",&lonmax);
  if(argc>8)sscanf(argv[8],"%lf",&latmin);
  if(argc>9)sscanf(argv[9],"%lf",&latmax);
  if(argc>10)sscanf(argv[10],"%lf",&maxdepth);
  if(argc>11){
    sscanf(argv[11],"%i",&itmp);
    use_aki = (BC_BOOLEAN)itmp;
  }
  if(argc>12){
    sscanf(argv[12],"%i",&weighting_method);
  }
  if(argc>13)sscanf(argv[13],"%lf",&mindepth);
  if(argc>14){
    sscanf(argv[14],"%i",&itmp);
    catalog->is_xy = (BC_BOOLEAN)itmp;
  }
  if(argc>15){
    sprintf(out_istring,"%s",argv[15]);
  }
  if(argc>16){
    sscanf(argv[16],"%lf",&dy);
  }
  /* output filename for Kostrov summations */
  sprintf(out_filename,"%s.%g.%g.%i",out_istring,dx,dy,(int)monte_carlo);
  
  if(use_aki){			/* aki with last column time */
    fprintf(stderr,"%s: assuming AKI format (last column is UNIX time)\n",argv[0]);
    read_catalog(argv[1],catalog,BC_AKI);
  }else{
    fprintf(stderr,"%s: assuming CMT format\n",argv[0]);
    read_catalog(argv[1],catalog,BC_CMT);
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
  /*  */
  setup_kostrov(catalog,lonmin,lonmax,latmin,latmax,mindepth,maxdepth,
		dx,dy,min_mag,max_mag,weighting_method);  
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
    sprintf(out_filename2,"%s.%g.%g",out_istring,dx,dy);
    calc_stress_tensor_for_kbins(catalog);
    print_stress_tensors(catalog,out_filename2);
  }
  return 0;
}

