#include "catalog.h"
/* 
   
   read a set of magnitudes from stdin and compute the b value,
   assuming some dm spacing of the magnitudes and magnitude of completeness
   
   calc_gr [dM,0.01] [Mcomplete, 2] [gr_mode, 1]

   if Mcomplete is smaller than -10, will try to determine best fit
   Mcomplete
*/

int main(int argc, char **argv)
{
  BC_CPREC dm = 0.01;		/* default precision spacing for
				   magnitude values */
  BC_CPREC mcomplete = 2;		/* default magnitude of
				   completeness */
  /*  */
  BC_CPREC *mag,b,sb,mmin,mmax;
  BC_CPREC *mbin = NULL;
  int *nentry = NULL;
  BC_CPREC bin_width = 0.2;	/* for histogram */
  int nm,nbin;
  int gr_mode = 1 ;		/* default mode */
  if(argc>1)			/* magnitude spacing */
    sscanf(argv[1],BC_PREC_FMT,&dm);
  if(argc>2)			/* mag completeness */
    sscanf(argv[2],BC_PREC_FMT,&mcomplete);
  if(argc>3)			/* mag completeness */
    sscanf(argv[3],"%i",&gr_mode);
  
  nm=0;mmin=1e20;mmax=-1e20;
  mag=(BC_CPREC *)malloc(sizeof(BC_CPREC));
  while(fscanf(stdin,BC_PREC_FMT,(mag+nm))==1){
   
    nm++;
    mag=(BC_CPREC *)realloc(mag,sizeof(BC_CPREC)*(nm+1));
  }
  /* compute extreme values of magnitudes and make a regular histogram
     with bin_width space */
  make_histogram(mag,nm,bin_width,&mmin,&mmax,&nentry,&mbin,&nbin);
  //print_histogram(nentry,mbin,nbin,stderr);
  
  if(mcomplete < -10){		/* determine from histogram */
    mcomplete = max_x_from_int_vector(mbin,nentry,nbin);
    fprintf(stderr,"%s: read %i events from %f to %f, computing GR b for dM=%g Mc=%g (from max hist), gr_mode: %i\n",
	    argv[0],nm,mmin,mmax,dm,mcomplete,gr_mode);
    
  }else{
    fprintf(stderr,"%s: read %i events from %f to %f, computing GR b for dM=%g Mc=%g, gr_mode: %i\n",
	    argv[0],nm,mmin,mmax,dm,mcomplete,gr_mode);
  }
  switch(gr_mode){
  case 0:					 /* not the best idea */
    calc_b_value_ml(mag,nm,dm,mcomplete,&b,&sb); /* simple Aki/ML
						    method wihtout dM
						    correction, old
						    errors */
    fprintf(stdout,"%7.5f %7.5f\n",b,sb);
    break;
  case 1:			/* pretty good idea */
    /* Utsu (1966), Bender (1983) dM corrected b value, Shi and Bolt
       (1982) sb approach */
    calc_b_value_thomas(mag,nm,dm,mcomplete,&b,&sb);
    fprintf(stdout,"%7.5f %7.5f\n",b,sb);
    break;
  case 2:
    calc_b_value_marzocci(mag,nm,dm,mcomplete,&b,&sb);
    fprintf(stdout,"%7.5f %7.5f\n",b,sb);
    break;

  case 3:
    calc_b_value_ml(mag,nm,dm,mcomplete,&b,&sb);fprintf(stdout,         "Aki:    %7.5f %7.5f\n",b,sb);
    calc_b_value_thomas(mag,nm,dm,mcomplete,&b,&sb);fprintf(stdout,  "Thomas: %7.5f %7.5f\n",b,sb);
    calc_b_value_marzocci(mag,nm,dm,mcomplete,&b,&sb);fprintf(stdout,"Marz:   %7.5f %7.5f\n",b,sb);
    break;
  default:
    fprintf(stderr,"%s: mode %i undefined\n",argv[0],gr_mode);
    exit(-1);
    break;
  }
  
  free(mag);
  return 0;
}
