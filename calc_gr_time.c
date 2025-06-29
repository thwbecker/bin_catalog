#include "catalog.h"
/* 
   
   read a set of time magnitude pairs from stdin and compute the b value,

   
  
*/
#define NGR_MODES 2

int main(int argc, char **argv)
{
  BC_CPREC dm = 0.01;		/* default precision spacing for
				   magnitude values or difference for b posiitve */
  BC_CPREC mcomplete = 2;		/* default magnitude of
					   completeness */
  /*  */
  BC_CPREC *mag,*time,b[NGR_MODES],sb,mmin,mmax,dt,tmin,tmax;
  BC_CPREC *mbin = NULL;
  int *nentry = NULL;
  BC_CPREC bin_width = 0.2;	/* for histogram */
  long nm,nm1;
  int nbin;

  int gr_mode[NGR_MODES]={1,3};		/* default modes */
  int j;
  long int ileft,iright,nuse,i;
  static int warned = 0;
  if((argc<2)||(argc>4)||(strcmp(argv[1],"-h")==0)){
    fprintf(stderr,"\ncalc_gr_time dt [dM,%g] [Mcomplete, %g] \n",
	    dm,mcomplete);
    exit(-1);
  }
  sscanf(argv[1],BC_PREC_FMT,&dt);
  if(argc>3)
    sscanf(argv[2],BC_PREC_FMT,&dm);
  if(argc>3)			/* mag completeness */
    sscanf(argv[3],BC_PREC_FMT,&mcomplete);

  nm=0;mmin=tmin=1e20;mmax=tmax=-1e20;
  mag=(BC_CPREC *)malloc(sizeof(BC_CPREC));
  time=(BC_CPREC *)malloc(sizeof(BC_CPREC));
  while(fscanf(stdin,BC_PREC2_FMT,(time+nm),(mag+nm))==2){
    if(nm)
      if(time[nm]-time[nm-1]<=0){
	fprintf(stderr,"%s: error, times should be monotonously increasing, %li %g %g\n",
		argv[0],nm,time[nm-1],time[nm]);
	exit(-1);
      }
    nm++;
    if(tmin>time[nm])
      tmin = time[nm];
    if(tmax<time[nm])
      tmax = time[nm];
    
    mag=(BC_CPREC *)realloc(mag,sizeof(BC_CPREC)*(nm+1));
    time=(BC_CPREC *)realloc(time,sizeof(BC_CPREC)*(nm+1));
    if((!mag)||(!time)){
      fprintf(stderr,"%s: mem error %li\n",argv[0],nm);
      exit(-1);
    }
  }
  nm1=nm-1;
  if(mcomplete < -10){		/* determine from histogram */
    /* compute extreme values of magnitudes and make a regular histogram
       with bin_width space */
    make_histogram(mag,nm,bin_width,&mmin,&mmax,&nentry,&mbin,&nbin);
    //print_histogram(nentry,mbin,nbin,stderr);
    
    mcomplete = max_x_from_int_vector(mbin,nentry,nbin);
    fprintf(stderr,"%s: read %li events t %g to %g, mag %g to %g, computing GR b for dt %g dM=%g Mc=%g (from max hist)\n",
	    argv[0],nm,tmin,tmax,mmin,mmax,dt,dm,mcomplete);
    
  }else{
    for(i=0;i<nm;i++){
      if(mag[i]>mmax)
	mmax = mag[i];
      if(mag[i]<mmin)
	mmin=mag[i];
    }
    fprintf(stderr,"%s: read %li events from %f to %f, mag %g to %g, computing b for dt %g dM=%g Mc=%g\n",
	    argv[0],nm,tmin,tmax,mmin,mmax,dt,dm,mcomplete);
  }
  ileft=0;
  iright=1;
  while((time[iright]-time[ileft]<dt)&&(iright<nm))
    iright++;
  nuse = iright-ileft+1;
  for(j=0;j<NGR_MODES;j++)
    calc_gr_switch((mag+ileft),nuse,dm,mcomplete,(b+j),&sb,gr_mode[j]);
  for(i=0;i<iright;i++)
    printf("%20.10e %8.5f %8.5f %6li\n",time[i],NAN,NAN,NAN);
  printf("%20.10e %8.5f %8.5f %6li\n",time[iright],b[0],b[1],nuse);

  //  fprintf(stderr,"%li %li N %li %g %g dt %g - %g\n",ileft,iright,nuse,time[ileft],time[iright],time[iright]-time[ileft],b);
  do{
    iright++;
    while((time[iright]-time[ileft] > dt)&&(ileft<iright))
      ileft++;
    nuse = iright-ileft+1;
    for(j=0;j<NGR_MODES;j++)
      calc_gr_switch((mag+ileft),nuse,dm,mcomplete,(b+j),&sb,gr_mode[j]);
    printf("%20.10e %8.5f %8.5f %6li\n",time[iright],b[0],b[1],nuse);
    //fprintf(stderr,"%li %li N %li %g %g dt %g - %g\n",ileft,iright,nuse,time[ileft],time[iright],time[iright]-time[ileft],b);
  }while(iright<nm1);
  free(mag);free(time);
}

void calc_gr_switch(BC_CPREC *mag, long nm, BC_CPREC dm,BC_CPREC mcomplete,BC_CPREC *b,BC_CPREC *sb,int mode)
{
  /* modes */
  /* 0: Aki ML */
  /* 1: Utsu (1966), Bender (1983) dM corrected b value, Shi and Bolt (1982) sb approach  */
  /* 2: Marzocci */
  /* 3: b-positive */
     
  static int warned=0;
  const int nmin=50;		/* warn if smaller  */
  if(!warned){
    if(nm < nmin){
      fprintf(stderr,"calc_gr_switch at least one bin with N %li smaller than desired Nmin %i\n",
	      nm,nmin);
      warned=1;
    }
  }
  switch(mode){
  case 0:					 /* not the best idea */
    calc_b_value_ml(mag,nm,mcomplete,b,sb); /* simple Aki/ML method
						 wihtout dM
						 correction, old
						 errors */
    break;
  case 1:			/* pretty good idea */
    /* Utsu (1966), Bender (1983) dM corrected b value, Shi and Bolt
       (1982) sb approach */
    calc_b_value_thomas(mag,nm,dm,mcomplete,b,sb);
    break;
  case 2:			/* Marzocci */
    calc_b_value_marzocci(mag,nm,dm,mcomplete,b,sb);
    break;
  case 3:			/* b positive */
    calc_b_value_bpos(mag,nm,dm,b);
    *sb = NAN;
    break;
  default:
    fprintf(stderr,"calc_gr_swithc: mode %i undefined\n",mode);
    exit(-1);
    break;
  }
}
