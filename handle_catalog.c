#include "catalog.h"

/* 

relocate catalog a with locations from catalog b and save to catalog c

*/
void relocate_catalog(struct cat *a,struct cat *b, struct cat *c)
{
  int j,yes,no,found,i,is_xy_sum;
  BC_CPREC hor_dist, temp_dist, vert_dist, mag_dist,
    hor_dist0,temp_dist0,vert_dist0,mag_dist0,
    hor_distm,temp_distm,vert_distm,mag_distm;
  is_xy_sum = (int)a->is_xy + (int)b->is_xy + (int)c->is_xy;
  if((is_xy_sum != 3)&&(is_xy_sum != 0)){
    fprintf(stderr,"relocate_catalog: catalog coordinate system inconsistent: %i/%i/%i\n",
	    (int)a->is_xy, (int)b->is_xy, (int)c->is_xy);
    exit(-1);
  }
  hor_distm = temp_distm = vert_distm = mag_distm = 0.0;
  yes = no = 0;
  for(i=0;i < a->n;i++){
    if(i % 500 == 0)
      fprintf(stderr,
	      "relocate_catalog: checking doubles: %7i out of %7i tests, reassigned %6i missed %6i\r",
	      i,a->n,yes,no);
    
    temp_dist0 = 200;		/* seconds */
    hor_dist0 = 150;		/* km horizontal */
    vert_dist0 = 100;		/* km vertical */
    mag_dist0 = 1.5;		/* magnitude distance */
    for(j=0, found = -1;j < b->n;j++){
      temp_dist = fabs(a->quake[i].tsec - b->quake[j].tsec);
      if(temp_dist < temp_dist0){	
	hor_dist = distance(a,b,i,j);
	vert_dist = fabs(a->quake[i].depth - b->quake[j].depth);
	//fprintf(stderr,"%g %g\n",a->quake[i].mag, b->quake[j].mag);
	mag_dist = fabs(a->quake[i].mag - b->quake[j].mag);
	
	if((hor_dist < hor_dist0) && 	/* horizontal distance in km */
	   (vert_dist < vert_dist0) && 
	   (mag_dist < mag_dist0)){ /* vertical  distance in km */

	  temp_dist0 = temp_dist;
	  mag_dist0  = mag_dist;
	  vert_dist0 = vert_dist;
	  hor_dist0  = hor_dist;

	  found = j;
	}
      }
    }
    copy_quake((a->quake+i),(c->quake+c->n));
    if(found < 0)
      no++;
    else{
      /* use the new locations */
      yes++;
      
      hor_distm  += hor_dist0;
      vert_distm += vert_dist0;
      temp_distm += temp_dist0;
      mag_distm  += mag_dist0;
      /* hmm */
      c->quake[c->n].dlon = b->quake[found].dlon; /* degree */
      c->quake[c->n].dlat = b->quake[found].dlat;

      c->quake[c->n].plon = b->quake[found].plon; /* degree */
      c->quake[c->n].plat = b->quake[found].plat;

      c->quake[c->n].lon = b->quake[found].lon; /* radian */
      c->quake[c->n].lat = b->quake[found].lat;
      
      c->quake[c->n].depth = b->quake[found].depth;
    }
    c->n += 1;
    make_room_for_quake(c);
  }
  fprintf(stderr,"\nrelocate_catalog: found %i new locations, %i old, %i new catalog total, %i old\n",
	  yes,no,c->n,a->n);
  hor_distm /= (BC_CPREC)yes;
  vert_distm /= (BC_CPREC)yes;
  temp_distm /= (BC_CPREC)yes;
  mag_distm /= (BC_CPREC)yes;
  fprintf(stderr,"relocate_catalog: mean h: %g v: %g t: %g mag: %g\n",
	  hor_distm,vert_distm,temp_distm,mag_distm);
}



/* 

   given an intialized set of bins and catalogs, sum over all quakes
   and assign weights to each quake for summing (and possible stress
   analysis)

*/
void sum_kostrov_bins(struct cat *catalog, BC_BOOLEAN do_remove_trace,
		      BC_BOOLEAN monte_carlo,BC_BOOLEAN verbose)
{
  int i,j,ix,iy,bind,n,nmonte,bc,ibin,iquake,jbq;
  BC_CPREC eps_angle = BC_D2R(BC_EPS_ANGLE_FOR_RANDOM_DEG); /* sigma for random */
  BC_CPREC pin[6],pout[6],me,dx,weight,tfracl,tmin,tmax,trange,dt,navg;
  int ndt,nndt_nz;		/* for time weighted binning */

  BC_CPREC dt0 = (60*60*24*30)*2;		/* for weighting, two months */

  int *nbc;
  char tbuf1[80];
  struct tm ts;
  time_t timet;
  struct kostrov_sum *kostrov;
  if(!catalog->sum->init){
    fprintf(stderr,"sum_kostrov_bins: sum bounds not initialized\n");
    exit(-1);
  }
  kostrov = catalog->sum;
  if(verbose)
    fprintf(stderr,"sum_kostrov_bins: summing with %g < mag < %g, max depth: %g, weighting method: %i\n",
	    kostrov->minmag,kostrov->maxmag,
	    kostrov->maxdepth,
	    kostrov->weighting_method);
  if(monte_carlo){		/* monte carlo */
    /* number of monte_carlos */
    nmonte = monte_carlo;
  }else{
    nmonte = 1;
  }
  if(!catalog->is_xy){		/* if not carteisan, move to 0...360 */
    for(i=0;i < catalog->n;i++)
      ranger(&(catalog->quake[i].dlon));
  }
  
  /* only clear once */
  clear_bins(catalog);		
  for(bc = 0;bc < nmonte;bc++){		/* loop through all monte_carlo sampling */
    /* 
       clear local pn 
    */
    for(ibin=0;ibin < kostrov->nxny;ibin++)
      for(j=0;j<6;j++)
	kostrov->bin[ibin].mnloc[j]=0.0;
     n=0;				/* how many are added  */
    for(i=0;i < catalog->n;i++){	/* loop through quakes */
      /* make sure within magnitude and depth range range */
      if((catalog->quake[i].mag >= kostrov->minmag) && 
	 (catalog->quake[i].mag <= kostrov->maxmag) && 
	 (catalog->quake[i].depth <= kostrov->maxdepth) && 
	 (catalog->quake[i].depth >= kostrov->mindepth)){
	/* 
	   move to grid-centered 
	*/
	dx = catalog->quake[i].dlon - kostrov->dlonmin;
	if(!catalog->is_xy){
	  if(dx > 360)
	    dx -= 360;
	}
	ix = (int)(dx/kostrov->ddlon + 0.5);
	iy = (int)((catalog->quake[i].dlat - kostrov->dlatmin)/kostrov->ddlat + 0.5);
	//fprintf(stderr,"%i %i\n",ix,iy);
	if((ix >= 0) && (ix < kostrov->nx) && (iy >= 0) && (iy < kostrov->ny)){ /* make sure within region */
	  /* 
	     within region
	     
	  */
	  bind = ix * kostrov->ny + iy; /* bin index */
	  if(kostrov->weighting_method == 1){
	    weight = 1 - fabs(catalog->quake[i].tsec - catalog->tcenter)/catalog->trange*2;
	    //fprintf(stderr,"w %g t %g tc %g tr %g\n",weight, catalog->quake[i].tsec, catalog->tcenter,catalog->trange);
	  }else{
	    weight = 1.0;
	  }
	  /* 
	     normalized 
	  */
	  if(bc == 0){
	    /* 
	       
	       only for regular/first time around, assign the quake
	       number and the quake's weight in the bin
	       
	    */
	    add_quake_to_bin_list((unsigned int)i,(kostrov->bin+bind),weight);
	    
	    kostrov->bin[bind].sumw += weight;	/* sum weights, increment */
	    /*  */
	    kostrov->mtot += catalog->quake[i].m0*weight; /* total potency */
	    for(j=0;j<6;j++){
	      kostrov->bin[bind].m[j]  += (catalog->quake[i].m[j] * weight * catalog->quake[i].m0); /* scaled sum */
	      kostrov->bin[bind].mnloc[j] += weight * catalog->quake[i].m[j]; /* normalized sum */
	    }
	  }else{
	    /* monte carlo realization for normalized */
	    for(j=0;j<6;j++)
	      pin[j] = catalog->quake[i].m[j]; /* rotate tensor */
	    rotate_vec6(pin,pout, gauss_ran(&catalog->seed,eps_angle),
			gauss_ran(&catalog->seed,eps_angle),gauss_ran(&catalog->seed,eps_angle));
	    for(j=0;j<6;j++)
	      kostrov->bin[bind].mnloc[j] += weight * pout[j];
	  }
	  n++;			/* increment total number of assigned quakes */
	}
      }
    } /* end data loop */
    if(verbose){
      if(bc == 0){
	fprintf(stderr,"sum_kostrov_bins: used %i out of %i events, summed moment: %.3e\n",
		n,catalog->n,kostrov->mtot);
      }else{
	fprintf(stderr,"sum_kostrov_bins: monte_carlo # %06i, used %06i \r",bc,n);
      }
    }
    if((bc==0)&&(kostrov->weighting_method==2)){
      /* 
	 
	 re-weigh over time

      */
      
      if(nmonte>1){
	fprintf(stderr,"sum_kostrov_bins: nmonte %i not compatible with weighting %i\n",
		nmonte,kostrov->weighting_method);
	exit(-1);
      }
      fprintf(stderr,"sum_kostrov_bins: reassigning weights for each quake given temporal sampling in bin\n");
      /* reset total weighting */
      kostrov->mtot = 0;
      for(ibin=0;ibin<kostrov->nxny;ibin++){ /* bin loop */
	/* reset bins  */
	if(kostrov->bin[ibin].n){
	  /* reset weights */
	  kostrov->bin[ibin].sumw = 0;
	  /* reset binned moments */
	  for(j=0;j < 6;j++){
	    kostrov->bin[ibin].m[j]=0;
	    kostrov->bin[ibin].mnloc[j]=0;
	  }
	  /* 
	     find time range covered in this bin
	  */
	  tmin = 1e20;tmax=-1e20;
	  for(jbq=0;jbq < kostrov->bin[ibin].n;jbq++){
	    iquake = kostrov->bin[ibin].quake[jbq];
	    if(catalog->quake[iquake].tsec > tmax)
	      tmax = catalog->quake[iquake].tsec;
	    if(catalog->quake[iquake].tsec < tmin)
	      tmin = catalog->quake[iquake].tsec;
	  }
	  trange = tmax - tmin;	/* range within this bin */
	  dt = dt0;		/* make the time interval one month if
				   possible */
	  if(trange < dt)
	    dt = trange;
	  if(dt == 0)
	    dt = dt0;
	  ndt = (int)(trange/dt) + 1;
	  nbc = (int *)malloc(ndt*sizeof(int));
	  for(i=0;i < ndt;i++)
	    nbc[i]=0;
	  if(trange > 0){	/* bin the event number over time */
	    for(jbq=0;jbq < kostrov->bin[ibin].n;jbq++){ /* bin over time */
	      iquake = kostrov->bin[ibin].quake[jbq];
	      i = (int)((catalog->quake[iquake].tsec-tmin)/dt);
	      if(i == ndt)
		i--;
	      nbc[i]++;
	    }
	    navg=0.;nndt_nz = 0;
	    for(i=0;i<ndt;i++)
	      if(nbc[i]){
		navg += (BC_CPREC)nbc[i];
		nndt_nz++;
	      }
	    navg /= (BC_CPREC)nndt_nz;
	  }else{
	    navg = 1;
	    nbc[0] = 1;

	  }
	  for(jbq=0;jbq < kostrov->bin[ibin].n;jbq++){
	    /* now redo the summation with the time-weighting */
	    iquake = kostrov->bin[ibin].quake[jbq];
	    /* obtain actual calendar time for debugging */
	    timet=(time_t)catalog->quake[iquake].tsec;ts = *localtime(&timet);
	    strftime(tbuf1, sizeof(tbuf1), "%a %Y-%m-%d %H:%M:%S %Z", &ts);
	    if(trange > 0){	/* more than one event */
	      tfracl =  (catalog->quake[iquake].tsec-tmin)/trange; /* local time  */
	      i = (int)((catalog->quake[iquake].tsec-tmin)/dt); /* which bin */
	      if(i==ndt)
		i--;
	      weight = navg/(BC_CPREC)nbc[i]/(BC_CPREC)kostrov->bin[ibin].n; /* scale with the inverse
									  of the number of
									  events */
	    }else{
	      i=0;
	      weight = 1;
	      tfracl = 1;
	    }
	    if(0)
	      fprintf(stderr,"%5i - %6i - %6.4f %3i/%3i - %3i/%6.1f - %6.4f - %s\n",
		      ibin+1,iquake+1,tfracl,i,ndt,nbc[i],navg,weight,tbuf1);
	    /*  */
	    kostrov->bin[ibin].weight[jbq] = weight; /* reassign the weight of this
							quake in the binning */
	    kostrov->bin[ibin].sumw += weight;	/* sum weights for bin, should total unity */
	    /*  */
	    kostrov->mtot += catalog->quake[iquake].m0 * weight; /* total potency */
	    for(j=0;j<6;j++){
	      kostrov->bin[ibin].m[j]  += (catalog->quake[iquake].m[j] * weight * catalog->quake[iquake].m0); /* scaled sum */
	      kostrov->bin[ibin].mnloc[j] += weight * catalog->quake[iquake].m[j]; /* normalized sum */
	    }
	  }
	  free(nbc);
	  //fprintf(stderr,"%3i %11g\n\n",ibin,kostrov->bin[ibin].sumw);
	}
      }
    } /* end time reweighting */
    
    /* 
       take mean and sum up 
    */
    for(i=0;i < kostrov->nxny;i++)
      if(kostrov->bin[i].n){
	//fprintf(stderr,"bin %i: ",i);
	//for(j=0;j < kostrov->bin[i].n;j++)fprintf(stderr,"%i ",kostrov->bin[i].quake[j]);
	//fprintf(stderr,"\n");
	for(j=0;j < 6;j++){
	  /* normalized tensor components */
	  kostrov->bin[i].mnloc[j]  /= kostrov->bin[i].sumw;
	  /* 
	     additions for mean and std of full tensor 
	  */
	  kostrov->bin[i].mn[j]  += kostrov->bin[i].mnloc[j];
	  kostrov->bin[i].smn[j] += kostrov->bin[i].mnloc[j]*kostrov->bin[i].mnloc[j];
	}
	/* addition for mean and std of mean hor strain */
	me = mean_hor_strain(kostrov->bin[i].mnloc);
	kostrov->bin[i].men  += me ;
	kostrov->bin[i].mens += me * me;
      }
  } /* end boot strap loops */
  if(verbose)
    if(bc)
      fprintf(stderr,"\n");
  
  /* 
     take mean and std 
  */
  for(i=0;i < kostrov->nxny;i++)
    if(kostrov->bin[i].n){
      /* get mean for scaled sum */
      kostrov->bin[i].me = mean_hor_strain(kostrov->bin[i].m);
      /* 
	 mean and std for mean strain 
      */
      kostrov->bin[i].mens = std_quick(bc,kostrov->bin[i].men,kostrov->bin[i].mens);
      /* mean */
      kostrov->bin[i].men = kostrov->bin[i].men/(BC_CPREC)bc;
      for(j=0;j < 6;j++){
	/* normalized tensor */
	/* standard deviation, the quick way */
	kostrov->bin[i].smn[j] = std_quick(bc, kostrov->bin[i].mn[j],kostrov->bin[i].smn[j]);
	/* mean */
	kostrov->bin[i].mn[j] = kostrov->bin[i].mn[j]/(BC_CPREC)bc;
      }
    }

  
  if(do_remove_trace){
    /* remove trace for all binx */
    fprintf(stderr,"sum_kostrov_bins: removing trace\n");
    for(i=0;i < kostrov->nxny;i++){
      remove_trace(kostrov->bin[i].m);
      remove_trace(kostrov->bin[i].mn);
      remove_trace(kostrov->bin[i].smn);
    }
  }
}
/* 

   given an intialized set of bins and catalogs, sum over all quakes
   given distance kernel 

*/
void sum_smoothed_seismicity(struct cat *catalog, int mode)
{
  int i,j;
  BC_CPREC weight,dist,tweight,wm0;
  struct kostrov_sum *kostrov;
  if(!catalog->sum->init){
    fprintf(stderr,"sum_smoothed_seismicity: sum bounds not initialized\n");
    exit(-1);
  }
  kostrov = catalog->sum;

  fprintf(stderr,"sum_smoothed_seismicity: summing with %g < mag < %g, max depth: %g, mode: %i\n",
	  kostrov->minmag,kostrov->maxmag,
	  kostrov->maxdepth,mode);
  /* only clear once */
  clear_bins(catalog);		

  /* 
     clear 
  */
  for(i=0;i < kostrov->nxny;i++)
    kostrov->bin[i].me = 0.0;
  
  for(i=0;i < catalog->n;i++){	/* loop through quakes */
    /* 
       make sure within magnitude and depth range range 
    */
    tweight = 0;
    if((catalog->quake[i].mag >= kostrov->minmag) && 
       (catalog->quake[i].mag <= kostrov->maxmag) && 
       (catalog->quake[i].depth <= kostrov->maxdepth) && 
       (catalog->quake[i].depth >= kostrov->mindepth)){
      

      
      for(j=0;j < kostrov->nxny;j++){ /* loop through all bins */
  
	if(catalog->is_xy){	/* cartesian, lon and lat are x and y
				   actually */
	  dist = distance_cart(catalog->quake[i].dlon,catalog->quake[i].dlat,
			       kostrov->bin[j].dlon,kostrov->bin[j].dlat);
	}else{
	  dist = distance_geo(catalog->quake[i].lon,catalog->quake[i].lat,
			      kostrov->bin[j].lon,kostrov->bin[j].lat,
			      catalog->quake[i].coslat,kostrov->bin[j].coslat);
	}

	  
	/* use equivalent length of quake for scale for moment */
	weight = quake_weight(catalog->quake[i].m0,catalog->quake[i].lkm,
			      dist,mode);
	wm0 = 0;
	switch(mode){
	case 0:
	  wm0 = catalog->quake[i].m0*weight/1e15;
	  break;
	case 1:
	  wm0 = sqrt(catalog->quake[i].m0)*weight/1e7; /* benioff strain */
	  break;
	case 2:
	  wm0 = weight; 
	  break;
	default:
	  fprintf(stderr,"mode %i undefined\n",mode);
	  break;
	}
	//fprintf(stderr,"%i %i %e %g %g\n",i,j,catalog->quake[i].m0,dist,weight);
	kostrov->mtot += wm0; 
	kostrov->bin[j].me  += wm0; /* scaled sum */
	tweight += weight;
      }
    } /* use this quake */
    fprintf(stderr,"quake %8i, tweight: %15.7f\r",i,tweight);
  }   /* loop through all quakes */
  fprintf(stderr,"\n");
}

/* pass 
   x[n] values
   bin_width
   nbin, xbin as NULL pointer

   returns:

   xmin,xmax
   nentry entries in each bin
   xbin mid value of each bin
   nbin number of bins

*/
void make_histogram(BC_CPREC *x, int n, BC_CPREC bin_width, BC_CPREC *xmin, BC_CPREC *xmax,
		    int **nentry,BC_CPREC **xbin, int *nbin)
{
  int i,j;
  BC_CPREC xbmin;
  *xmin = 1e20;*xmax=-1e20;
  for(i=0;i<n;i++){		/* find extrema */
    if(x[i] > *xmax)
      *xmax = x[i];
    if(x[i] < *xmin)
      *xmin = x[i];
  }
  *nbin = (int)((*xmax- *xmin)/bin_width) + 1;	/* number of bins */

  *nentry  = (int *)realloc(*nentry,sizeof(int)*(*nbin));
  *xbin  = (BC_CPREC *)realloc(*xbin,sizeof(BC_CPREC)*(*nbin));

  xbmin = ((int)(*xmin/bin_width+0.5)+0.5)*bin_width;
  
  for(i=0;i < (*nbin);i++){		/* initialize */
    *(*nentry+i) = 0;
    *(*xbin+i)   = xbmin + bin_width *i;
  }
  fprintf(stderr,"make_histogram: xmin: %g xbmin: %g xmax: %g bin_width: %g nbin: %i\n",
	  *xmin,xbmin,*xmax,bin_width,*nbin);
 
  /* bin */
  for(i=0;i < n;i++){
    j = (int)((x[i]- *xmin)/bin_width);
    *(*nentry+j) += 1;
  }
}

void print_histogram(int *nentry, BC_CPREC *xbin, int nbin, FILE *out_stream)
{
  int i,ntotal;
  for(i=ntotal=0;i<nbin;i++)
    ntotal += nentry[i];
  for(i=0;i < nbin;i++)
    fprintf(out_stream,"%12g %12i %8.6f %12g\n",
	    xbin[i],nentry[i],(BC_CPREC)nentry[i]/(BC_CPREC)ntotal,log10((BC_CPREC)nentry[i]));
	    
}



/* 

   call after initializing lon, lat range and sum->dlon/sum->dlat ! 

*/
void setup_kostrov(struct cat *catalog,
		   BC_CPREC dlonmin,BC_CPREC dlonmax, /* in degrees */
		   BC_CPREC dlatmin,BC_CPREC dlatmax,
		   BC_CPREC mindepth,BC_CPREC maxdepth,
		   BC_CPREC ddlon, BC_CPREC ddlat,
		   BC_CPREC minmag,BC_CPREC maxmag,
		   int weighting_method)
{
  int i,j,ind;
  double xmin,ymin,darea;
  struct kostrov_sum *kostrov; 
  kostrov = catalog->sum;

  kostrov->mtot = 0.0;		/* total moment of summation */
  kostrov->weighting_method = weighting_method;
  
  kostrov->dlonmin = dlonmin;
  kostrov->dlonmax = dlonmax;
  kostrov->dlatmin = dlatmin;
  kostrov->dlatmax = dlatmax;

  kostrov->mindepth = mindepth;
  kostrov->maxdepth = maxdepth;

  kostrov->minmag = minmag;
  kostrov->maxmag = maxmag;
  
  kostrov->ddlon = ddlon;
  kostrov->ddlat = ddlat;
  
  kostrov->nx = (kostrov->dlonmax - kostrov->dlonmin)/
    kostrov->ddlon;
  kostrov->ny = (kostrov->dlatmax - kostrov->dlatmin)/
    kostrov->ddlat;
  kostrov->nxny = kostrov->nx * kostrov->ny;
  if(kostrov->nx < 1 || kostrov->ny < 1){
    fprintf(stderr,"setup_bins: error: nx: %i ny: %i, lon: %g - %g - %g, lat: %g - %g - %g\n",
	    kostrov->nx,kostrov->ny,
	    kostrov->dlonmin,kostrov->ddlon,kostrov->dlonmax,
	    kostrov->dlatmin,kostrov->ddlat,kostrov->dlatmax);
    exit(-1);
  }
  fprintf(stderr,"setup_kostrov: using magnitudes from %g to %g, depths from %g to %g\n",
	  kostrov->minmag,kostrov->maxmag,
	  kostrov->mindepth,kostrov->maxdepth);
  
  fprintf(stderr,"setup_bins: setting up %i bins for -R%g/%g/%g/%g -I%g/%g nx: %i ny %i\n",
	  kostrov->nxny,
	  kostrov->dlonmin,
	  kostrov->dlonmax,
	  kostrov->dlatmin,
	  kostrov->dlatmax,
	  kostrov->ddlon,kostrov->ddlat,kostrov->nx,kostrov->ny);
  
  
  kostrov->bin = (struct bn *)
    realloc(kostrov->bin,kostrov->nxny * sizeof(struct bn));
  if(!kostrov->bin)
    BC_MEMERROR("setup_bins");
  for(i=0;i < kostrov->nxny;i++){
    kostrov->bin[i].quake = (unsigned int *)malloc(sizeof(unsigned int));
    kostrov->bin[i].weight = (BC_CPREC *)malloc(sizeof(BC_CPREC));
  }
  /* area without latitude correction (done below) */
  darea = BC_D2R(kostrov->ddlon) * BC_D2R(kostrov->ddlat) * BC_RADIUS * BC_RADIUS;

  /* in center of bin */
  xmin = kostrov->dlonmin + kostrov->ddlon/2.;
  ymin = kostrov->dlatmin + kostrov->ddlat/2.;
  
  /* set all to zero */
  clear_bins(catalog);		/*  */
  for(i=0;i < kostrov->nx;i++){
    for(j=0;j < kostrov->ny;j++){
      /* 
	 setup geography 
      */
      ind = i * kostrov->ny + j;
      /* center coordinates */
      kostrov->bin[ind].dlon = xmin + kostrov->ddlon * (BC_CPREC)i;
      kostrov->bin[ind].dlat = ymin + kostrov->ddlat * (BC_CPREC)j;
      /* in radians */
      kostrov->bin[ind].lon =  BC_D2R(kostrov->bin[ind].dlon);
      kostrov->bin[ind].lat =  BC_D2R(kostrov->bin[ind].dlat);
      kostrov->bin[ind].coslat = cos((double)kostrov->bin[ind].lat);
      /* 
	 spherical approximation for area, in km^2
      */
      kostrov->bin[ind].area =   kostrov->bin[ind].coslat  * darea;
    }
  }
  kostrov->init = BC_TRUE;
}

void clear_bins(struct cat *catalog)
{
  int i,k;
  struct kostrov_sum *kostrov; 
  kostrov = catalog->sum;
  
  kostrov->mtot = 0.0;		/* total moment of summation */
  for(i=0;i < kostrov->nxny;i++){
    kostrov->bin[i].sumw = 0.;	 /* sum of weights, number of
					    entries if no weights */
    kostrov->bin[i].n = 0;
    kostrov->bin[i].quake = (unsigned int *)
      realloc(kostrov->bin[i].quake,sizeof(unsigned int));

    kostrov->bin[i].weight = (BC_CPREC *)
      realloc(kostrov->bin[i].weight,sizeof(BC_CPREC));
    
    kostrov->bin[i].me = kostrov->bin[i].men = kostrov->bin[i].mens = 0.0;
    for(k=0;k < 6;k++){
      kostrov->bin[i].m[k] = kostrov->bin[i].mn[k]  = kostrov->bin[i].mnloc[k] =
	kostrov->bin[i].smn[k] = 0.0;
      kostrov->bin[i].s[k] = kostrov->bin[i].ds[k] =NAN;
    }
  }
}

/* 
   
   print filled cells, in spherical system rr rt rp tt tp pp

*/
void print_kostrov_bins(struct cat *catalog, char *filename,BC_BOOLEAN monte_carlo)
{
  FILE *out1,*out2;
  int i,j,ind,k;
  BC_CPREC mn,mn2,m;
  char outname1[300],outname2[300];
  BC_CPREC mscale;
  struct kostrov_sum *kostrov;
  if(!catalog->sum->init){
    fprintf(stderr,"sum_kostrov_bins: sum bounds not initialized\n");
    exit(-1);
  }
  kostrov = catalog->sum;

  /* 

     normalized

  */
  m=mn=mn2=0.;
  sprintf(outname1,"%s.norm.dat",filename);
  out1 = myopen(outname1,"w","print_kostrov_bins");
  for(i=0;i < kostrov->nx;i++)
    for(j=0;j < kostrov->ny;j++){
      ind =  i * kostrov->ny + j;
      if(kostrov->bin[ind].n){
	mn += (BC_CPREC)kostrov->bin[ind].n;
	mn2 += (BC_CPREC)(kostrov->bin[ind].n * kostrov->bin[ind].n);
	m += 1.0;			/* number filled */
	/* normalize */
	normalize_tens6(kostrov->bin[ind].mn);
	for(k=0;k < 6;k++){
	  //fprintf(stderr,"%g\n",tensor6_norm(kostrov->bin[ind].mn));
	  fprintf(out1,"%20.8e ",kostrov->bin[ind].mn[k]); /* norm,
							      don't
							      take out
							      area */
	}
	/* print bin centers */
	fprintf(out1,"%8.3f %8.3f %12i %12.5e %12.5e\n",
		kostrov_bdlon(ind,kostrov),kostrov_bdlat(ind,kostrov),
		kostrov->bin[ind].n,kostrov->bin[ind].men,kostrov->bin[ind].mens);
      }
    }
  fclose(out1);
  fprintf(stderr,"print_kostrov_bins: written normalized summations to %s,   %i out of %i cells, on avg %g +/- %g events\n",
	  outname1,(int)m,
	  kostrov->nxny,mn/m,sqrt ((m * mn2 - mn * mn) / (m*(m-1.))));
  if(monte_carlo){
    /* 

       normalized sigma tensor

    */
    m=mn=0.;
    sprintf(outname1,"%s.norm.sigma.dat",filename);
    out1 = myopen(outname1,"w","print_kostrov_bins");
    for(i=0;i < kostrov->nx;i++)
      for(j=0;j < kostrov->ny;j++){
	ind =  i * kostrov->ny + j;
	if(kostrov->bin[ind].n){
	  mn += (BC_CPREC)kostrov->bin[ind].n;
	  m += 1.0;
	  for(k=0;k < 6;k++)
	    fprintf(out1,"%20.8e ",kostrov->bin[ind].smn[k]); /* norm, don't take out area */
	  fprintf(out1,"%8.3f %8.3f %12i\n",
		  kostrov_bdlon(ind,kostrov),kostrov_bdlat(ind,kostrov),
		  kostrov->bin[ind].n);
	}
      }
    fclose(out1);
    fprintf(stderr,"print_kostrov_bins: written sigma of normalized summations to %s,   %i out %i cells, on avg %g entries\n",
	    outname1,(int)m,kostrov->nxny,mn/m);
  }
  /* 


     scaled tensor
     
  */
  mscale = kostrov->mtot/m;

  sprintf(outname2,"%s.scaled.dat",filename);
  out2 = myopen(outname2,"w","print_kostrov_bins");
  for(i=0;i < kostrov->nx;i++)
    for(j=0;j < kostrov->ny;j++){
      ind =  i * kostrov->ny + j;
      if(kostrov->bin[ind].n){
	for(k=0;k<6;k++)
	  fprintf(out2,"%20.8e ",
		  kostrov->bin[ind].m[k]/kostrov->bin[ind].area/mscale*1e6);	/* scaled p */
	fprintf(out2,"%8.3f %8.3f %6i %12.5e\n",
		kostrov_bdlon(ind,kostrov),kostrov_bdlat(ind,kostrov),
		kostrov->bin[ind].n,kostrov->bin[ind].me/kostrov->bin[ind].area/mscale*1e6);
      }
    }
  fclose(out2);
  fprintf(stderr,"print_kostrov_bins: written scaled     summations to %s, using mean moment %.4e for scale (%.3e/%i)\n",
	  outname2,mscale,kostrov->mtot,(int)m);

}

/* 
   output of stress tensors for each bin

*/
void print_stress_tensors(struct cat *catalog, char *filename)
{
  FILE *out1;
  int i,k,m;
  char outname1[300];
  struct kostrov_sum *kostrov;
  BC_CPREC t1[6],t2[6],dt[6];
  if(!catalog->sum->init){
    fprintf(stderr,"print_stress_tensor: sum not initialized\n");
    exit(-1);
  }
  kostrov = catalog->sum;

  sprintf(outname1,"%s.s.dat",filename);
  out1 = myopen(outname1,"w","print_kostrov_bins");
  for(m=i=0;i < kostrov->nxny;i++)
    if(kostrov->bin[i].n > BC_NQUAKE_LIM_FOR_STRESS){
      m++;
      for(k=0;k < 6;k++)	/* stress tensor */
	fprintf(out1,"%8.4f ",kostrov->bin[i].s[k]);
      fprintf(out1,"\t%8.3f %8.3f %12i",
	      kostrov_bdlon(i,kostrov),kostrov_bdlat(i,kostrov),kostrov->bin[i].n);
      for(k=0;k < 6;k++)	/* stress tensor uncertainty */
	fprintf(out1,"%8.4f ",kostrov->bin[i].ds[k]);
      /* sigma norm */
      fprintf(out1,"\t%8.4f\n",tensor6_norm(kostrov->bin[i].ds));
    }
  fclose(out1);
  fprintf(stderr,"print_stress_tensors: written stress tensors from Michael inversion to %s, %i out of %i cells filled\n",
	  outname1,m,kostrov->nxny);

  /* difference with normalized strain */
  sprintf(outname1,"%s.smn.dat",filename);
  out1 = myopen(outname1,"w","print_kostrov_bins");
  for(m=i=0;i < kostrov->nxny;i++)
    if(kostrov->bin[i].n > BC_NQUAKE_LIM_FOR_STRESS){
      m++;
      for(k=0;k < 6;k++){	
	t1[k] = kostrov->bin[i].s[k];/* stress tensor */
	t2[k] = kostrov->bin[i].mn[k]; /* normalized strain */
      }
      normalize_tens6(t1);normalize_tens6(t2);
      for(k=0;k < 6;k++){	/* difference tensor */
	dt[k] = t1[k] - t2[k];
      }
      remove_trace(dt);
      for(k=0;k < 6;k++){	/* difference tensor */
	fprintf(out1,"%8.4f ",dt[k]);
      }
      fprintf(out1,"\t%8.3f %8.3f %12i",
	      kostrov_bdlon(i,kostrov),kostrov_bdlat(i,kostrov),kostrov->bin[i].n);
      for(k=0;k < 6;k++)	/* stress tensor uncertainty */
	fprintf(out1,"%8.4f ",kostrov->bin[i].ds[k]);
      /* sigma norm */
      fprintf(out1,"\t%8.4f\n",tensor6_norm(dt));
    }
  fclose(out1);
  fprintf(stderr,"print_stress_tensors: written stress tensors minus normalized to %s\n",outname1);


  
}

/* print coordinates in degree */
BC_CPREC kostrov_bdlon(int i, struct kostrov_sum *kostrov)
{
  return kostrov->bin[i].dlon - kostrov->ddlon/2.;
}
BC_CPREC kostrov_bdlat(int i, struct kostrov_sum *kostrov)
{
  return kostrov->bin[i].dlat - kostrov->ddlat/2.;
}

/* 
   

   in units of 10^15

*/
void print_summed_moment(struct cat *catalog, char *filename)
{
  FILE *out1;
  int i,j,ind;
  char outname1[300];
  BC_CPREC me_min=1e30,me_max=-1e30;
  struct kostrov_sum *kostrov;
  if(!catalog->sum->init){
    fprintf(stderr,"print_summed_moment: sum bounds not initialized\n");
    exit(-1);
  }
  kostrov = catalog->sum;
  sprintf(outname1,"%s.dat",filename);
  out1 = myopen(outname1,"w","print_summed_moment");
  for(i=0;i < kostrov->nx;i++)
    for(j=0;j < kostrov->ny;j++){
      ind =  i * kostrov->ny + j;
      fprintf(out1,"%g %g %g\n",kostrov_bdlon(ind,kostrov),kostrov_bdlat(ind,kostrov),
	      kostrov->bin[ind].me);
      if(kostrov->bin[ind].me < me_min)
	me_min = kostrov->bin[ind].me;
      if(kostrov->bin[ind].me > me_max)
	me_max = kostrov->bin[ind].me;
    }
  fclose(out1);
  fprintf(stderr,"print_summed_moment: printed %i bins to %s with me between %g/%g \n",
	  kostrov->nxny,outname1,me_min,me_max);
}



/* 

catalog handling routines

eps is in [km]

first catalog wins
*/

void merge_catalog(struct cat *a,struct cat *b, struct cat *c,BC_CPREC eps,BC_CPREC vert_eps)
{
  int i,j,d;
  BC_CPREC dist;


  d = 0;
  for(i=0;i<a->n;i++){
    if(i % 500 == 0)
      fprintf(stderr,"merge_catalog: checking doubles: %7i out of %7i tests with %8i\r",i,a->n,b->n);
    for(j=0;j < b->n;j++){
      dist = distance(a,b,i,j);
      if((dist < eps) && 	/* horizontal distance */
	 (fabs(a->quake[i].depth - b->quake[j].depth) < vert_eps) && /* vertical  distance */
	 (fabs(a->quake[i].tsec - b->quake[j].tsec) < 0.000125) /* temporal <~ two minutes */
	 ){
	/* found BC_CPREC */
	d++;
	if(a->quake[i].mag >  b->quake[j].mag)
	  b->quake[j].deleted = BC_TRUE; /* large event wins */
	else
	  a->quake[i].deleted = BC_TRUE; /* large event wins */
      }
    }
  }
  fprintf(stderr,"\nmerge_catalog: deleted %i quakes out of %i + %i = %i, merging catalog\n",
	  d,a->n,b->n,a->n+b->n);
  for(i=0;i < a->n;i++)
    if(!a->quake[i].deleted){
      copy_quake((a->quake+i),(c->quake+c->n));
      c->n += 1;
      make_room_for_quake(c);
    }
  for(i=0;i < b->n;i++)
    if(!b->quake[i].deleted){
      copy_quake((b->quake+i),(c->quake+c->n));
      c->n += 1;
      make_room_for_quake(c);
    }
  fprintf(stderr,"merge_catalog: merged catalog has %i events\n",c->n);
}


int print_catalog(char *filename, struct cat *catalog, int mode)
{
  FILE *out;
  int i;
  /* open file */
  out = myopen(filename,"w","print_catalog");
  /* read quakes */
  for(i=0;i < catalog->n;i++)
    print_quake(out,catalog->quake[i],mode);
  fprintf(stderr,"print_catalog: printed %8i events to %30s, mode: %s\n",
	  catalog->n,filename,mode_name(mode));
  fclose(out);
  return 0;
}


int read_catalog(char *filename, struct cat *catalog, int mode)
{
  FILE *in;
  int i,ndup,hit,ilim;
  long int init_random_seed = -1; /* change to create new numbers */
  static BC_BOOLEAN  check_duplicates = BC_TRUE;
  BC_CPREC dist;
  struct tm ts;
  char tbuf1[80],tbuf2[80];
  time_t timet;
  /* 
     init catalog 
  */
  create_catalog(catalog,init_random_seed);
  /* open file */
  in = myopen(filename,"r","handle_catalog");
  ndup = 0;
  /* read quakes */
  while(read_quake(in,(catalog->quake + catalog->n),mode)){
    hit = 0;
    if(check_duplicates){
      /* look for duplicates within last few */
      ilim = catalog->n - 30;
      if(ilim < 0)ilim = 0;
      for(i=ilim;(i<catalog->n-1)&&(!hit);i++){
	dist = distance(catalog,catalog,i,catalog->n);
	if((dist < .200) && 	/* horizontal distance, 500 m */
	   (fabs(catalog->quake[i].depth - catalog->quake[catalog->n].depth) <  BC_DEP_CLOSE) && /* vertical  distance, km */
	   (fabs(catalog->quake[i].tsec - catalog->quake[catalog->n].tsec) <  BC_TIME_CLOSE)){ /* 0.00025 is <~1 min */
	  /* 
	     we found a duplicate 
	  */
	  ndup++;
	  if(fabs(catalog->quake[i].mag - catalog->quake[catalog->n].mag) > 2){
	    /* notify for large mag difference */
	    fprintf(stderr,"read_catalog: from %s: duplicate: %i matches %i, hdist: %g vdist: %g tdist: %g\n",
		    filename,catalog->n+1,i+1,dist,
		    fabs(catalog->quake[i].depth - catalog->quake[catalog->n].depth),
		    fabs(catalog->quake[i].tsec - catalog->quake[catalog->n].tsec)/BC_TIME_CLOSE);
	    fprintf(stderr,"%5i\t",i+1);print_quake(stderr,catalog->quake[i],BC_AKI);
	    fprintf(stderr,"%5i\t",catalog->n+1);print_quake(stderr,catalog->quake[catalog->n],BC_AKI);
	  }
	  if(catalog->quake[catalog->n].mag > catalog->quake[i].mag){
	    /* new magnitude is larger than old, use new one */
	    copy_quake((catalog->quake+catalog->n),(catalog->quake+i));
	  }
	  hit=1;
	}
      }
    }
    if(!hit){
      /* add a new quake to the catalog */
      catalog->n += 1;
      make_room_for_quake(catalog);
    }
  }
  catalog->tmin = catalog->tmax = catalog->quake[0].tsec;
  for(i=0;i < catalog->n;i++){
    //fprintf(stderr,"%g %g %g %g\n",catalog->quake[i].lon,catalog->quake[i].lat,catalog->quake[i].depth, catalog->quake[i].mag);
    catalog->quake[i].lkm = quake_scale_lkm(catalog->quake[i].m0); 
     
    if(catalog->quake[i].tsec > catalog->tmax)
      catalog->tmax = catalog->quake[i].tsec;
    if(catalog->quake[i].tsec < catalog->tmin)
      catalog->tmin = catalog->quake[i].tsec;

    if(catalog->quake[i].mag < catalog->minmag)
      catalog->minmag = catalog->quake[i].mag; 
    if(catalog->quake[i].mag > catalog->maxmag)
      catalog->maxmag = catalog->quake[i].mag; 

    if(catalog->quake[i].depth < catalog->mindepth)
      catalog->mindepth = catalog->quake[i].depth;
    if(catalog->quake[i].depth > catalog->maxdepth)
      catalog->maxdepth = catalog->quake[i].depth;

    if(catalog->quake[i].dlon < catalog->minlond)
      catalog->minlond = catalog->quake[i].dlon;
    if(catalog->quake[i].dlon > catalog->maxlond)
      catalog->maxlond = catalog->quake[i].dlon;

    if(catalog->quake[i].dlat < catalog->minlatd)
      catalog->minlatd = catalog->quake[i].dlat;
    if(catalog->quake[i].dlat > catalog->maxlatd)
      catalog->maxlatd = catalog->quake[i].dlat;

    if(catalog->quake[i].lkm > catalog->lkm_max)
      catalog->lkm_max = catalog->quake[i].lkm;
    if(catalog->quake[i].lkm < catalog->lkm_min)
      catalog->lkm_min = catalog->quake[i].lkm;
    
  }
  catalog->tcenter = (catalog->tmin+catalog->tmax)/2;
  catalog->trange = (catalog->tmax-catalog->tmin);
  fprintf(stderr,"read_catalog: read %8i events from %s, mode: %s\n",
	  catalog->n,filename,mode_name(mode));
  fprintf(stderr,"read_catalog: input times are between %15is and %15is, with center at %15is\n",
	  (int)catalog->tmin,(int)catalog->tmax,(int)(catalog->tcenter+.5));

  // Format time, "ddd yyyy-mm-dd hh:mm:ss zzz"
  timet=(time_t)catalog->tmin;ts = *localtime(&timet);
  strftime(tbuf1, sizeof(tbuf1), "%a %Y-%m-%d %H:%M:%S %Z", &ts);
  timet=(time_t)catalog->tmax;ts = *localtime(&timet);
  strftime(tbuf2, sizeof(tbuf2), "%a %Y-%m-%d %H:%M:%S %Z", &ts);
  fprintf(stderr,"read_catalog: if in UNIX times, from %s to %s\n",tbuf1,tbuf2);

	  
  fprintf(stderr,"read_catalog: range is -R%g/%g/%g/%g, depths within %g/%g, magnitude within %g/%g, inferred L: %g/%g km (%s)\n",
	  catalog->minlond,catalog->maxlond,catalog->minlatd,catalog->maxlatd,
	  catalog->mindepth,catalog->maxdepth,
	  catalog->minmag,catalog->maxmag,
	  catalog->lkm_min,catalog->lkm_max,
	  (catalog->is_xy)?("Cartesian"):("geographic"));

  if(check_duplicates)
    fprintf(stderr,"read_catalog: found %i duplicates\n",ndup);
  fclose(in);
  return 0;
}

void copy_quake(struct qke *a, struct qke *b)
{
  memcpy(b,a,sizeof(struct qke));
}

char *mode_name(int mode)
{
  switch(mode){
  case BC_AKI:
    return("AKI");
    break;
  case BC_CMT:
    return("CMT");
    break;
  case BC_ENG:
    return("ENG");
    break;
  case BC_CMT_FP:
    return("CMT_FP");
    break;
  default:
    fprintf(stderr,"mode_name: mode %i undefined\n",mode);
    exit(-1);
  }
}
int read_quake(FILE *in, struct qke *quake, int mode)
{
  switch(mode){
  case BC_AKI:
    return(read_quake_aki(in,quake));
    break;
  case BC_CMT:
    return(read_quake_cmt(in,quake));
    break;
  case BC_ENG:
    return(read_quake_eng(in,quake));
    break;
  default:
    fprintf(stderr,"read_quake: mode %i undefined\n",mode);
    exit(-1);
  }
}
void print_quake(FILE *out, struct qke quake, int mode)
{
 switch(mode){
 case BC_AKI:
    print_quake_aki(out, quake); /* Aki style */
    break;
 case BC_CMT:
   print_quake_cmt(out,quake);	/* gCMT moment tensors */
   break;
 case BC_CMT_FP:
   print_quake_cmt_fp(out,quake); /* best fit fault planes */
   break;
 default:
    fprintf(stderr,"print_quake: mode %i undefined\n",mode);
    exit(-1);
  }
}

/* 
   read AKI & Richards format, assume all angles and coordinates are in degree 
   and last column is time
*/
//
// print(lon,lat,dep,str[1],dip[1],rake[1],mag,lon,lat,tsec);
//
int read_quake_aki(FILE *in, struct qke *quake)
{
  BC_CPREC tmp1,tmp2;
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    &(quake->dlon),&(quake->dlat),&(quake->depth),
	    &(quake->strike),&(quake->dip),&(quake->rake),
	    &(quake->mag),&tmp1,&tmp2,&quake->tsec)==10){
    /* convert */
    quake->strike = BC_D2R(quake->strike); /* strike/dip/rake in radians */
    quake->dip    = BC_D2R(quake->dip);
    quake->rake   = BC_D2R(quake->rake);

    quake->plon=quake->dlon;
    quake->plat=quake->dlat;
  
    quake->lon = BC_D2R(quake->dlon); /* radian versions */
    quake->lat = BC_D2R(quake->dlat);
    /* get the other plane */
    find_alt_plane(quake->strike, quake->dip, quake->rake,
		   &(quake->strike2), &(quake->dip2), &(quake->rake2));
    /* those won't make sense for x-y, but won't hurt */
    quake->coslat = cos(quake->lat);
    /*  */
    /* compute potency */
    aki2mom(quake->strike,quake->dip,quake->rake,quake->mag,
	    quake->m,&(quake->m0));
    return 1;
  }else{
    return 0;
  }
}
/* 

   modified mdat format where eevnt id is replaced by time_sec_since_epoch

   X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp[dyn cm], newX, newY, tsec_unix

   also compute best fit nodal planes

*/
int read_quake_cmt(FILE *in, struct qke *quake)
{
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i %lf %lf %lf\n",
	    &(quake->dlon),&(quake->dlat),&(quake->depth),
	    (quake->m+BC_RR),(quake->m+BC_TT),(quake->m+BC_PP),
	    (quake->m+BC_RT),(quake->m+BC_RP),(quake->m+BC_TP),
	    &(quake->exp),&(quake->plon),&(quake->plat),
	    &quake->tsec) == 13){
    quake->lon = BC_D2R(quake->dlon); /* radian versions */
    quake->lat = BC_D2R(quake->dlat);
    quake->coslat = cos(quake->lat);  

    quake->exp -= 7;		/* convert to Nm */

    quake->m0 = pow(10,(BC_CPREC)quake->exp) * scalar_mom(quake->m);
    quake->mag = mom2mag(quake->m0);
    /* 
       compute fault planes (this I checked against psmeca)
    */
    tensor2fpangle(quake->m, &(quake->strike), &(quake->dip), &(quake->rake),
		   &(quake->strike2), &(quake->dip2), &(quake->rake2));
    
    //fprintf(stderr,"%i %g %g %g\n",quake->exp,scalar_mom(quake->m),quake->m0,quake->mag);
    return 1;
  }else{
    return 0;
  }
}
/* FORMAT: lon lat depth mw tsec_unix */
int read_quake_eng(FILE *in, struct qke *quake)
{
  if(fscanf(in,"%lf %lf %lf %lf %lf\n",
	    &(quake->dlon),&(quake->dlat),&(quake->depth),
	    &(quake->mag),&quake->tsec) == 5){
    quake->plon=quake->dlon;
    quake->plat=quake->dlat;
    
    quake->lon = BC_D2R(quake->dlon); /* radian versions */
    quake->lat = BC_D2R(quake->dlat);
   
    quake->coslat = cos(quake->lat);
    quake->m0 = mag2mom(quake->mag);
    return 1;
  }else{
    return 0;
  }
}





/*  */
void print_quake_aki(FILE *out, struct qke quake)
{
  fprintf(out,"%9.4f %9.4f %6.2f %8.2f  %8.2f %8.2f %5.2f %8.2f %8.2f %.8e\n",
	  quake.dlon,quake.dlat,quake.depth,
	  BC_R2D(quake.strike),BC_R2D(quake.dip),BC_R2D(quake.rake),
	  quake.mag,quake.plon,quake.plat,quake.tsec);
}
//
// X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp[dyn cm], newX, newY, time
//
void print_quake_cmt(FILE *out, struct qke quake)
{
  fprintf(out,"%9.4f %9.4f %6.2f\t%10.5f %10.5f %10.5f %10.5f %10.5f  %10.5f %3i %9.4f %9.4f %f\n",
	  quake.dlon,quake.dlat,quake.depth,
	  quake.m[BC_RR],quake.m[BC_TT],quake.m[BC_PP],
	  quake.m[BC_RT],quake.m[BC_RP],quake.m[BC_TP],
	  quake.exp + 7,quake.plon,quake.plat,quake.tsec);
}

/* 

   X , Y, depth, strike1, dip1, rake1, strike2, dip2, rake2, moment, newX, newY, event_title

   moment in units of dyn-cm
*/
void print_quake_cmt_fp(FILE *out, struct qke quake)
{
  int iexp;
  iexp = (int)(log10(quake.m0)+0.5);
  
  fprintf(out,"%9.4f %9.4f %6.2f\t%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\t%10.8f %i\t%9.4f %9.4f %f\n",
	  quake.dlon,quake.dlat,quake.depth,
	  BC_R2D(quake.strike),BC_R2D(quake.dip),BC_R2D(quake.rake),
	  BC_R2D(quake.strike2),BC_R2D(quake.dip2),BC_R2D(quake.rake2),
	  quake.m0/pow(10,iexp),iexp+7,quake.plon,quake.plat,quake.tsec);
}



/* make a fresh catalog */
void create_catalog(struct cat *catalog,long int init_seed)
{
  static BC_BOOLEAN init = BC_FALSE;
  catalog->n=0;
  catalog->quake = NULL;
  catalog->sum->bin = NULL;
  catalog->sum->init = BC_FALSE;
  make_room_for_quake(catalog);	/* make room for one quake */
  catalog->minmag = 100;
  catalog->maxmag = -100;
  catalog->mindepth = 1000;
  catalog->maxdepth = -1000;
  catalog->minlond = 1000;
  catalog->maxlond = -1000;
  catalog->minlatd = 1000;
  catalog->maxlatd = -1000;
  catalog->lkm_min = 6371;
  catalog->lkm_max = 0;
  if(!init){
    catalog->seed = init_seed;		/* default random number seed */
    BC_RGEN(&catalog->seed);
  }
  init = BC_TRUE;
}


/*  */
void make_room_for_quake(struct cat *catalog)
{
  catalog->quake = (struct qke *)
    realloc(catalog->quake,sizeof(struct qke)*(catalog->n+1));
  if(!catalog->quake)
    BC_MEMERROR("make_room_for_quake");
  catalog->quake[catalog->n].deleted = BC_FALSE;
}




/* 
   compute the weight of this quake given a distance kernel and
   possibly m0 (not used right now) for different modes, right now
   only one
*/
BC_CPREC quake_weight(BC_CPREC m0, BC_CPREC lkm, BC_CPREC dist,int mode)
{
  const double sqrt_two_pi = 2.506628274631;
  BC_CPREC sigma,weight=0,frac,frac2;
  switch(mode){
  case 0:
  case 1:
  case 2:
    /* 
       use estimated l_km as 6 sigma for a Gaussian
    */
    if(lkm < 1)
      lkm = 1;
    sigma = 5*lkm;

    frac = dist/sigma;
    if(frac < 10){
      frac2 = frac*frac;
      weight = 1/(sigma*sqrt_two_pi) * exp(-frac2/2);
    } /* else, zero */
    break;
  default:
    fprintf(stderr,"quake_weight: out of bounds with mode %i\n",mode);
    exit(-1);
  }

  return weight;
}

BC_CPREC quake_scale_lkm(BC_CPREC m0)
{
  BC_CPREC lkm,radius;
  /* M0 = 16/7 dsigma r^3 from Brune */
  const BC_CPREC dsigma = 2.9e6;	/* stress drop, 3 MPa */
  const BC_CPREC scale=7./16./dsigma;
  const BC_CPREC one_third = 1./3.;

  radius = pow(m0*scale,one_third)/1e3; /* radius in km */
  if(radius > 15){			/* account for strike slip setting */
    lkm = radius*radius/15;
  }else{
    lkm = radius;
  }
  return lkm;
}



/* 

   add quake to bin, increment counter

 */
void add_quake_to_bin_list(unsigned int iquake, struct bn *bin,BC_CPREC weight)
{
  bin->quake[bin->n] = iquake; /* add to list */
  bin->weight[bin->n] = weight;
  
  bin->n++;	       /* inc counter */
  bin->quake = (unsigned int *)realloc(bin->quake,(bin->n+1)*sizeof(unsigned int));/* make room */
  if(!bin->quake)
    BC_MEMERROR("add_quake_to_bin_list: quake");
  bin->weight = (BC_CPREC *)realloc(bin->weight,(bin->n+1)*sizeof(BC_CPREC));/* make room */
  if(!bin->weight)
    BC_MEMERROR("add_quake_to_bin_list: weight");
}

void assign_quake_angles(struct qke *quake,BC_CPREC *angles)
{
  angles[0] = quake->strike;
  angles[1] = quake->dip;
  angles[2] = quake->rake;
  
  angles[3] = quake->strike2;
  angles[4] = quake->dip2;
  angles[5] = quake->rake2;


}
