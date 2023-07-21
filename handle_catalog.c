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
      //fprintf(stderr,"%11g -> %11g, %11g -> %11g, %11g -> %11g\n",c->quake[c->n].lon, b->quake[found].lon,c->quake[c->n].lat, b->quake[found].lat,c->quake[c->n].depth, b->quake[found].depth);
      
      hor_distm  += hor_dist0;
      vert_distm += vert_dist0;
      temp_distm += temp_dist0;
      mag_distm  += mag_dist0;
      
      c->quake[c->n].lon = b->quake[found].lon;
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
  BC_CPREC eps_angle = BC_DEG2RAD(30.0); /* sigma for random */
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
	dx = catalog->quake[i].lon - kostrov->lonmin;
	if(!catalog->is_xy){
	  if(dx > 360)
	    dx -= 360;
	}
	ix = (int)(dx/kostrov->dlon + 0.5);
	iy = (int)((catalog->quake[i].lat - kostrov->latmin)/kostrov->dlat + 0.5);
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

	if(catalog->is_xy){
	  dist = distance_cart(catalog->quake[i].lon,catalog->quake[i].lat,
			       kostrov->bin[j].lon,kostrov->bin[j].lat);
	}else{
	  dist = distance_geo(catalog->quake[i].lon,catalog->quake[i].lat,
			      kostrov->bin[j].lon,kostrov->bin[j].lat,
			      catalog->quake[i].coslat,
			      kostrov->bin[j].coslat);
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

BC_CPREC max_x_from_int_vector(BC_CPREC *x, int *y, int n)
{
  int i;
  int max = -1e8;
  BC_CPREC xmax;
  xmax = y[0];
  for(i=0;i<n;i++){
    if(y[i] > max){
      max = y[i];
      xmax = x[i];
    }
  }
  return xmax;
}

/* 

   call after initializing lon, lat range and sum->dlon/sum->dlat ! 

*/
void setup_kostrov(struct cat *catalog,
		   BC_CPREC lonmin,BC_CPREC lonmax,
		   BC_CPREC latmin,BC_CPREC latmax,
		   BC_CPREC mindepth,BC_CPREC maxdepth,
		   BC_CPREC dlon, BC_CPREC dlat,
		   BC_CPREC minmag,BC_CPREC maxmag,
		   int weighting_method)
{
  int i,j,ind;
  double xmin,ymin,darea;
  struct kostrov_sum *kostrov; 
  kostrov = catalog->sum;

  kostrov->mtot = 0.0;		/* total moment of summation */
  kostrov->weighting_method = weighting_method;
  
  kostrov->lonmin = lonmin;
  kostrov->lonmax = lonmax;
  kostrov->latmin = latmin;
  kostrov->latmax = latmax;

  kostrov->mindepth = mindepth;
  kostrov->maxdepth = maxdepth;

  kostrov->minmag = minmag;
  kostrov->maxmag = maxmag;
  
  kostrov->dlon = dlon;
  kostrov->dlat = dlat;
  
  kostrov->nx = (kostrov->lonmax - kostrov->lonmin)/
    kostrov->dlon;
  kostrov->ny = (kostrov->latmax - kostrov->latmin)/
    kostrov->dlat;
  kostrov->nxny = kostrov->nx * kostrov->ny;
  if(kostrov->nx < 1 || kostrov->ny < 1){
    fprintf(stderr,"setup_bins: error: nx: %i ny: %i, lon: %g - %g - %g, lat: %g - %g - %g\n",
	    kostrov->nx,kostrov->ny,
	    kostrov->lonmin,kostrov->dlon,kostrov->lonmax,
	    kostrov->latmin,kostrov->dlat,kostrov->latmax);
    exit(-1);
  }
  fprintf(stderr,"setup_kostrov: using magnitudes from %g to %g, depths from %g to %g\n",
	  kostrov->minmag,kostrov->maxmag,
	  kostrov->mindepth,kostrov->maxdepth);
  
  fprintf(stderr,"setup_bins: setting up %i bins for -R%g/%g/%g/%g -I%g/%g nx: %i ny %i\n",
	  kostrov->nxny,
	  kostrov->lonmin,
	  kostrov->lonmax,
	  kostrov->latmin,
	  kostrov->latmax,
	  kostrov->dlon,kostrov->dlat,kostrov->nx,kostrov->ny);
  
  
  kostrov->bin = (struct bn *)
    realloc(kostrov->bin,kostrov->nxny * sizeof(struct bn));
  if(!kostrov->bin)
    BC_MEMERROR("setup_bins");
  for(i=0;i < kostrov->nxny;i++){
    kostrov->bin[i].quake = (unsigned int *)malloc(sizeof(unsigned int));
    kostrov->bin[i].weight = (BC_CPREC *)malloc(sizeof(BC_CPREC));
  }
  /* area without latitude correction (done below) */
  darea = BC_DEG2RAD(kostrov->dlon) * BC_DEG2RAD(kostrov->dlat) * BC_RADIUS * BC_RADIUS;

  /* in center of bin */
  xmin = kostrov->lonmin + kostrov->dlon/2.;
  ymin = kostrov->latmin + kostrov->dlat/2.;
  
  /* set all to zero */
  clear_bins(catalog);		/*  */
  for(i=0;i < kostrov->nx;i++){
    for(j=0;j < kostrov->ny;j++){
      /* 
	 setup geography 
      */
      ind = i * kostrov->ny + j;
      /* center coordinates */
      kostrov->bin[ind].lon = xmin + kostrov->dlon * (double)i;
      kostrov->bin[ind].lat = ymin + kostrov->dlat * (double)j;

      kostrov->bin[ind].coslat = (float)cos((double)kostrov->bin[ind].lat/BC_PIF);

      /* 
	 spherical approximation for area, in km^2
      */
      kostrov->bin[ind].area = (float)
	cos(BC_DEG2RAD((double)kostrov->bin[ind].lat)) * darea;
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
  int i,j,ind,k,m,mn,mn2;
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
  m=mn=mn2=0;
  sprintf(outname1,"%s.norm.dat",filename);
  out1 = myopen(outname1,"w","print_kostrov_bins");
  for(i=0;i < kostrov->nx;i++)
    for(j=0;j < kostrov->ny;j++){
      ind =  i * kostrov->ny + j;
      if(kostrov->bin[ind].n){
	mn += kostrov->bin[ind].n;
	mn2 += kostrov->bin[ind].n * kostrov->bin[ind].n;
	m++;
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
		kostrov_blon(ind,kostrov),kostrov_blat(ind,kostrov),
		kostrov->bin[ind].n,kostrov->bin[ind].men,kostrov->bin[ind].mens);
      }
    }
  fclose(out1);
  fprintf(stderr,"print_kostrov_bins: written normalized summations to %s,   %i out of %i cells, on avg %g +/- %g entries\n",
	  outname1,m,kostrov->nxny,(double)mn/(double)m,
	  sqrt((((double)m) * ((double)mn2) - (((double)mn) * ((double)mn))) /
	       ((double)m*((double)m-1.))));
  if(monte_carlo){
    /* 

       normalized sigma tensor

    */
    m=mn=0;
    sprintf(outname1,"%s.norm.sigma.dat",filename);
    out1 = myopen(outname1,"w","print_kostrov_bins");
    for(i=0;i < kostrov->nx;i++)
      for(j=0;j < kostrov->ny;j++){
	ind =  i * kostrov->ny + j;
	if(kostrov->bin[ind].n){
	  mn += kostrov->bin[ind].n;
	  m++;
	  for(k=0;k < 6;k++)
	    fprintf(out1,"%20.8e ",kostrov->bin[ind].smn[k]); /* norm, don't take out area */
	  fprintf(out1,"%8.3f %8.3f %12i\n",
		  kostrov_blon(ind,kostrov),kostrov_blat(ind,kostrov),
		  kostrov->bin[ind].n);
	}
      }
    fclose(out1);
    fprintf(stderr,"print_kostrov_bins: written sigma of normalized summations to %s,   %i out %i cells, on avg %g entries\n",
	    outname1,m,kostrov->nxny,(double)mn/(double)m);
  }
  /* 


     scaled tensor
     
  */
  mscale = kostrov->mtot/(BC_CPREC)m;

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
		kostrov_blon(ind,kostrov),kostrov_blat(ind,kostrov),
		kostrov->bin[ind].n,kostrov->bin[ind].me/kostrov->bin[ind].area/mscale*1e6);
      }
    }
  fclose(out2);
  fprintf(stderr,"print_kostrov_bins: written scaled     summations to %s, using mean moment %.4e for scale (%.3e/%i)\n",
	  outname2,mscale,kostrov->mtot,m);

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
      fprintf(out1,"\t%8.3f %8.3f %12i",kostrov_blon(i,kostrov),kostrov_blat(i,kostrov),kostrov->bin[i].n);
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
      fprintf(out1,"\t%8.3f %8.3f %12i",kostrov_blon(i,kostrov),kostrov_blat(i,kostrov),kostrov->bin[i].n);
      for(k=0;k < 6;k++)	/* stress tensor uncertainty */
	fprintf(out1,"%8.4f ",kostrov->bin[i].ds[k]);
      /* sigma norm */
      fprintf(out1,"\t%8.4f\n",tensor6_norm(dt));
    }
  fclose(out1);
  fprintf(stderr,"print_stress_tensors: written stress tensors minus normalized to %s\n",outname1);


  
}


BC_CPREC kostrov_blon(int i, struct kostrov_sum *kostrov)
{
  return kostrov->bin[i].lon - kostrov->dlon/2;
}
BC_CPREC kostrov_blat(int i, struct kostrov_sum *kostrov)
{
  return kostrov->bin[i].lat - kostrov->dlat/2;
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
      fprintf(out1,"%g %g %g\n",kostrov_blon(ind,kostrov),kostrov_blat(ind,kostrov),
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
 
convert event in Aki convention (angles in deg) into moment, magnitude
is Ml

*/

void aki2mom(BC_CPREC strike, BC_CPREC dip, BC_CPREC rake, BC_CPREC mag, 
	     BC_CPREC *m, BC_CPREC *m0)
{
  BC_CPREC phi,delta,gamma,phi2,delta2,sin_2delta,cos_gamma,cos_2delta;
  BC_CPREC sin_phi,cos2_phi,sin2_phi,sin_delta,sin_gamma,cos_delta,
    cos_phi,cos_2phi,sin_2phi;
  // angles
  phi =   BC_DEG2RAD(strike);		        // strike
  delta = BC_DEG2RAD(dip);			// dip
  gamma = BC_DEG2RAD(rake);			// rake

  /* use moment */
  *m0 = mag2mom(mag);

  //
  // trigonometry
  phi2 = phi * 2.0;
  delta2 = delta * 2.0;

  sin_phi   = sin(phi);  
  cos_phi   = cos(phi);  
  sin_2phi  = sin(phi2);  
  cos_2phi  = cos(phi2);
  sin2_phi = sin_phi*sin_phi;
  cos2_phi = cos_phi*cos_phi;
  sin_delta = sin(delta);
  cos_delta = cos(delta);
  sin_2delta= sin(delta2);
  cos_2delta= cos(delta2);
  sin_gamma = sin(gamma);
  cos_gamma = cos(gamma);
  //
  // moment components, x, y, z, = East, North, Up = P, -T, R
  //
  m[BC_RR] =  (sin_2delta * sin_gamma); /*                                                     UU */
  m[BC_RT] = -(cos_delta  * cos_gamma * cos_phi  +     cos_2delta * sin_gamma * sin_phi); /*  -NU */
  m[BC_RP] =  (cos_delta  * cos_gamma * sin_phi  -     cos_2delta * sin_gamma * cos_phi); /*   EU */
  m[BC_TT] = -(sin_delta  * cos_gamma * sin_2phi +     sin_2delta * sin_gamma * sin2_phi); /*  NN */
  m[BC_TP] = -(sin_delta  * cos_gamma * cos_2phi + 0.5*sin_2delta * sin_gamma * sin_2phi); /* -EN */
  m[BC_PP] =  (sin_delta  * cos_gamma * sin_2phi -     sin_2delta * sin_gamma * cos2_phi); /*  EE */
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
	   (fabs(catalog->quake[i].depth - catalog->quake[catalog->n].depth) < .400) && /* vertical  distance, km */
	   (fabs(catalog->quake[i].tsec - catalog->quake[catalog->n].tsec) < 0.0003)){ /* 0.00025 is <~1 min */
	  /* 
	     we found a duplicate 
	  */
	  ndup++;
	  if(fabs(catalog->quake[i].mag - catalog->quake[catalog->n].mag) > 2){
	    /* notify for large mag difference */
	    fprintf(stderr,"read_catalog: from %s: duplicate: %i matches %i, hdist: %g vdist: %g tdist: %g\n",
		    filename,catalog->n+1,i+1,dist,
		    fabs(catalog->quake[i].depth - catalog->quake[catalog->n].depth),
		    fabs((catalog->quake[i].tsec - catalog->quake[catalog->n].tsec)/0.00025));
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
  for(i=1;i < catalog->n;i++){
    
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
    if(catalog->quake[i].lon < catalog->minlon)
      catalog->minlon = catalog->quake[i].lon;
    if(catalog->quake[i].lon > catalog->maxlon)
      catalog->maxlon = catalog->quake[i].lon;
    if(catalog->quake[i].lat < catalog->minlat)
      catalog->minlat = catalog->quake[i].lat;
    if(catalog->quake[i].lat > catalog->maxlat)
      catalog->maxlat = catalog->quake[i].lat;
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
	  catalog->minlon,catalog->maxlon,catalog->minlat,catalog->maxlat,
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
  BC_CPREC rlat,tmp1,tmp2;
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	    &(quake->lon),&(quake->lat),&(quake->depth),
	    &(quake->strike),&(quake->dip),&(quake->rake),
	    &(quake->mag),&tmp1,&tmp2,&quake->tsec)==10){
    if(quake->lon <0)
      quake->lon += 360;
    /* get the other plane */
    find_alt_plane(quake->strike, quake->dip, quake->rake,
		   &(quake->strike2), &(quake->dip2), &(quake->rake2));
    /* those won't make sense for x-y, but won't hurt */
    rlat = BC_DEG2RAD(quake->lat);
    quake->coslat = cos((double)rlat);
    /*  */
    quake->plon = quake->lon;	/* plotting location */
    quake->plat = quake->lat;
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
  BC_CPREC rlat;
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i %lf %lf %lf\n",
	    &(quake->lon),&(quake->lat),&(quake->depth),
	    (quake->m+BC_RR),(quake->m+BC_TT),(quake->m+BC_PP),
	    (quake->m+BC_RT),(quake->m+BC_RP),(quake->m+BC_TP),
	    &(quake->exp),&(quake->plon),&(quake->plat),
	    &quake->tsec) == 13){
    if(quake->lon <0)
      quake->lon += 360;
    quake->exp -= 7;		/* convert to Nm */
    rlat = BC_DEG2RAD(quake->lat);
    quake->coslat = cos((double)rlat);
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
  BC_CPREC rlat;
  if(fscanf(in,"%lf %lf %lf %lf %lf\n",
	    &(quake->lon),&(quake->lat),&(quake->depth),
	    &(quake->mag),&quake->tsec) == 5){
    if(quake->lon <0)
      quake->lon += 360;
    rlat = BC_DEG2RAD(quake->lat);
    quake->coslat = cos((double)rlat);
    quake->m0 = mag2mom(quake->mag);
    return 1;
  }else{
    return 0;
  }
}




/* tin, tou[6] symmetric storage  */
void rotate_vec6(BC_CPREC *tin, BC_CPREC *tout, BC_CPREC alpha, BC_CPREC beta, BC_CPREC gamma)
{
  BC_CPREC r[3][3],xin[3][3],xout[3][3];
  get_gen_rot(r,alpha,beta,gamma); /* get rotation matrix */
  sixsymtomat(tin,xin);		   /* convert to 3x3 */
  rotate_mat(xin,xout,r);	   /* rotate 3x3 */
  mattosixsym(xout,tout);	   /* convert to [6] */
}
/* convert symmetric [6] storage to [3][3] matrix */
void sixsymtomat(BC_CPREC *in6,BC_CPREC out[3][3])
{
  out[BC_R][BC_R] = in6[BC_RR];
  out[BC_R][BC_THETA]= out[BC_THETA][BC_R] = in6[BC_RT];
  out[BC_R][BC_PHI] = out[BC_PHI][BC_R] = in6[BC_RP];
  out[BC_THETA][BC_THETA] = in6[BC_TT];
  out[BC_THETA][BC_PHI] = out[BC_PHI][BC_THETA] = in6[BC_TP];
  out[BC_PHI][BC_PHI] = in6[BC_PP];
}
void mattosixsym(BC_CPREC in[3][3], BC_CPREC *out6)
{
  out6[BC_RR] = in[BC_R][BC_R] ;
  out6[BC_RT] = in[BC_R][BC_THETA];
  out6[BC_RP] = in[BC_R][BC_PHI];
  out6[BC_TT] = in[BC_THETA][BC_THETA];
  out6[BC_TP] = in[BC_THETA][BC_PHI];
  out6[BC_PP] = in[BC_PHI][BC_PHI];
}

/*
  
  obtain a general rotation matrix with angles alpha, beta, and gamma
  (given in degrees) as defined in Dahlen and Tromp, p. 921
  
*/
void get_gen_rot(BC_CPREC r[3][3],BC_CPREC alpha,
		 BC_CPREC beta, BC_CPREC gamma)
{
  BC_CPREC ralpha,sin_alpha,cos_alpha;
  BC_CPREC rbeta,sin_beta,cos_beta;
  BC_CPREC rgamma,sin_gamma,cos_gamma;

  ralpha=BC_DEG2RAD(alpha);
  sin_alpha = sin(ralpha);
  cos_alpha = cos(ralpha);

  rbeta=BC_DEG2RAD(beta);
  sin_beta = sin(rbeta);
  cos_beta = cos(rbeta);

  rgamma=BC_DEG2RAD(gamma);
  sin_gamma = sin(rgamma);
  cos_gamma = cos(rgamma);
  
  r[BC_R][BC_R] = cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma; 
  r[BC_R][BC_THETA] = sin_alpha*cos_beta*cos_gamma + cos_alpha*sin_gamma;
  r[BC_R][BC_PHI] = -sin_beta*cos_gamma;
  r[BC_THETA][BC_R] = -cos_alpha*cos_beta*sin_gamma -sin_alpha*cos_gamma;
  r[BC_THETA][BC_THETA] = -sin_alpha*cos_beta*sin_gamma +cos_alpha*cos_gamma;
  r[BC_THETA][BC_PHI] = sin_beta*sin_gamma;
  r[BC_PHI][BC_R] = cos_alpha*sin_beta;
  r[BC_PHI][BC_THETA] = sin_alpha*sin_beta;
  r[BC_PHI][BC_PHI] = cos_beta;
}

/* 
   rotate a 3x3 tensor using a general rotation matrix r 
   xout = r . xin . r^T
*/
void rotate_mat(BC_CPREC xin[3][3],BC_CPREC xout[3][3],
		BC_CPREC r[3][3])
{
  BC_CPREC tmp[3][3];
  int i,j,k;
  // calculate xin . r^T
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      tmp[i][j] = 0.0;
      for(k=0;k<3;k++)
	tmp[i][j] += xin[i][k] * r[j][k];
    }
  // calculate r . tmp
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      xout[i][j] = 0.0;
      for(k=0;k<3;k++)
	xout[i][j] += r[i][k] * tmp[k][j];
    }
}
/*  */
void print_quake_aki(FILE *out, struct qke quake)
{
  fprintf(out,"%9.4f %9.4f %6.2f %8.2f  %8.2f %8.2f %5.2f %8.2f %8.2f %.8e\n",
	  quake.lon,quake.lat,quake.depth,
	  quake.strike,quake.dip,quake.rake,
	  quake.mag,quake.lon,quake.lat,quake.tsec);
}
//
// X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp[dyn cm], newX, newY, time
//
void print_quake_cmt(FILE *out, struct qke quake)
{
  fprintf(out,"%9.4f %9.4f %6.2f\t%10.5f %10.5f %10.5f %10.5f %10.5f  %10.5f %3i %9.4f %9.4f %f\n",
	  quake.lon,quake.lat,quake.depth,
	  quake.m[BC_RR],quake.m[BC_TT],quake.m[BC_PP],quake.m[BC_RT],quake.m[BC_RP],quake.m[BC_TP],
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
	  quake.lon,quake.lat,quake.depth,
	  quake.strike,quake.dip,quake.rake,
	  quake.strike2,quake.dip2,quake.rake2,
	  quake.m0/pow(10,iexp),iexp+7,quake.lon,quake.lat,quake.tsec);
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
  catalog->minlon = 1000;
  catalog->maxlon = -1000;
  catalog->minlat = 1000;
  catalog->maxlat = -1000;
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

FILE *myopen(char *filename, char *rwmode, char *program)
{
  FILE *stream;
  stream = fopen(filename,rwmode);
  if(!stream){
    fprintf(stderr,"%s: can not open file \"%s\" for mode %s, exiting\n",
	    program,filename,rwmode);
    exit(-1);
  }
  return stream;
}
/* 

   compute distance on sphere, input in degree, return in km 

   cosines for lat1 and lat2 are precomputed

*/
BC_CPREC distance_geo(BC_CPREC lon1,BC_CPREC lat1,
		   BC_CPREC lon2,BC_CPREC lat2,
		   BC_CPREC coslat1, BC_CPREC coslat2)
{
  BC_CPREC tmp1,tmp2,tmp3;

  lon1/=BC_PIF;
  lat1/=BC_PIF;
  lon2/=BC_PIF;
  lat2/=BC_PIF;
  tmp1 = sin((lat1 - lat2)/2.0);
  tmp1 = tmp1 * tmp1;
  tmp2 = sin((lon1 - lon2)/2.0);
  tmp2 = tmp2 * tmp2;
  tmp2 *= coslat1;
  tmp2 *= coslat2;
  tmp3 = sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3)*BC_RADIUS;
}

BC_CPREC distance_cart(BC_CPREC x1,BC_CPREC y1,BC_CPREC x2, BC_CPREC y2)
{
  BC_CPREC tmp1,tmp2;
  tmp1 = x1-x2;
  tmp1*=tmp1;
  tmp2 = y1-y2;
  tmp2*=tmp2;
  return sqrt(tmp1+tmp2);
}
/* compute the distance between quake i of catalog one and quake j (or
   bin j) of catalog two, depending on cartesian coordinates or not */

BC_CPREC distance(struct cat *c1,struct cat *c2,int i,int j)
{
  if(c1->is_xy){
    return distance_cart(c1->quake[i].lon,c1->quake[i].lat,
			 c2->quake[j].lon,c2->quake[j].lat);
    
  }else{
    return distance_geo(c1->quake[i].lon,c1->quake[i].lat,
			c2->quake[j].lon,c2->quake[j].lat,
			c1->quake[i].coslat,c2->quake[j].coslat);
    
  }
}

void remove_trace(BC_CPREC *a)
{
  BC_CPREC t3;
  t3 = (a[BC_RR]  + a[BC_TT] + a[BC_PP])/3.0;
  a[BC_RR] -= t3;
  a[BC_PP] -= t3;
  a[BC_TT] -= t3;
}
/*  */
BC_CPREC mean_hor_strain(BC_CPREC *a)
{
  BC_CPREC s11,s12,s22,r,x1,x2,fms,sms;
  s11 =  a[BC_PP];
  s12 = -a[BC_TP];
  s22 =  a[BC_TT];
  x1 = (s11 + s22)/2.0;
  x2 = (s11 - s22)/2.0;
  r = x2 * x2 + s12 * s12;
  if(r > 0.0){
    r = sqrt(r);
    fms = x1 + r;
    sms = x1 - r;
  }else{
    fms = sms = x1;
  }
  return ((fms+sms)/2.0);
}


void get_index_vector(int **ind, int n, int random,long *seed)
{
  int i;
  *ind = (int *)realloc(*ind,n * sizeof(int));
  if(!ind){BC_MEMERROR("index");}
  if(random > 0 ){
    for(i=0;i<n;i++)
      *(*ind+i) = (int)(BC_RGEN(seed) * n);
  }else{
    for(i=0;i < n;i++)
      *(*ind+i) = i;
  }
}

BC_CPREC gauss_ran(long int *seed, BC_CPREC sigma)
{
  return sigma * gasdev(seed);
}

/* get gaussian distribution */
BC_CPREC gasdev(long int *seed)
{
  static int iset=0;
  static BC_CPREC gset;
  BC_CPREC fac,rsq,v1,v2;
  if  (iset == 0) {
    do {
      v1=2.0*BC_RGEN(seed)-1.0;
      v2=2.0*BC_RGEN(seed)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 5e-15
#define RNMX (1.0-EPS)
BC_CPREC ran2(long int *idum)
{
  int j;
  long k;
  static int ntabp7 = NTAB + 7;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  BC_CPREC temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=ntabp7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) 
    *idum += IM1;
  k=idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) 
    idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) 
    iy += IMM1;
  if ((temp=AM*iy) > RNMX) 
    return RNMX;
  else 
    return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/* 

moment conversions
 from magnitude to moment 
*/
BC_CPREC mag2mom(BC_CPREC mag)
{				     /* Hanks and Kanamori (1979) */
  return pow(10.0,3./2.*mag + 9.1); /* M0 in Nm */
}
BC_CPREC mom2mag(BC_CPREC mom)
{
  return ((2./3.)*(log10(mom) - 9.1)); /* in Nm */
}

BC_CPREC mag2pot(BC_CPREC mag)
{
  /* scalar potency */
  if(mag < 3.5){		// lower scaling
    return pow(10.0,1.45 * mag - 5.69);
  }else{ //
    return pow(10.0,1.08 * mag - 4.87);
  } 
}

/* input is a [6] tensor */
BC_CPREC scalar_mom(BC_CPREC *m)	/* compute scalar moment from tensor */
{
  return tensor6_norm(m) * M_SQRT1_2; /* norm/sqrt(2) */
}

void normalize_tens6(BC_CPREC *m6)	/* normalize a tensor given in 0...5
				   format */
{
  BC_CPREC norm;
  int i;
  norm = tensor6_norm(m6);
  for(i=0;i<6;i++)
    m6[i] /= norm;
  /* test */
  //fprintf(stderr,"old norm: %12.5e new norm: %12.5e\n",norm,tensor6_norm(m6));
}

/* tensor in [6] notation */
BC_CPREC tensor6_norm(BC_CPREC *m6)	/* plug in upper triangle */
{
  BC_CPREC t[3][3], ret;
  int i,j;
  tens6to3by3(m6,t);		/* convert to 3x3 */
  ret = 0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      ret += t[i][j] * t[i][j];
  return sqrt(ret);
}

void tens6to3by3(BC_CPREC *m, BC_CPREC t[3][3])
{
  t[BC_R][BC_R] = m[BC_RR];
  t[BC_THETA][BC_R] = t[BC_R][BC_THETA] = m[BC_RT];
  t[BC_THETA][BC_THETA] = m[BC_TT];
  t[BC_PHI][BC_THETA] = t[BC_THETA][BC_PHI] = m[BC_TP];
  t[BC_PHI][BC_PHI] = m[BC_PP];
  t[BC_R][BC_PHI] = t[BC_PHI][BC_R] = m[BC_RP];
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
   calculate gutenberg richter following Marzocci 2013

   m[nm] magnitudes, with dm spacing (e.g. 0.1)
   mmin: magnitude threshold (e.g. 2)

   returns b and error thereof, sb
*/
void calc_b_value_marzocci(BC_CPREC *m, int nm, BC_CPREC dm, BC_CPREC mmin,
			   BC_CPREC *b, BC_CPREC *sb)
{
  // M_LN10 = 2.30258509299405
  BC_CPREC mu,sum,std,p; 
  int i,nm_ac;
  for(nm_ac=0,mu=0.,i=0;i<nm;i++){ /* mean for magnitudes above
				   completeness */
    if(m[i] >= mmin){
      mu += m[i];
      nm_ac++;
    }
  }
  mu /= (BC_CPREC)nm_ac;
  
  if(dm > 0){			/* for magnitudes with dm spacing */
    /* Tinti and Mulargia (1987) approach */
    /* this doesn't quite work for sb, why? */
    p = 1. + dm/(mu-mmin); 	/* eq. 3.10 */
    *b = 1/(M_LN10*dm) * log(p);		    /* b estimate, eq. 3.9 */
    *sb = -(1-p)/(M_LN10*dm*sqrt((BC_CPREC)nm_ac*p)); /* error */
  }else{
    for(std=0.,i=0;i < nm_ac;i++){
      sum = m[i]-mu;
      std += sum*sum;
    }
    std /= ((BC_CPREC)nm_ac)*(((BC_CPREC)nm_ac)-1.);
    std = sqrt(std);
    *b = 1/(M_LN10*(mu-mmin));
    *sb = 2.3*(*b)*(*b)*std;
  }
}
/* 
   calculate gutenberg richter 

   m[nm] magnitudes, with dm spacing (e.g. 0.1)
   mmin: magnitude threshold (e.g. 2)

   returns b and error thereof, sb
*/
void calc_b_value_thomas(BC_CPREC *m, int nm, BC_CPREC dm, BC_CPREC mmin,
			 BC_CPREC *b, BC_CPREC *sb)
{
  // M_LN10 = 2.30258509299405
  BC_CPREC mu,sum,sum2,std; 
  int i,nm_ac;
  for(nm_ac=0.,sum=sum2=0.,i=0;i<nm;i++){	/* mean for magnitudes above
					   completeness */
    if(m[i] >= mmin){
      sum += m[i];
      sum2 += m[i]*m[i];
      nm_ac++;
    }
  }
  mu = sum/(BC_CPREC)nm_ac;		/* mean magnitude */
  std = std_quick(nm_ac, sum,sum2);
  /* Utsu (1966), Bender (1983) dm corrected equation, eq. 3.1 of
     Marzocci */
  *b =  ( 1. / ( mu - (mmin - dm/2)) )/M_LN10; /* b value (1/ln(10) = log10(e)) */
  /* standard deviation */
  *sb = 2.3 * std/sqrt((BC_CPREC)nm_ac) * (*b)*(*b); /* Shi and Bolt (1982), eq. 2.5 */
}
/* 
   m[nm] magnitudes, mmin: magnitude threshold

   simple max likelihood from Aki, this is not a good idea

   dm unused!
*/
void  calc_b_value_ml(BC_CPREC *m, int nm, BC_CPREC dm, BC_CPREC mmin,BC_CPREC *b, BC_CPREC *sb)
{
  // M_LN10 = 2.30258509299405
  BC_CPREC mu; 
  int i,nm_ac;
  for(nm_ac=0,mu=0.,i=0;i<nm;i++){	/* mean for magnitudes above
					   completeness */
    if(m[i] > mmin){
      mu += m[i];
      nm_ac++;
    }
  }
  mu /= (BC_CPREC)nm_ac;
  
  *b = 1/(mu-mmin)/M_LN10;	/* eq. 2.3 from Marzocci */
  *sb = *b/sqrt((BC_CPREC)nm_ac);	/* eq. 2.4  */
}

/* 
   quick and dirty way to compute standard deviationq 

*/
BC_CPREC std_quick(int n, BC_CPREC sum, BC_CPREC sum2)
{
  double nf;
  nf = (double)n;
  return sqrt((nf * (double)sum2 - (double)sum*(double)sum) / (nf*(nf-1)));
}
/* code from Andy Michael, as the SLICK package */

/* switcher routine 

   find alternate fault plane from given strike(dir, CW from N), dip
   and rake
*/
void find_alt_plane(BC_CPREC ddir1,BC_CPREC dip1,BC_CPREC rake1,
		    BC_CPREC *ddir2,BC_CPREC *dip2,BC_CPREC *rake2)
{
  BC_CPREC z,z2,z3,s1,s2,s3,sin_z,cos_z,sin_z2,cos_z2,sin_z3,cos_z3;
  BC_CPREC  n1,n2,h1,h2;
  
  z=ddir1/BC_PIF;
  if(dip1==90)
    dip1=89.99999;
  z2=dip1/BC_PIF;
  z3=rake1/BC_PIF;

  sin_z  = sin(z); cos_z  = cos(z);
  sin_z2 = sin(z2);cos_z2 = cos(z2);
  sin_z3 = sin(z3);cos_z3 = cos(z3);
  
  /* slick vector in plane 1 */
  s1= -cos_z3*cos_z-sin_z3*sin_z*cos_z2;
  s2=  cos_z3*sin_z-sin_z3*cos_z*cos_z2;
  s3= sin_z3*sin_z2;

  n1=sin_z*sin_z2;  /* normal vector to plane 1 */
  n2=cos_z*sin_z2;

  h1= -s2; /* strike vector of plane 2 */
  h2= s1;
  /* note h3=0 always so we leave it out */
  stridip(s2,s1,s3,&z,&z2);
  z+= 90.;
  *ddir2=z;
  ranger(ddir2);
  *dip2=z2;
  z= h1*n1 + h2*n2;
  z/= sqrt(h1*h1 + h2*h2);
  z=acos(z);
  if(s3>=0)
    *rake2= z*BC_PIF;
  else
    *rake2= -z*BC_PIF;
}

void ranger(BC_CPREC *z)
/* makes z in 0 to 360 */
{
  while(*z >= 360)
    *z-= 360;
  while(*z < 0)
    *z += 360;
}

void stridip(BC_CPREC n,BC_CPREC e,BC_CPREC u,BC_CPREC *strike,BC_CPREC *dip)
/* finds the strike and dip of a plane given its normal */
/* vector, output is in degrees north of east and then  */
/* uses a right hand rule for the dip of the plane */
{
  BC_CPREC x;
  if(u <0.) {
    n= -n;
    e= -e;
    u= -u;
  }
  *strike=atan2(e,n)*BC_PIF;
  *strike= *strike-90.;
  if(*strike < 0.)
    *strike+= 360.;
  if(*strike > 360.)
    *strike-= 360.;
  x=sqrt(n*n+e*e);   /* x is the horizontal magnitude */
  *dip=atan2(x,u)*BC_PIF;
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

/* 
   
   compute the stress tensor as expected from a Michael (1984, 1987)
   inversion for each bin, using the assigned weights
   
*/
void calc_stress_tensor_for_kbins(struct cat *catalog)
{
  struct kostrov_sum *kostrov;
  int i,j,n,iq;
  BC_CPREC *angles,*weights;
  kostrov = catalog->sum;
  angles = (BC_CPREC *)malloc(6*sizeof(BC_CPREC));
  weights = (BC_CPREC *)malloc(sizeof(BC_CPREC));

  for(i=0;i < kostrov->nxny;i++){	
    /* 
       
       for each bin, pass all events with their two fault planes

    */
    if(kostrov->bin[i].n > BC_NQUAKE_LIM_FOR_STRESS){	/* five parameters, at least two events */
      n=0;
      for(j=0;j < kostrov->bin[i].n;j++){ /* assign */

	weights[n] = kostrov->bin[i].weight[j]; /* weight of this
						   earthquake in the
						   bin */
	
	iq = (int)kostrov->bin[i].quake[j]; /* number of this earthquake */
	
	angles[n*6+0] = catalog->quake[iq].strike/BC_PIF;
	angles[n*6+1] = catalog->quake[iq].dip/BC_PIF;
	angles[n*6+2] = catalog->quake[iq].rake/BC_PIF;

	angles[n*6+3] = catalog->quake[iq].strike2/BC_PIF;
	angles[n*6+4] = catalog->quake[iq].dip2/BC_PIF;
	angles[n*6+5] = catalog->quake[iq].rake2/BC_PIF;
	n++;
	angles = (BC_CPREC *)realloc(angles,6*(n+1)*sizeof(BC_CPREC));
	weights = (BC_CPREC *)realloc(weights,(n+1)*sizeof(BC_CPREC));
      }
      solve_stress_michael(n,angles,weights,
			   kostrov->bin[i].s,kostrov->bin[i].ds,&catalog->seed);
    }
  }
  free(angles);free(weights);
}
/* 
   solve for the five component stress tensor a la Andy Michael angles
   are given in radians
   
*/
void solve_stress_michael(int nquakes, BC_CPREC *angles,BC_CPREC *weights,
			  BC_CPREC *stress, BC_CPREC *sig_stress,
			  long int *seed)
{
  const int npar = BC_MICHAEL_NPAR;
  const int ndim = BC_NDIM;
  const int nrandom_limit  = BC_MICHAEL_NMC; /* monte carlo simulations (2000
					     good number for 95%
					     confidence?)*/
  BC_BOOLEAN proceed;
  int nobs,iquake,nrandom,i;
  BC_CPREC ind_stress[6],*slick,*amat,tot_stress[6],tot_stress2[6],snorm;

  slick = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim);
  amat = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim*npar);
  proceed  = BC_TRUE;nrandom=0;
  for(i=0;i < 6;i++)
    tot_stress[i] = tot_stress2[i] = 0.0;
  do{
    /* converted from Andy Michael's slick routine */
    nobs = 0;
    for(iquake=0;iquake < nquakes;iquake++){
      /* randomly assign fault planes */
      if(BC_RGEN(seed) >= 0.5)
	michael_assign_to_matrix(angles[iquake*6+0],angles[iquake*6+1],angles[iquake*6+2],
				 &nobs,&slick,&amat);
      else			/* alternate FP */
	michael_assign_to_matrix(angles[iquake*6+3],angles[iquake*6+4],angles[iquake*6+5],
				 &nobs,&slick,&amat);
    }  /* end of data assignment loop */
    /* 
       michael least squares solve, augmented by weights (will
       override amat and slick)
    */
    michael_solve_lsq(npar,ndim,nobs,amat,slick,weights,ind_stress);
    for(i=0;i<6;i++){
      tot_stress[i]  += ind_stress[i];
      tot_stress2[i] += ind_stress[i] * ind_stress[i];
    }
    nrandom++;
    if(nrandom>nrandom_limit)
      proceed = BC_FALSE;
  }while(proceed);
  for(i=0;i<6;i++){
    sig_stress[i] = std_quick(nrandom,tot_stress[i],tot_stress2[i]);
    stress[i] = tot_stress[i]/(BC_CPREC)nrandom;
  }
  /* normalize tensor */
  snorm = tensor6_norm(stress);
  for(i=0;i<6;i++){
    stress[i]     /= snorm;
    sig_stress[i] /= snorm;
  }
  if(0)
    fprintf(stderr,"stress_solve: avg: %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g)\n",
	    stress[0],sig_stress[0],stress[1],sig_stress[1],stress[2],sig_stress[2],
	    stress[3],sig_stress[3],stress[4],sig_stress[4],stress[5],sig_stress[5]);
  free(amat);
  free(slick);
}

/* 
   given assigned slick and amat vectors and matrices, use Andy
   Michael's least squares solver

   this will overwrite amat and slick

*/
void michael_solve_lsq(int npar,int ndim, int nobs, BC_CPREC *amat, BC_CPREC *slick,BC_CPREC *weights,
		       BC_CPREC *stress)
{
  int m,i,j,k,moff;
  BC_CPREC *a2,*cc,sigma,lstress[6],w;

  m = nobs * ndim;
  a2 = (BC_CPREC *)malloc(sizeof(BC_CPREC)*npar*npar);
  cc = (BC_CPREC *)malloc(sizeof(BC_CPREC)*npar);

  /* rescale with weights */
  for(i=0;i < nobs;i++){
    w = weights[i];		/* weight for each observation */
    for(j=0;j < ndim;j++){
      moff = i*ndim+j;
      slick[moff] *= w;
      for(k=0;k < npar;k++)
	amat[moff*npar+k] *= w;
    }
  }
    
  /* solve equations via linear least squares */
  /*               0    1   2  3    4   5   */
  /* lstress is in xx, xy, xz, yy, yz, zz format E, N, U */
  /* 
     solve amat.ltress = slick in a least squares sense 
  */
  michael_leasq(amat,npar,m,lstress,slick,a2,cc,&sigma); 
  /* 
     fix zz element by using trace = 0 
  */
  lstress[5]= -(lstress[0]+lstress[3]);
  /* 
     reassign from ENU to my spherical coordinates 
  */
  stress[BC_RR] =  lstress[5];	/*  zz */
  stress[BC_RT] = -lstress[4];	/* -zy */
  stress[BC_RP] =  lstress[2];	/*  zx */
  stress[BC_TT] =  lstress[3];	/*  yy */
  stress[BC_TP] = -lstress[1];	/* -yx */
  stress[BC_PP] =  lstress[0];	/*  xx */
  
  free(a2);free(cc);
}

/* 

   build the design matrix and data vector (slickenside)

   z,z2,z3 are 

   strike, dip, rake angles (aki and richards) in radians

   from Andy Michael's slick.c routine

   NOTE: 
   
   this is converted from Michael's direction (i.e. strike + pi/2),
   dip, rake to Aki and Richards
   
*/
void michael_assign_to_matrix(BC_CPREC strike_rad,BC_CPREC dip_rad,BC_CPREC rake_rad,
			      int *nobs,BC_CPREC **slick,BC_CPREC **amat)
{
  int j;
  BC_CPREC sin_z,cos_z,sin_z2,cos_z2,sin_z3,cos_z3,n1,n2,n3,n12,n22,n32,z,z2,z3;
  const int ndim = BC_NDIM;
  const int npar = BC_MICHAEL_NPAR;
  
  z  = strike_rad + M_PI_2;	/* convert strike of fault to
				   direction of dip */
  
  z2 = dip_rad;			/* rest is same */
  z3 = rake_rad;
  
  /*  */
  sin_z =  sin(z);  cos_z  = cos(z);
  sin_z2 = sin(z2); cos_z2 = cos(z2);
  sin_z3 = sin(z3); cos_z3 = cos(z3);
    
  n1=sin_z*sin_z2;  /* normal vector to fault plane */
  n2=cos_z*sin_z2;
  n3=cos_z2;

  n12 = n1*n1;
  n22 = n2*n2;
  n32 = n3*n3;

  j = (*nobs) * ndim;
  
  /* slickenside vector calculation */
  *(*slick+j)=  -cos_z3*cos_z-sin_z3*sin_z*cos_z2;
  *(*slick+j+1)= cos_z3*sin_z-sin_z3*cos_z*cos_z2;
  *(*slick+j+2)= sin_z3*sin_z2;
  
  /* find the matrix elements */
  *(*amat+j*npar+0)    = n1-n12*n1+n1*n32;
  *(*amat+j*npar+1)    = n2-2.*n12*n2;
  *(*amat+j*npar+2)    = n3-2.*n12*n3;
  *(*amat+j*npar+3)    = -n1*n22+n1*n32;
  *(*amat+j*npar+4)    = -2.*n1*n2*n3;
  
  *(*amat+(j+1)*npar+0)= -n2*n12+n2*n32;
  *(*amat+(j+1)*npar+1)= n1-2.*n1*n22;
  *(*amat+(j+1)*npar+2)= -2.*n1*n2*n3;
  *(*amat+(j+1)*npar+3)= n2-n22*n2+n2*n32;
  *(*amat+(j+1)*npar+4)= n3-2.*n22*n3;
  
  *(*amat+(j+2)*npar+0)= -n3*n12-n3+n32*n3;
  *(*amat+(j+2)*npar+1)= -2.*n1*n2*n3;
  *(*amat+(j+2)*npar+2)= n1-2.*n1*n32;
  *(*amat+(j+2)*npar+3)= -n3*n22-n3+n32*n3;
  *(*amat+(j+2)*npar+4)= n2-2.*n2*n32;

  /* increment counter and make room for mroe */
  *nobs = *nobs + 1;
  /*  */
  *slick = (BC_CPREC *)realloc(*slick,sizeof(BC_CPREC)*ndim*     (*nobs+1));
  if(!slick)BC_MEMERROR("assign to Michael mat, slick");
  *amat = (BC_CPREC *) realloc(*amat, sizeof(BC_CPREC)*ndim*npar*(*nobs+1));
  if(!amat)BC_MEMERROR("assign to Michael mat, amat");
}
