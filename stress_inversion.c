#include "catalog.h"
/* 

   bunch of stress inversion stuff

   some routines modified from code by Andy Michael 

   Michael, A.J., 1984. Determination of stress from slip data: Faults
   and folds, J. Geophys. Res. 89, 11.517-11.526.

   Michael, A.J., 1987. Use of focal mechanisms to determine stress: A
   control study, J. Geophys. Res. 92, 357-368.

   some routines are modified from STRESSINVERSE
   https://www.ig.cas.cz/en/stress-inverse/ 

   Vavryčuk, V., 2014. Iterative joint inversion for stress and fault
   orientations from focal mechanisms, Geophysical Journal
   International, 199, 69-77, doi: 10.1093/gji/ggu224.


 */

/* 
   
   compute the stress tensor as expected from a Michael (1984, 1987)
   inversion for each bin, using the assigned weights, and possibly add Vavrychuk style
   
*/
void calc_stress_tensor_for_kbins(struct cat *catalog)
{
  struct kostrov_sum *kostrov;
  int i,j,j6,iq,opt;
  BC_CPREC *angles,*weights,ainst[2],shape_ratio,fopt,minst,fmin,fmax,finc;
#ifdef DEBUG
  int k;
  BC_CPREC sdev[2];
#endif
  kostrov = catalog->sum;
  angles = (BC_CPREC *)malloc(6*sizeof(BC_CPREC));
  weights = (BC_CPREC *)malloc(sizeof(BC_CPREC));

  for(i=0;i < kostrov->nxny;i++){	
    /* 
       
       for each bin, pass all events with their two fault planes

    */
    if(kostrov->bin[i].n >= kostrov->nmin){	/* five parameters, at least nmin events */
#ifdef DEBUG
      fprintf(stderr,"bin %05i: %g, %g working on stress inversion, N %i\n",
	      i, kostrov_bdlon(i,kostrov),kostrov_bdlat(i,kostrov),kostrov->bin[i].n);
#endif
      angles = (BC_CPREC *)realloc(angles,kostrov->bin[i].n*6*sizeof(BC_CPREC));
      weights = (BC_CPREC *)realloc(weights,kostrov->bin[i].n*sizeof(BC_CPREC));

      for(j=j6=0;j < kostrov->bin[i].n;j++,j6+=6){ /* assign */
	weights[j] = kostrov->bin[i].weight[j]; /* weight of this
						   earthquake in the
						   bin */
	iq = (int)kostrov->bin[i].quake[j]; /* number of this earthquake */
	assign_quake_angles((catalog->quake+iq),(angles+j6));
#ifdef DEBUG
	for(k=0;k<6;k++)
	  if(!finite(angles[j6+k])){
	    fprintf(stderr,"calc_stress_tensor_for_kbins: ERROR: quake %i angle %i not finite: %g\n",
		    iq,k,angles[j6+k]);
	    exit(-1);
	  }
#endif
      }
      /* randomized Michael: point estimate s and Monte Carlo uncertainty ds */
      solve_stress_michael_random_sweep(kostrov->bin[i].n,angles,weights,kostrov->bin[i].s,
					kostrov->bin[i].ds,&catalog->seed,BC_MICHAEL_RSWEEP_MAX);
      /* instability of the Michael stress, MATLAB/Vavrycuk convention,
	 consistent with stress_inversion_mstyle */
      mstyle_average_instability(kostrov->bin[i].n,angles,weights,BC_FRIC_DEF,kostrov->bin[i].s,ainst);
      kostrov->bin[i].inst[0] = ainst[0];
      
#ifdef DEBUG
      calc_misfits_from_single_angle_set(kostrov->bin[i].s,angles,kostrov->bin[i].n, sdev);
      fprintf(stderr,"srandom:          %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tsdev: %7.5f %7.5f\tinst: %7.5f %7.5f\n",
	      kostrov->bin[i].s[0],kostrov->bin[i].s[1],kostrov->bin[i].s[2],kostrov->bin[i].s[3],
	      kostrov->bin[i].s[4],kostrov->bin[i].s[5],
	      1-sdev[0],1-sdev[1],ainst[0],ainst[1]);
#endif
      if(catalog->use_friction_solve){
	/* Vavrycuk (mstyle), default fixed friction -> def_s, inst[1] */
	stress_inversion_mstyle(kostrov->bin[i].n,angles,weights,BC_FRIC_DEF,BC_FRIC_DEF,1.0,
				BC_VI_ITER,BC_VI_NREAL,&catalog->seed,
				kostrov->bin[i].def_s,&shape_ratio,&fopt,&minst,NULL,
				BC_STRESS_NORM_TENSOR);
	kostrov->bin[i].inst[1] = minst;
	/* best friction by scan -> best_s, best_fric, inst[2].
	   scan range mirrors the old optimize flag (use_friction_solve-1):
	   1 -> 0..1, else -> 0.2..0.8 */
	opt = catalog->use_friction_solve - 1;
	if(opt == 1){
	  fmin = 0.0; fmax = 1.0; finc = BC_FRIC_SCAN_INC;
	}else{
	  fmin = 0.2; fmax = 0.8; finc = BC_FRIC_SCAN_INC2;
	}
	stress_inversion_mstyle(kostrov->bin[i].n,angles,weights,fmin,fmax,finc,
				BC_VI_ITER,BC_VI_NREAL,&catalog->seed,
				kostrov->bin[i].best_s,&shape_ratio,&(kostrov->bin[i].best_fric),
				&minst,NULL,BC_STRESS_NORM_TENSOR);
	kostrov->bin[i].inst[2] = minst;
      }
    }
  }
  free(angles);free(weights);
}


/* 
   solve for the five component stress tensor a la Andy Michael 

   angles [6 = 2 times 3] are given in radians TWO SETS OF ANGLES


   michael_rsweep_max is by default something like 50000, but can be
   zero to not do any random sweeps (there is also
   solve_stress_michael_specified_plane which does ~the same thing)
   


   input:

   nquakes number of events
   angles strike1/dip1/rake1 strik2/dip2/rake2
   weights

   output:

   stress[6] tensor
   sig_stress[6] tensor
 
*/
void solve_stress_michael_random_sweep(int nquakes, BC_CPREC *angles,BC_CPREC *weights,
				       BC_CPREC *stress, BC_CPREC *sig_stress,
				       long int *seed, int michael_rsweep_max)
{
  const int npar = BC_MICHAEL_NPAR;
  const int ndim = BC_NDIM;
  const int nrandom_limit  = BC_MICHAEL_NMC; /* monte carlo simulations (2000 good
						number for 95%
						confidence?)*/
  BC_BOOLEAN proceed,acc_bail,iter_warned;
  static BC_BOOLEAN warned = BC_FALSE;
  int nobs,iquake,nrandom,i,iquake6,icheck,j;
  BC_CPREC ind_stress[6],*slick,*amat,tot_stress[6],tot_stress2[6],snorm,last_stress[6],this_stress[6],ds,tmp;
  size_t ssize = 6*sizeof(BC_CPREC);
  BC_CPREC bail_acc_squared = BC_MICHAEL_RACC*BC_MICHAEL_RACC;
 
  if(michael_rsweep_max == 0){	/* running the stuff below only once
				   gives the same answer, do like so
				   to avoid having an if statement */
    solve_stress_michael_specified_plane(nquakes,angles,weights,stress,BC_TRUE);
    for(i=0;i < 6;i++)
      sig_stress[i]=NAN;
  }else{
    slick = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim);
    amat = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim*npar);
    for(i=0;i < 6;i++){
      tot_stress[i] = tot_stress2[i] = last_stress[i] = 0.0;
    }
    iter_warned = BC_FALSE;
    proceed  = BC_TRUE;
    nrandom=0;
    
    do{
      acc_bail = BC_FALSE;
      /* converted from Andy Michael's slick routine */
      nobs = 0;
      for(iquake=iquake6=0;iquake < nquakes;iquake++,iquake6+=6){
	/* randomly assign fault planes */
	if(BC_RGEN(seed) >= 0.5)
	  michael_assign_to_matrix((angles+iquake6),  &nobs,&slick,&amat);
	else			/* alternate FP */
	  michael_assign_to_matrix((angles+iquake6+3),&nobs,&slick,&amat);
      }  /* end of data assignment loop */
      /* 
	 michael least squares solve, augmented by weights (will
	 override amat and slick)
      */
      michael_solve_lsq(npar,ndim,nobs,amat,slick,weights,ind_stress);
      for(i=icheck=0;i<6;i++){
	if(!finite(ind_stress[i])){
	  if(!iter_warned){
	    fprintf(stderr,"ssm random_sweep: WARNING: ind_stress not finite\n");
	    for(iquake=iquake6=0;iquake < nquakes;iquake++,iquake6+=6){
	      for(j=0;j<6;j++)
		fprintf(stderr,"%6.3f ",(float)*(angles+iquake6+j));
	      fprintf(stderr,"\n");
	    }
	  }
	  iter_warned = BC_TRUE;
	}else{
	  icheck++;
	}      
      }
      if(icheck == 6){		/* all finite */
	for(i=0;i<6;i++){
	  tot_stress[i]  += ind_stress[i];
	  tot_stress2[i] += ind_stress[i] * ind_stress[i];
	}
	nrandom++;
      }
      if(nrandom%100==0){
	/* every 100 iterations, check for finiteness and convergence */
	ds = 0;
	for(i=0;i<6;i++){
	  if(!finite(tot_stress[i])){
	    fprintf(stderr,"ssm random_sweep: ERROR: tot_stress not finite\n");
	    exit(-1);
	  }
	  this_stress[i] = tot_stress[i]/(BC_CPREC)nrandom;
	  tmp = this_stress[i] - last_stress[i];
	  ds += tmp*tmp;
	}
	if(ds < bail_acc_squared)
	  acc_bail = BC_TRUE;
#ifdef SUPER_DEBUG
	fprintf(stderr,"ssm random_sweep: iter %05i s %12g %12g %12g %12g %12g %12g: diff: %12.5e\n",nrandom,
		this_stress[0], this_stress[1], this_stress[2], this_stress[3], this_stress[4],this_stress[5],sqrt(ds));
#endif
	memcpy(last_stress,this_stress,ssize);
      }
      if(acc_bail || (nrandom > BC_MICHAEL_RSWEEP_MAX))
	proceed = BC_FALSE;
    }while(proceed);
    if(nrandom > BC_MICHAEL_RSWEEP_MAX){
      if(!warned){
	fprintf(stderr,"ssm random_sweep: WARNING: bailed on max sweeps (%i) not accuracy (%12.5e) on at least one bin\n",
		nrandom,sqrt(ds));
	warned = BC_TRUE;
      }
    }
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
    /*
      fprintf(stderr,"stress_solve: avg: %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g)\n",
      stress[0],sig_stress[0],stress[1],sig_stress[1],stress[2],sig_stress[2],
      stress[3],sig_stress[3],stress[4],sig_stress[4],stress[5],sig_stress[5]);
    */
    free(amat);
    free(slick);
  }
 
 
}



/* 

   Michael, just one set of angles to use - like slick.c (but note
   that we are using the strike azimuth, not dip azimuth, for the
   angles)

*/
void solve_stress_michael_specified_plane(int nquakes, BC_CPREC *angles,
					  BC_CPREC *weights,BC_CPREC *stress,
					  BC_BOOLEAN normalize)
{
  const int npar = BC_MICHAEL_NPAR;
  const int ndim = BC_NDIM;
  int nobs,iquake,nrandom,i,iquake6;
  BC_CPREC *slick,*amat,lsdev[2];
  BC_CPREC m_smat[3][3],snorm;
  
  slick = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim);
  amat = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim*npar);
  /* converted from Andy Michael's slick routine */                                        
  nobs = 0;
  for(iquake=iquake6=0;iquake < nquakes;iquake++,iquake6+=6)
    michael_assign_to_matrix((angles+iquake6),&nobs,&slick,&amat);
  /*  */
  michael_solve_lsq(npar,ndim,nobs,amat,slick,weights,stress);
  if(normalize)
    normalize_tens6(stress,stress);
  free(amat);                    
  free(slick);
}


/* 
   given assigned slick and amat vectors and matrices, use Andy
   Michael's least squares solver

   this will overwrite amat and slick

*/
// debug the matrix assemble and solution?
//#define BC_SDD_DEBUG
void michael_solve_lsq(int npar,int ndim, int nobs, BC_CPREC *amat,
		       BC_CPREC *slick,BC_CPREC *weights,
		       BC_CPREC *stress)
{
  int m,i,j,k,moff;
  BC_CPREC *a2,*cc,sigma,lstress[6],w;
#ifdef BC_SDD_DEBUG
  FILE *out1,*out2;
#endif
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
#ifdef BC_SDD_DEBUG
  out1 = fopen("a.dat","w");
  out2 = fopen("b.dat","w");
  for(i=0;i < nobs;i++){
    for(j=0;j < ndim;j++){
      moff = i*ndim+j;
      fprintf(out2,"%15.7e\n",slick[moff]);
      for(k=0;k < npar;k++)
	fprintf(out1,"%15.7e ",amat[moff*npar+k]);
      fprintf(out1,"\n");
    }
  }
  fclose(out1);  fclose(out2);
#endif
    
  /* solve equations via linear least squares */
  /*               0    1   2  3    4   5   */
  /* lstress is in xx, xy, xz, yy, yz, zz format E, N, U */
  /* 
x     solve amat.ltress = slick in a least squares sense 
  */
  michael_leasq(amat,npar,m,lstress,slick,a2,cc,&sigma); 
#ifdef BC_SDD_DEBUG
  out1 = fopen("x.dat","w");
  for(k=0;k < npar;k++)
    fprintf(out1,"%15.7e\n",lstress[k]);
  fclose(out1);
#endif 
  /* 
     we solve for five components, assign sixth (zz) element by using
     the trace = 0 constraint
  */
  lstress[5]= -(lstress[0]+lstress[3]);

  /* 
     reassign from Michael cooridnate system (xyz = ENU) to my
     spherical coordinates
  */
  stress[BC_RR] =  lstress[5];	/*  zz =  UU */
  stress[BC_RT] = -lstress[4];	/* -yz = -NU */
  stress[BC_RP] =  lstress[2];	/*  xz =  EU */
  stress[BC_TT] =  lstress[3];	/*  yy =  NN */
  stress[BC_TP] = -lstress[1];	/* -yx = -EN */
  stress[BC_PP] =  lstress[0];	/*  xx =  EE */
  
  free(a2);free(cc);
}
/* 
   convert my polar stress vector to a matrix in Michael format, ENU 
*/
void my6stress2m3x3(BC_CPREC *svec,BC_CPREC msmat[3][3])
{
  msmat[0][0]            =  svec[BC_PP];
  msmat[0][1]=msmat[1][0]= -svec[BC_TP];
  msmat[0][2]=msmat[2][0]=  svec[BC_RP];
  msmat[1][1]            =  svec[BC_TT];
  msmat[1][2]=msmat[2][1]= -svec[BC_RT];
  msmat[2][2]            =  svec[BC_RR];
}
/* 

   build the design matrix and data vector (slickenside)

   z,z2,z3 are 

   angles[3] strike, dip, rake angles (Aki and Richards) convention IN RADIANS

   from Andy Michael's slick.c routine

   NOTE: 
   
   this is converted from Michael who uses 

   dip direction (i.e. strike + pi/2), dip, and rake 
   
   From Michael's description:

   "The first number on each line is the dip direction of the fault
    plane in degrees East of North.  The second number on each line is the
    dip of the fault plane.  The third number on each line is the rake of
    the fault, such that 0 is left lateral motion, 90 is thrust motion, 
    180 is right lateral motion, and 270 is normal faulting."

   else, this checks out against SATSI
   
*/
void michael_assign_to_matrix(BC_CPREC *angles,int *nobs,BC_CPREC **slick,BC_CPREC **amat)
{
  BC_CPREC sin_z,cos_z,sin_z2,cos_z2,sin_z3,cos_z3,n1,n2,n3,n12,n22,n32;
  int ind,j;
  const int ndim = BC_NDIM;
  const int npar = BC_MICHAEL_NPAR;
  
  /* strike */
  /* 
     using this convention, A, b, and x are identical with Michael
     slick.c and x solution checks out for both as actual least
     squares solution
  */
  //sincos(angles[0],         &sin_z, &cos_z); 
  sincos(angles[0]+ M_PI_2, &sin_z, &cos_z); /* convert from A&R
						strike azimuth to dip
						azimuth which is
						strike + 90 = strike + pi/2 */
  /* dip */
  sincos(angles[1],         &sin_z2,&cos_z2);
  /* rake */
  sincos(angles[2],         &sin_z3,&cos_z3);
    
    
  n1=sin_z*sin_z2;  /* normal vector to fault plane */
  n2=cos_z*sin_z2;
  n3=cos_z2;
  /* squared normal vector components */
  n12 = n1*n1;
  n22 = n2*n2;
  n32 = n3*n3;

  j = (*nobs) * ndim;
  
  /* slickenside vector calculation */
  *(*slick+j)=  -cos_z3*cos_z-sin_z3*sin_z*cos_z2;
  *(*slick+j+1)= cos_z3*sin_z-sin_z3*cos_z*cos_z2;
  *(*slick+j+2)= sin_z3*sin_z2;
 
  /* find the matrix elements */
  ind = j*npar;
  *(*amat+ind+0)    = n1-n12*n1+n1*n32;
  *(*amat+ind+1)    = n2-2.*n12*n2;
  *(*amat+ind+2)    = n3-2.*n12*n3;
  *(*amat+ind+3)    = -n1*n22+n1*n32;
  *(*amat+ind+4)    = -2.*n1*n2*n3;
  //fprintf(stderr,"%8.4f %8.4f %8.4f %8.4f %8.4f\n", *(*amat+ind+0), *(*amat+ind+1), *(*amat+ind+2), *(*amat+ind+3), *(*amat+ind+4));
  ind += npar;
  *(*amat+ind+0)= -n2*n12+n2*n32;
  *(*amat+ind+1)= n1-2.*n1*n22;
  *(*amat+ind+2)= -2.*n1*n2*n3;
  *(*amat+ind+3)= n2-n22*n2+n2*n32;
  *(*amat+ind+4)= n3-2.*n22*n3;
  //fprintf(stderr,"%8.4f %8.4f %8.4f %8.4f %8.4f\n", *(*amat+ind+0), *(*amat+ind+1), *(*amat+ind+2), *(*amat+ind+3), *(*amat+ind+4));
  ind += npar;
  *(*amat+ind+0)= -n3*n12-n3+n32*n3;
  *(*amat+ind+1)= -2.*n1*n2*n3;
  *(*amat+ind+2)= n1-2.*n1*n32;
  *(*amat+ind+3)= -n3*n22-n3+n32*n3;
  *(*amat+ind+4)= n2-2.*n2*n32;
  //fprintf(stderr,"%8.4f %8.4f %8.4f %8.4f %8.4f\n", *(*amat+ind+0), *(*amat+ind+1), *(*amat+ind+2), *(*amat+ind+3), *(*amat+ind+4));
  /* increment counter and make room for more */
  *nobs = *nobs + 1;
  /*  */
  *slick = (BC_CPREC *)realloc(*slick,sizeof(BC_CPREC)*ndim*     (*nobs+1));
  if(!slick)BC_MEMERROR("assign to Michael mat, slick");
  *amat = (BC_CPREC *) realloc(*amat, sizeof(BC_CPREC)*ndim*npar*(*nobs+1));
  if(!amat)BC_MEMERROR("assign to Michael mat, amat");
}
