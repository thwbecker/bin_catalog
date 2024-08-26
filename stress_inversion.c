#include "catalog.h"
/* 

   bunch of stress inversion stuff


   some routines are modified from STRESSINVERSE
   https://www.ig.cas.cz/en/stress-inverse/ 

   Vavryčuk, V., 2014. Iterative joint inversion for stress and fault
   orientations from focal mechanisms, Geophysical Journal
   International, 199, 69-77, doi: 10.1093/gji/ggu224.


 */

/* 
   
   compute the stress tensor as expected from a Michael (1984, 1987)
   inversion for each bin, using the assigned weights
   
*/
void calc_stress_tensor_for_kbins(struct cat *catalog)
{
  struct kostrov_sum *kostrov;
  int i,j,j6,iq,rsweep,osweep;
  BC_CPREC *angles,*weights,ainst[2];
  BC_BOOLEAN optim_verbose;
#ifdef DEBUG
  BC_CPREC sdev[2];
#endif
  kostrov = catalog->sum;
  angles = (BC_CPREC *)malloc(6*sizeof(BC_CPREC));
  weights = (BC_CPREC *)malloc(sizeof(BC_CPREC));

  for(i=0;i < kostrov->nxny;i++){	
    /* 
       
       for each bin, pass all events with their two fault planes

    */
    if(kostrov->bin[i].n > BC_NQUAKE_LIM_FOR_STRESS){	/* five parameters, at least two events */
#ifdef DEBUG
      fprintf(stderr,"bin %05i: working on stress inversion, N %i\n",i,kostrov->bin[i].n);
#endif
      angles = (BC_CPREC *)realloc(angles,kostrov->bin[i].n*6*sizeof(BC_CPREC));
      weights = (BC_CPREC *)realloc(weights,kostrov->bin[i].n*sizeof(BC_CPREC));

      for(j=j6=0;j < kostrov->bin[i].n;j++,j6+=6){ /* assign */
	weights[j] = kostrov->bin[i].weight[j]; /* weight of this
						   earthquake in the
						   bin */
	iq = (int)kostrov->bin[i].quake[j]; /* number of this earthquake */
	assign_quake_angles((catalog->quake+iq),(angles+j6));
      }
      solve_stress_michael_random_sweep(kostrov->bin[i].n,angles,weights,kostrov->bin[i].s,kostrov->bin[i].ds,&catalog->seed);
      calc_average_instability(kostrov->bin[i].n,angles,weights,BC_FRIC_DEF,kostrov->bin[i].s,ainst); /* instability of default stress */
      kostrov->bin[i].inst[0] = ainst[0];
      
#ifdef DEBUG
      calc_misfits_from_single_angle_set(kostrov->bin[i].s,angles,kostrov->bin[i].n, sdev);
      fprintf(stderr,"srandom:          %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tsdev: %7.5f %7.5f\tinst: %7.5f %7.5f\n",
	      kostrov->bin[i].s[0],kostrov->bin[i].s[1],kostrov->bin[i].s[2],kostrov->bin[i].s[3],kostrov->bin[i].s[4],kostrov->bin[i].s[5],
	      1-sdev[0],1-sdev[1],ainst[0],ainst[1]);
      optim_verbose = BC_TRUE;
#else
      optim_verbose = BC_FALSE;
#endif
      if(catalog->use_friction_solve){
	adjust_stress_for_friction(kostrov->bin[i].n,angles,weights,kostrov->bin[i].s,
				   kostrov->bin[i].def_s,kostrov->bin[i].best_s,&(kostrov->bin[i].best_fric),
				   (kostrov->bin[i].inst+1),(kostrov->bin[i].inst+2),optim_verbose,
				   ((catalog->use_friction_solve>1)?(BC_TRUE):(BC_FALSE)),&rsweep,&osweep);
	
	fprintf(stderr,"bin %05i: adjust_stress_for_friction: df sweep %i total sweep %i\n",i,rsweep,osweep);
      }
    }
  }
  free(angles);free(weights);
}


/* 
   solve for the five component stress tensor a la Andy Michael 

   angles [6 = 2 times 3] are given in radians TWO SETS OF ANGLES
   
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
				       long int *seed)
{
  const int npar = BC_MICHAEL_NPAR;
  const int ndim = BC_NDIM;
  const int nrandom_limit  = BC_MICHAEL_NMC; /* monte carlo simulations (2000
						good number for 95%
						confidence?)*/
  BC_BOOLEAN proceed;
  int nobs,iquake,nrandom,i,iquake6;
  BC_CPREC ind_stress[6],*slick,*amat,tot_stress[6],tot_stress2[6],snorm;

  slick = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim);
  amat = (BC_CPREC *)malloc(sizeof(BC_CPREC)*ndim*npar);
  proceed  = BC_TRUE;nrandom=0;
  for(i=0;i < 6;i++)
    tot_stress[i] = tot_stress2[i] = 0.0;
  do{
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
    for(i=0;i<6;i++){
      tot_stress[i]  += ind_stress[i];
      tot_stress2[i] += ind_stress[i] * ind_stress[i];
    }
    nrandom++;
    if(nrandom > nrandom_limit)
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
  /*
  fprintf(stderr,"stress_solve: avg: %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g) %12g(%12g)\n",
	  stress[0],sig_stress[0],stress[1],sig_stress[1],stress[2],sig_stress[2],
	  stress[3],sig_stress[3],stress[4],sig_stress[4],stress[5],sig_stress[5]);
  */
  
 
  free(amat);
  free(slick);
}

/* istress = iniitial stress

   dstress = stress for default friction
   bstress = stress for best friction, given deviation
   best_fric
 */
void adjust_stress_for_friction(int n, BC_CPREC *iangles, BC_CPREC *weights,
				BC_CPREC *istress,BC_CPREC *dstress,BC_CPREC *bstress,
				BC_CPREC *best_fric, BC_CPREC *dinst,BC_CPREC *binst,
				BC_BOOLEAN verbose, BC_BOOLEAN optimize,int *rsweep,int *osweep)
{
  BC_CPREC max[2],friction;
  BC_CPREC tstress[6],inst[2],sdev[2],mdev[2];
  BC_CPREC *angles;
  size_t asize,ssize;
  int i,isweep;
  asize = sizeof(BC_CPREC)*n*6;
  ssize = sizeof(BC_CPREC)*6;

  angles=(BC_CPREC *)malloc(asize);

  
  /* 
     default friction 
  */
  memcpy(dstress,istress,ssize);
  memcpy(angles,iangles,asize);
  optimize_angles_via_instability(n,angles,weights,BC_FRIC_DEF,dstress,inst,rsweep); /* for default friction */
  if(inst[1] > inst[0]){	/* messed up */
    for(i=0;i<n;i++)		/* swap all angles */
      swap_angles((angles+i*6));
    optimize_angles_via_instability(n,angles,weights,BC_FRIC_DEF,dstress,inst,&isweep);
    *rsweep += isweep;
  }
  /* compute the slip misgit angles */
  calc_misfits_from_single_angle_set(dstress,angles,n, sdev);
  if(verbose)
    fprintf(stderr,"default: s(%4.2f): %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tsdev: %7.5f %7.5f\tinst: %7.5f %7.5f\n",
	    BC_FRIC_DEF,dstress[0],dstress[1],dstress[2],dstress[3],dstress[4],dstress[5],1-sdev[0],1-sdev[1],inst[0],inst[1]);
  *dinst = inst[0];		/* keep default instability */

  *osweep =0;
  if(optimize){
    /* optimize */
    max[0] = -10;
    for(friction=BC_FRIC_SCAN_INC;friction < 1;friction += BC_FRIC_SCAN_INC){ /* scan
										 friction
										 for
										 maximum
										 instability */
      memcpy(tstress,istress,ssize);
      memcpy(angles,iangles,asize);
      optimize_angles_via_instability(n,angles,weights,friction,tstress,inst,&isweep);
      *osweep += isweep;
#ifdef DEBUG
      fprintf(stderr,"         s(%4.2f): %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tsdev:               \tinst: %7.5f %7.5f\n",
	      friction,tstress[0],tstress[1],tstress[2],tstress[3],tstress[4],tstress[5],inst[0],inst[1]);
#endif
      if(inst[1] > inst[0]){
	for(i=0;i<n;i++)		/* swap all angles */
	  swap_angles((angles+i*6));
	optimize_angles_via_instability(n,angles,weights,friction,tstress,inst,&isweep);
	*osweep += isweep;
      }
      if(inst[0] > max[0]){
	memcpy(bstress,tstress,ssize);
	calc_misfits_from_single_angle_set(bstress,angles,n, sdev);
	max[0] = inst[0];
	max[1] = inst[1];
	mdev[0] = sdev[0];
	mdev[1] = sdev[1];
	*best_fric = friction;
      }
    }
    *binst = max[0];
    if(verbose)
      fprintf(stderr,"optim:   s(%4.2f): %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tsdev: %7.5f %7.5f\tinst: %7.5f %7.5f\n",
	      *best_fric,bstress[0],bstress[1],bstress[2],bstress[3],bstress[4],bstress[5],1-mdev[0],1-mdev[1],max[0],max[1]);
  }else{
    *binst = NAN;
  }
  free(angles);
}







/* 

   JUST ONE SET OF ANGLES TO USE 

 */
void solve_stress_michael_specified_plane(int nquakes, BC_CPREC *angles,BC_CPREC *weights,BC_CPREC *stress)
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
  michael_solve_lsq(npar,ndim,nobs,amat,slick,weights,stress);
  /* normalize tensor */
  snorm = tensor6_norm(stress);
  for(i=0;i<6;i++)
    stress[i]     /= snorm;
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
  
  stress[BC_PP] =  lstress[0];	/*  xx */
  stress[BC_TP] = -lstress[1];	/* -yx */
  stress[BC_RP] =  lstress[2];	/*  xz */
  stress[BC_TT] =  lstress[3];	/*  yy */
  stress[BC_RT] = -lstress[4];	/* -yz */
  stress[BC_RR] =  lstress[5];	/*  zz */
  
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

   angles[3] strike, dip, rake angles (aki and richards) in radians

   from Andy Michael's slick.c routine

   NOTE: 
   
   this is converted from Michael's direction (i.e. strike + pi/2),
   dip, rake to Aki and Richards
   
*/
void michael_assign_to_matrix(BC_CPREC *angles,int *nobs,BC_CPREC **slick,BC_CPREC **amat)
{
  int j;
  BC_CPREC sin_z,cos_z,sin_z2,cos_z2,sin_z3,cos_z3,n1,n2,n3,n12,n22,n32;
  int ind;
  const int ndim = BC_NDIM;
  const int npar = BC_MICHAEL_NPAR;
  
  
  /*  */
  sincos(angles[0]+ M_PI_2, &sin_z, &cos_z); /* this needs to be flipped it seems */
  sincos(angles[1],         &sin_z2,&cos_z2);
  sincos(angles[2],         &sin_z3,&cos_z3);
    
    
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
  ind = j*npar;
  *(*amat+ind+0)    = n1-n12*n1+n1*n32;
  *(*amat+ind+1)    = n2-2.*n12*n2;
  *(*amat+ind+2)    = n3-2.*n12*n3;
  *(*amat+ind+3)    = -n1*n22+n1*n32;
  *(*amat+ind+4)    = -2.*n1*n2*n3;

  ind += npar;
  *(*amat+ind+0)= -n2*n12+n2*n32;
  *(*amat+ind+1)= n1-2.*n1*n22;
  *(*amat+ind+2)= -2.*n1*n2*n3;
  *(*amat+ind+3)= n2-n22*n2+n2*n32;
  *(*amat+ind+4)= n3-2.*n22*n3;

  ind += npar;
  *(*amat+ind+0)= -n3*n12-n3+n32*n3;
  *(*amat+ind+1)= -2.*n1*n2*n3;
  *(*amat+ind+2)= n1-2.*n1*n32;
  *(*amat+ind+3)= -n3*n22-n3+n32*n3;
  *(*amat+ind+4)= n2-2.*n2*n32;

  /* increment counter and make room for mroe */
  *nobs = *nobs + 1;
  /*  */
  *slick = (BC_CPREC *)realloc(*slick,sizeof(BC_CPREC)*ndim*     (*nobs+1));
  if(!slick)BC_MEMERROR("assign to Michael mat, slick");
  *amat = (BC_CPREC *) realloc(*amat, sizeof(BC_CPREC)*ndim*npar*(*nobs+1));
  if(!amat)BC_MEMERROR("assign to Michael mat, amat");
}
