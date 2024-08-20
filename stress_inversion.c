#include "catalog.h"
/* 

   bunch of stress inversion stuff


 */

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
