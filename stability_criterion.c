#include "catalog.h"
/*

these routines are modified from STRESSINVERSE
https://www.ig.cas.cz/en/stress-inverse/ 

Vavryčuk, V., 2014. Iterative joint inversion for stress and fault
orientations from focal mechanisms, Geophysical Journal International,
199, 69-77, doi: 10.1093/gji/ggu224.


the logic is that re here rely on a principal component system, so the
coordinate convention for the stress tensor should not matter

the routines get a set of angles (first strike/dip/rake , second
strike/dip/rake and return the degree of instability for each

then, routines will swap angles to make the first more unstable

*/
/* 
   with a given stress field, optimize via instability, swapping
   angles - this does not change the ordering in the quakes structure

   this gets a stress tensor and returns a new one, overwriting the old
*/
void optimize_angles_via_instability(int n,BC_CPREC *angles,
				     BC_CPREC *weights,BC_CPREC friction,
				     BC_CPREC *stress,BC_CPREC *ainst, int *nsweep)
{
  BC_CPREC svec[9],sigma[3],saved_stress[6];
  BC_CPREC inst[2],*bangles,inst_max[2];
  int j,j6,nswap,nmswap,max_sweep;
  BC_BOOLEAN *toswap;
  long int seed;
  size_t  asize,ssize;
  seed = -(long int)n;
  asize = sizeof(BC_CPREC)*n*6;
  ssize = 6*sizeof(BC_CPREC);
  /*  */
  bangles = (BC_CPREC *)malloc(asize);
  toswap = (BC_BOOLEAN *)malloc(sizeof(BC_BOOLEAN)*n);
  if((!bangles)||(!toswap))
    BC_MEMERROR("optimize_angles_via_instability");
  max_sweep = BC_MAX_SWEEP_FAC * n;
  if(max_sweep < 10)
    max_sweep = 10;
  if(max_sweep > 100)
    max_sweep = 100;


  nswap=1;
  nmswap=0;
  inst_max[0]=-10;
  while(nswap){
    *nsweep =0;
    while(nswap && (*nsweep < max_sweep)){
      /* solve stresses anew with current planes */
      solve_stress_michael_specified_plane(n,angles,weights,stress);
      /* compute the eigensystem for this set of stresses */
      calc_eigensystem_vec6(stress,sigma,svec,BC_TRUE,BC_TRUE);
      /* compute stability and take note of those planes to swap */
      nswap=0;
      ainst[0]=ainst[1]=0.0;
      for(j=j6=0;j < n;j++,j6+=6){
	/* 
	   check each pair for stability assume stress has changed
	   only for the first call, others are for different angles
	   but same stress
	*/
	stability_criterion_eig(sigma,svec,friction,(angles+j6),(j==0)?(BC_TRUE):(BC_FALSE),inst);
	if(inst[1] > inst[0]){	/* swap to make first more unstable */
	  toswap[j]=BC_TRUE;
	  nswap++;
	}else{
	  toswap[j]=BC_FALSE;
	}
	ainst[0]+=inst[0];
	ainst[1]+=inst[1];
      }
      /* average instability on first and second*/
      ainst[0]/=(BC_CPREC)n;
      ainst[1]/=(BC_CPREC)n;
      /* this particular set is better than all before */
      if(ainst[0] > inst_max[0]){
	inst_max[0] = ainst[0];
	inst_max[1] = ainst[1];
	memcpy(bangles,angles,asize);
	memcpy(saved_stress,stress,ssize);
      }
#ifdef SUPERDEBUG    
      fprintf(stderr,"f: %.2f s: %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tainsta: %12g %12g nswap %i nsweep %i\n",
	      friction,stress[0],stress[1],stress[2],stress[3],stress[4],stress[5],ainst[0],ainst[1],nswap,*nsweep);
#endif
      /* now acutally swap */
      for(j=j6=0;j < n;j++,j6+=6)
	if(toswap[j])
	  swap_angles((angles+j6));
      *nsweep += 1;
    }
    if(nswap){			/* we have reached the max number of iterations */
      if(nmswap < 2){
#ifdef SUPERDEBUG
	fprintf(stderr,"optimize_angles_via_instability: still swapping at %i out of %i, major loop %i (max inst: %g), randomizing\n",
		*nsweep,max_sweep,nmswap,inst_max[0]);
#endif
	for(j=j6=0;j < n;j++,j6+=6){
	  if(BC_RGEN(&seed) >= 0.5){
	    swap_angles((angles+j6));
	  }
	}
      }else{			/* bailing with best */
#ifdef SUPERDEBUG
	fprintf(stderr,"optimize_angles_via_instability: still swapping at %i out of %i, major loop %i (max inst: %g), bailing with best\n",
		*nsweep,max_sweep,nmswap,inst_max[0]);
#endif
	memcpy(angles,bangles,asize);
	memcpy(stress,saved_stress,6*sizeof(BC_CPREC));
	ainst[0] = inst_max[0];
	ainst[1] = inst_max[1];
	nswap=0;
      }
    }
    nmswap++;
  }
  free(toswap);
  free(bangles);
}
/* compute average instability of set */
void calc_average_instability(int n,BC_CPREC *angles,BC_CPREC *weights,BC_CPREC friction, BC_CPREC *stress,BC_CPREC *ainst)
{
  int i,i6;
  BC_CPREC sigma[3],svec[9],inst[2];

  calc_eigensystem_vec6(stress,sigma,svec,BC_TRUE,BC_TRUE);
  ainst[0]=ainst[1]=0.0;
  for(i=i6=0;i < n;i++,i6+=6){
    stability_criterion_eig(sigma,svec,friction,(angles+i6),(i==0)?(BC_TRUE):(BC_FALSE),inst);
    ainst[0]+=inst[0];
    ainst[1]+=inst[1];
  }
  /* average instability on first and second*/
  ainst[0]/=(BC_CPREC)n;
  ainst[1]/=(BC_CPREC)n;
}

/*

modified from STRESSINVERSE 
https://www.ig.cas.cz/en/stress-inverse/
Vavryčuk, V., 2014. Iterative joint inversion for stress and fault orientations from focal mechanisms, Geophysical Journal International, 199, 69-77, doi: 10.1093/gji/ggu224.

%*************************************************************************%
%                                                                         %
%  function STABILITY_CRITERION                                           %
%                                                                         %
%  function calculates the fault instability and and identifies faults    %
%  with unstable nodal planes                                             %
%                                                                         %
%  input:  stress                                                         %
%          friction                                                       %
%          complementary focal mechanisms                                 %
%                                                                         %
%  output: focal mechanisms with correct fault orientations               %
%          instability of faults                                          %
%                                                                         %
%*************************************************************************%
*/


/* given eigensystem, ordered largest to smallest, compute stability
   on both fault planes 

   takes two sets of angles
*/
void  stability_criterion_eig(BC_CPREC *sigma, BC_CPREC *svec,
			      BC_CPREC friction,BC_CPREC *angles,
			      BC_BOOLEAN stress_has_changed,
			      BC_CPREC *instability)
{
  BC_CPREC n1[3],n2[3],n1_[3],n2_[3];
  BC_CPREC sin_strike1,cos_strike1,sin_strike2,cos_strike2;
  BC_CPREC sin_dip1,cos_dip1,sin_dip2,cos_dip2;

  BC_CPREC n1_0_p2, n1_1_p2, n1_2_p2;
  BC_CPREC n2_0_p2, n2_1_p2, n2_2_p2;
  BC_CPREC tmp_pow,tau_shear_n1_norm,tau_normal_n1_norm;
  BC_CPREC tau_shear_n2_norm,tau_normal_n2_norm;
  static BC_CPREC shape_ratio,fric_fac,shape_fac,shape_facp2;
  static BC_CPREC sigma_vector_1[3],sigma_vector_2[3],sigma_vector_3[3];
  const size_t nsize3=3*sizeof(BC_CPREC);
  if(stress_has_changed){
    /* recompute stress dependent factors */
    fric_fac = friction + sqrt(1.0+friction*friction); 
    shape_ratio = (sigma[0]-sigma[1])/(sigma[0]-sigma[2]);
    
    shape_fac = (1.0-2.0*shape_ratio);
    shape_facp2 = shape_fac * shape_fac;
    memcpy(sigma_vector_1,(svec),  nsize3);
    memcpy(sigma_vector_2,(svec+3),nsize3);
    memcpy(sigma_vector_3,(svec+6),nsize3);
  }
  //--------------------------------------------------------------------------
  //  two alternative fault normals
  //--------------------------------------------------------------------------
  // first fault normal
  sincos(angles[0],&sin_strike1,&cos_strike1);
  sincos(angles[1],&sin_dip1,&cos_dip1);
  /* second */
  sincos(angles[3],&sin_strike2,&cos_strike2);
  sincos(angles[4],&sin_dip2,&cos_dip2);

  n1[0] = -sin_dip1 * sin_strike1;
  n1[1] =  sin_dip1 * cos_strike1;
  n1[2] = -cos_dip1;
  
  // second fault normal
  n2[0] = -sin_dip2 * sin_strike2;
  n2[1] =  sin_dip2 * cos_strike2;
  n2[2] = -cos_dip2;
  
  //--------------------------------------------------------------------------
  // notation: sigma1 = 1; sigma2 = 1-2*shape_ratio; sigma3 = -1
  //--------------------------------------------------------------------------
  // fault plane normals in the coordinate system of the principal stress axes
  n1_[0] = n1[0]*sigma_vector_1[0] + n1[1]*sigma_vector_1[1] + n1[2]*sigma_vector_1[2];
  n1_[1] = n1[0]*sigma_vector_2[0] + n1[1]*sigma_vector_2[1] + n1[2]*sigma_vector_2[2];
  n1_[2] = n1[0]*sigma_vector_3[0] + n1[1]*sigma_vector_3[1] + n1[2]*sigma_vector_3[2];
  
  n2_[0] = n2[0]*sigma_vector_1[0] + n2[1]*sigma_vector_1[1] + n2[2]*sigma_vector_1[2];
  n2_[1] = n2[0]*sigma_vector_2[0] + n2[1]*sigma_vector_2[1] + n2[2]*sigma_vector_2[2];
  n2_[2] = n2[0]*sigma_vector_3[0] + n2[1]*sigma_vector_3[1] + n2[2]*sigma_vector_3[2];

  n1_0_p2 = n1_[0]*n1_[0];
  n1_1_p2 = n1_[1]*n1_[1];
  n1_2_p2 = n1_[2]*n1_[2];
  
  n2_0_p2 = n2_[0]*n2_[0];
  n2_1_p2 = n2_[1]*n2_[1];
  n2_2_p2 = n2_[2]*n2_[2];
    
  //--------------------------------------------------------------------------
  // 1. alternative
  //--------------------------------------------------------------------------
  tmp_pow = n1_0_p2 + shape_fac*n1_1_p2 -n1_2_p2;
  tmp_pow *= tmp_pow;
  
  tau_shear_n1_norm   = sqrt(n1_0_p2 + shape_facp2*n1_1_p2 + n1_2_p2 - tmp_pow);
  tau_normal_n1_norm = n1_0_p2 + shape_fac*n1_1_p2 - n1_2_p2;
  
  //--------------------------------------------------------------------------
  // 2. alternative
  //--------------------------------------------------------------------------
  tmp_pow = n2_0_p2 + shape_fac*n2_1_p2 - n2_2_p2;
  tmp_pow *= tmp_pow;
  
  tau_shear_n2_norm   = sqrt(n2_0_p2 + shape_facp2*n2_1_p2 + n2_2_p2 - tmp_pow);
  tau_normal_n2_norm = n2_0_p2 + shape_fac*n2_1_p2 - n2_2_p2;
  
  //--------------------------------------------------------------------------
  // instability
  //--------------------------------------------------------------------------
  instability[0] = (tau_shear_n1_norm - friction*(tau_normal_n1_norm-1))/fric_fac;
  instability[1] = (tau_shear_n2_norm - friction*(tau_normal_n2_norm-1))/fric_fac;
  
  
  
}

