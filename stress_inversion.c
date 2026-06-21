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
  BC_CPREC *angles,*weights,adotp[2],shape_ratio,fopt,mdotp,fmin,fmax,finc;
  BC_CPREC fr_best, fr_mean, fr_std, fr_16, fr_84;
  BC_BOOLEAN calc_fric_uncertainty;
#ifdef DEBUG
  int k;
  BC_CPREC dotp[2];
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
      /* instability of the Michael stress, MATLAB/Vavrycuk convention  */
      vavrycuk_average_instability(kostrov->bin[i].n,angles,weights,BC_FRIC_DEF,kostrov->bin[i].s,adotp);
      kostrov->bin[i].dotp[0] = adotp[0];
      
#ifdef DEBUG
      calc_misfits_from_single_angle_set(kostrov->bin[i].s,angles,kostrov->bin[i].n, dotp);
      fprintf(stderr,"srandom:          %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tdotp: %7.5f a%7.5f\tdotp: %7.5f %7.5f\n",
	      kostrov->bin[i].s[0],kostrov->bin[i].s[1],kostrov->bin[i].s[2],kostrov->bin[i].s[3],
	      kostrov->bin[i].s[4],kostrov->bin[i].s[5],
	      dotp[0],dotp[1],adotp[0],adotp[1]);
#endif
      if(catalog->use_friction_solve){
	/* Vavrycuk , default fixed friction -> def_s, dotp[1] */
	stress_inversion_vavrycuk(kostrov->bin[i].n,angles,weights,BC_FRIC_DEF,BC_FRIC_DEF,1.0,
				BC_VI_ITER,BC_VI_NREAL,&catalog->seed,
				kostrov->bin[i].def_s,&shape_ratio,&fopt,&mdotp,NULL,
				BC_STRESS_NORM_TENSOR);
	kostrov->bin[i].dotp[1] = mdotp;
	/* best friction by scan -> best_s, best_fric, dotp[2].
	   scan range mirrors the old optimize flag (use_friction_solve-1):
	   1 -> 0..1, else -> 0.2..0.8 */
	opt = catalog->use_friction_solve - 1;
	if(opt == 1){		/* use_friction == 2 */
	  fmin = 0.0; fmax = 1.0; finc = BC_FRIC_SCAN_INC;calc_fric_uncertainty = BC_FALSE;
	}else if(opt == 2){	/* use_friction == 3 */
	  fmin = 0.2; fmax = 0.8; finc = BC_FRIC_SCAN_INC2;;calc_fric_uncertainty = BC_FALSE;
	}else if(opt == 3){	/* use_friction == 3 */
	  fmin = 0.0; fmax = 1.0; finc = BC_FRIC_SCAN_INC;calc_fric_uncertainty = BC_TRUE;
	}
	stress_inversion_vavrycuk(kostrov->bin[i].n,angles,weights,fmin,fmax,finc,
				  BC_VI_ITER,BC_VI_NREAL,&catalog->seed,
				  kostrov->bin[i].best_s,&shape_ratio,&(kostrov->bin[i].best_fric),
				  &mdotp,NULL,BC_STRESS_NORM_TENSOR);
	if(calc_fric_uncertainty){
	  vavrycuk_friction_error(kostrov->bin[i].n, angles, weights,fmin,fmax,finc,
				  BC_VI_ITER, BC_VI_NREAL, BC_VI_MC_ERR, &catalog->seed,
				  &fr_best, &fr_mean, &fr_std, &fr_16, &fr_84);
	  kostrov->bin[i].best_fric = fr_mean; /* replace with mean */
	  kostrov->bin[i].std_fric = fr_std; 
	}else{
	  kostrov->bin[i].std_fric = NAN; 
	}
	kostrov->bin[i].dotp[2] = mdotp;
      }
    }
  }
  free(angles);free(weights);
}

