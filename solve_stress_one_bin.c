#include "catalog.h"
/* 
   read 
   - MY AKI style events, or
   - strike, dip, rake format 
   
   and perform stress inversion, Michael or Vavrychuk type

   (c) 2020 - 2026, Thorsten Becker, thwbecker@post.harvard.edu, see
   README
*/
int main(int argc, char **argv)
{
 struct cat *catalog;
 BC_CPREC *weights,*angles,sdev[2],use_friction,fopt,shape_ratio,
   sig_stress[6],stress[6],nstress[6],dstress[6],bstress[6],best_fric,mean_instability;
 BC_SWITCH *select;
 int i,i6,rsweep,tsweep,optimize;
 long int seed = -1;
 struct qke quake[1];
 catalog=(struct cat *)calloc(1,sizeof(struct cat));
 if(argc < 2){
   fprintf(stderr,"%s catalog.aki\n",argv[0]);
   fprintf(stderr,"\tif the catalog filename is stdin_sdr will read strike dip rake in degree from stdin instead\n"); 
   exit(-1);
 }
 if(strcmp(argv[1], "stdin_sdr") == 0){
   catalog->n = i6=0;
   angles = (BC_CPREC *)malloc(6*sizeof(BC_CPREC));
   weights = (BC_CPREC *)malloc(sizeof(BC_CPREC));
   select = (BC_SWITCH *)malloc(sizeof(BC_SWITCH));
   while(fscanf(stdin,"%lf %lf %lf",
		&(quake->strike),&(quake->dip),&(quake->rake))==3){
     quake->strike = BC_D2R(quake->strike); /* strike/dip/rake in
					       radians internally */
     quake->dip    = BC_D2R(quake->dip);
     quake->rake   = BC_D2R(quake->rake);
     find_alt_plane(quake->strike, quake->dip, quake->rake,
		    &(quake->strike2), &(quake->dip2), &(quake->rake2));
     weights[catalog->n]=1.0;
     select[catalog->n]=0;
     assign_quake_angles(quake,(angles+i6));
     
     catalog->n++;
     i6+=6;
     angles = (BC_CPREC *)realloc(angles,6*sizeof(BC_CPREC)*(catalog->n+1));
     weights = (BC_CPREC *)realloc(weights,sizeof(BC_CPREC)*(catalog->n+1));
     select = (BC_SWITCH *)realloc(select,sizeof(BC_SWITCH)*(catalog->n+1));
 
   }
   
 }else{
   read_catalog(argv[1],catalog,BC_AKI,BC_FALSE);
   angles = (BC_CPREC *)malloc(6*sizeof(BC_CPREC)*catalog->n);
   weights = (BC_CPREC *)malloc(sizeof(BC_CPREC)*catalog->n);
   select = (BC_SWITCH *)malloc(sizeof(BC_SWITCH)*catalog->n);
   for(i=i6=0;i<catalog->n;i++,i6+=6){
     weights[i]=1.0;
     assign_quake_angles((catalog->quake+i),(angles+i6));
     select[i]=0;
   }
 }
 /* solve Michael style for those particular planes */
 solve_stress_michael_specified_plane(catalog->n, angles, weights,stress,BC_FALSE);
 fprintf(stderr,"single, no norm:  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\ttr: %7.4f\n",
	 stress[BC_RR],stress[BC_RT],stress[BC_RP],stress[BC_TT],stress[BC_TP],stress[BC_PP],trace6(stress));

 /* max abs eigen value normalized */
 max_ev_normalize_tens6(stress,nstress);
 fprintf(stderr,"evmax_norm:       %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	 nstress[BC_RR],nstress[BC_RT],nstress[BC_RP],nstress[BC_TT],nstress[BC_TP],nstress[BC_PP]);

 /* proper normalized */
 normalize_tens6(stress,nstress);
 fprintf(stderr,"s norm:           %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	 nstress[BC_RR],nstress[BC_RT],nstress[BC_RP],nstress[BC_TT],nstress[BC_TP],nstress[BC_PP]);

 /* needed for instability and varychuk inversion */
 use_friction = 0.6;
  
 /* 
    randomized Michael
 */
 solve_stress_michael_random_sweep(catalog->n, angles,weights,stress, sig_stress,&seed,BC_MICHAEL_RSWEEP_MAX);
 calc_misfits_from_single_angle_set(stress,angles,catalog->n, sdev);
 fprintf(stderr,"srandom:          %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\tsdev: %7.5f %7.5f\n\n",
	 stress[0],stress[1],stress[2],stress[3],stress[4],stress[5],1-sdev[0],1-sdev[1]);
 /* 
    Varychuk for default and best friction 
 */
 /* 
    default 
 */
 stress_inversion_mstyle(catalog->n, angles, weights, use_friction, use_friction, 0.1,
			 BC_VI_ITER,  BC_VI_NREAL, &seed,
			 dstress, &shape_ratio, &best_fric, &mean_instability, NULL,
			 BC_STRESS_NORM_TENSOR);
 fprintf(stderr, "Vav fixed R = %.4f   friction_opt = %.3f   mean_instability = %.5f\n",shape_ratio, best_fric, mean_instability);
 fprintf(stderr,"s norm:           %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	 dstress[BC_RR],dstress[BC_RT],dstress[BC_RP],dstress[BC_TT],dstress[BC_TP],dstress[BC_PP]);
 
 max_ev_normalize_tens6(dstress,nstress); 
 fprintf(stderr,"evmax_norm:       %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\ttr: %7.4f\n",
	 nstress[BC_RR],nstress[BC_RT],nstress[BC_RP],nstress[BC_TT],nstress[BC_TP],nstress[BC_PP],trace6(nstress));
 /* 
    sweep
 */
 stress_inversion_mstyle(catalog->n, angles, weights, 0.01, 0.99, 0.001,BC_VI_ITER,  BC_VI_NREAL, &seed,
			 dstress, &shape_ratio, &best_fric, &mean_instability, NULL,
			 BC_STRESS_NORM_TENSOR);
 fprintf(stderr, "Vav sweep R = %.4f   friction_opt = %.3f   mean_instability = %.5f\n",shape_ratio, best_fric, mean_instability);
 fprintf(stderr,"s norm:           %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	 dstress[BC_RR],dstress[BC_RT],dstress[BC_RP],dstress[BC_TT],dstress[BC_TP],dstress[BC_PP]);
 



 return 0;
}
