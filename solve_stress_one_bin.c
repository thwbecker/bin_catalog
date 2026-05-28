#include "catalog.h"
/* 
   read AKI style events or in strike, dip, rake format and perform
   stress inversion

   (c) 2020 - 2026, Thorsten Becker, thwbecker@post.harvard.edu, see README
*/
int main(int argc, char **argv)
{
 struct cat *catalog;
 BC_CPREC *weights,*angles,sdev[2],
   sig_stress[6],stress[6],dstress[6],bstress[6],best_fric,bdev,ddev,inst[2];
 BC_SWITCH *select;
 int i,i6,rsweep,tsweep;
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
 /* solve for those particular planes */
 solve_stress_michael_specified_plane(catalog->n, angles, weights,stress);
 fprintf(stderr,"single1:          %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
	 stress[0],stress[1],stress[2],stress[3],stress[4],stress[5]);
 /* randomized  */
 solve_stress_michael_random_sweep(catalog->n, angles,weights,stress, sig_stress,&seed,BC_MICHAEL_RSWEEP_MAX);
 calc_misfits_from_single_angle_set(stress,angles,catalog->n, sdev);
 calc_average_instability(catalog->n,angles,weights,BC_FRIC_DEF, stress,inst);
 fprintf(stderr,"srandom:          %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\tsdev: %7.5f %7.5f\tinst: %7.5f %7.5f\n",
	 stress[0],stress[1],stress[2],stress[3],stress[4],stress[5],1-sdev[0],1-sdev[1],inst[0],inst[1]);
 /* for default and best friction */
 adjust_stress_for_friction(catalog->n,angles,weights,stress,dstress,bstress,&best_fric,&ddev,&bdev,
			    BC_TRUE,BC_TRUE,&rsweep,&tsweep);

 return 0;
}
