#include "catalog.h"
/* 
   read strike, dip, rake and compute auxiliary plane

   (c) 2020 - 2026, Thorsten Becker, thwbecker@post.harvard.edu, see README
*/

int main(int argc, char **argv)
{
  struct qke quake[1];
  int n=0;
  fprintf(stderr,"%s: readin strike dip rake from stdin...\n",argv[0]);
  while(fscanf(stdin,"%lf %lf %lf",
	       &(quake->strike),&(quake->dip),&(quake->rake))==3){
    quake->strike = BC_D2R(quake->strike); /* strike/dip/rake in
					      radians internally */
    quake->dip    = BC_D2R(quake->dip);
    quake->rake   = BC_D2R(quake->rake);
    find_alt_plane(quake->strike, quake->dip, quake->rake,
		   &(quake->strike2), &(quake->dip2), &(quake->rake2));
    
    fprintf(stdout,"%9.3f %9.3f %9.3f\t%9.3f %9.3f %9.3f\n",
	    BC_R2D(quake->strike),  BC_R2D(quake->dip),  BC_R2D(quake->rake),
	    BC_R2D(quake->strike2),  BC_R2D(quake->dip2),  BC_R2D(quake->rake2));
    n++;
  }
  fprintf(stderr,"%s: computed %i aux planes\n",argv[0],n);
  return 0;
}

