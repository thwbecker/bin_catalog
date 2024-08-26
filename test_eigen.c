#include "catalog.h"
#include "eigen.h"


int main(int charc,char **argv)
{
  BC_CPREC a[3][3],eval[3],evec[9];
  int i;
  while(fscanf(stdin, "%lf %lf %lf %lf %lf %lf",
	       &a[0][0],&a[0][1],&a[0][2],&a[1][1],&a[1][2],&a[2][2])==6){
    /* eigenbalues are sorted bottom up */
    calc_eigensystem_sym_3x3(a,eval,evec,BC_TRUE,BC_TRUE);
    for(i=0;i<3;i++)
      fprintf(stdout,"%i\t%12g\t%12.6f %12.6f %12.6f\t%12.6f %12.6f %12.6f\t%11g\n",
	      i+1,eval[i],evec[i*3+0],evec[i*3+1],evec[i*3+2],
	      evec[i*3+0]/evec[i*3+2],evec[i*3+1]/evec[i*3+2],evec[i*3+2]/evec[i*3+2],
	      sqrt(evec[i*3+0]*evec[i*3+0]+evec[i*3+1]*evec[i*3+1]+evec[i*3+2]*evec[i*3+2]));

  }
}
