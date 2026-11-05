#include "catalog.h"
/* 


read GMT psmeca moment format with time at end and print magnitude


*/

int main(int argc, char **argv)
{
  int i;
  struct cat *catalog;
  if(argc<2){
    fprintf(stderr,"%s catalog.cmt\n",argv[0]);
    exit(-1);
  }
  catalog=(struct cat *)calloc(1,sizeof(struct cat)); 
  if(!catalog)BC_MEMERROR(argv[0]);
  read_catalog(argv[1],catalog,BC_CMT,BC_FALSE);
  for(i=0;i<catalog->n;i++){
    fprintf(stdout,"%10.4f %10.4f %5.1f %12.5e %6.3f %012i\n",
	    catalog->quake[i].dlon,
	    catalog->quake[i].dlat,
	    catalog->quake[i].depth,
	    catalog->quake[i].m0,
	    catalog->quake[i].mag,
	    (int)catalog->quake[i].tsec);
  }
}
