#include "catalog.h"
/* 


read GMT psmeca moment tensor format and print as best fit double
couple fault plane


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
  if(!catalog)
    BC_MEMERROR(argv[0]);
  read_catalog(argv[1],catalog,BC_CMT);
  for(i=0;i<catalog->n;i++){
    print_quake_cmt_fp(stdout,catalog->quake[i]);
  }
}
