#include "catalog.h"
/* 


merge two aki and richards formatted earthquake catalogs, trying to remove duplicates


*/

int main(int argc, char **argv)
{
  int i;
  struct cat *catalog;

  
  if(argc != 4){
    fprintf(stderr,"%s catalog.1.aki catalog.2.aki catalog.merged.aki\n",argv[0]);
    fprintf(stderr,"\tmerges two AKI catalogs; doubles are entries within %g km horizontally, %g km in depth,\n",
	    BC_MERGE_HDIST_CLOSE,BC_MERGE_DEP_CLOSE);
    fprintf(stderr,"\tand %g s in origin time, of which the larger magnitude is kept\n",BC_DUP_TIME_CLOSE);
    exit(-1);
  }
  fprintf(stderr,"%s: merging two catalogs into %s\n",argv[0],argv[3]);
  
  catalog=(struct cat *)calloc(3,sizeof(struct cat));
  if(!catalog)
    BC_MEMERROR(argv[0]);
  
  for(i=0;i < 2;i++)
    read_catalog(argv[i+1],(catalog+i),BC_AKI,BC_FALSE);
  /* merge catalog */
  create_catalog((catalog+2),-1);
  /* thresholds from catalog.h, larger magnitude wins */
  merge_catalog((catalog),(catalog+1),(catalog+2),BC_MERGE_HDIST_CLOSE,BC_MERGE_DEP_CLOSE);

  /* print  */
  print_catalog(argv[3],(catalog+2),BC_AKI);

  return 0;
}

