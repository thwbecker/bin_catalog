#include "catalog.h"
/* 


merge two aki and richards formatted earthquake catalogs, trying to remove duplicates


*/

int main(int argc, char **argv)
{
  int i;
  struct cat *catalog;

  
  if(argc != 4){
    fprintf(stderr,"%s catalog.1.aki catalog.2.aki catalog.merged.aki for doubles, first catalog wins\n",argv[0]);
    exit(-1);
  }
  fprintf(stderr,"%s: merging two catalogs into %s\n",argv[0],argv[3]);
  
  catalog=(struct cat *)malloc(3  * sizeof(struct cat));
  if(!catalog)
    BC_MEMERROR(argv[0]);
  
  for(i=0;i < 2;i++)
    read_catalog(argv[i+1],(catalog+i),BC_AKI);
  /* merge catalog */
  create_catalog((catalog+2),-1);
  /* 1 km slop in horizontal, 15 km in vertical, largest event wins */
  merge_catalog((catalog),(catalog+1),(catalog+2),1,15);

  /* print  */
  print_catalog(argv[3],(catalog+2),BC_AKI);

  return 0;
}

