#include "catalog.h"

/* 
   tester for 
   KD tree type implementation for geographic data, partially based on
   code from Claude

 */

// Demo function
void demo() {
  geo_tree_t *tree;
  FILE *in,*out;
  BC_CPREC dlon,dlat,dist;
  unsigned int code,mode;
  int ndist;
  char fname[300];
  result_array_t *results;
  fprintf(stderr,"Geographic Tree Search Demo\n");
  fprintf(stderr,"===========================\n");
  
  // Create tree
  tree = geo_tree_create(2000000);
  if (!tree) {
    printf("Failed to create tree\n");
    return;
  }
  in = fopen("cheng25.aki","r");
  code=0;
  while(fscanf(in,"%lf %lf %*f %*f %*f %*f %*f %*f %*f %*s",&dlon,&dlat)==2 ){
    geo_tree_add_point(tree, dlat,dlat,code);
    code++;
  }
  fclose(in);
  fprintf(stderr,"Building tree with %d events...\n", tree->num_points);
  clock_t start = clock();
  geo_tree_build(tree);
  clock_t end = clock();
  BC_CPREC build_time = ((BC_CPREC)(end - start)) / CLOCKS_PER_SEC;
  fprintf(stderr,"Tree built in %.3f seconds\n", build_time);
  
  // Test 0: Radius query
  // test 1: neighbor query
  //dlon = -119;dlat = 36;
  mode=0;
  
  dist = 25;ndist = 10;
  //for(dlon=-126;dlon<=-112;dlon+=2)
  //for(dlat=30;dlat<=42;dlat+=2){
  for(dlon=232;dlon<=250;dlon+=0.5)
    for(dlat=30;dlat<=45;dlat+=0.5){
      if(mode==0)
	results = geo_tree_query_radius(tree,  dlat,  dlon, dist, 0);
      else
	results = geo_tree_query_k_nearest(tree, dlat,dlon, ndist);
      if(results->count)
	fprintf(stdout,"%g %g %i\n",dlon,dlat,results->count);
      //if(results->count){
      //snprintf(fname,sizeof(fname),"p.%g.%g.dat",dlon,dlat);
      //out=fopen(fname,"w");print_results_file(results,out);fclose(out);
      //printf("printed %i events to %s\n",results->count,fname);
      //}
      result_array_destroy(results);
    }
  
  
  
 
  
  

  geo_tree_destroy(tree);
}

int main() {
    demo();
    return 0;
}
