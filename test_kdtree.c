#include "catalog.h"

/* 
   tester for 
   KD tree type implementation for geographic data, partially based on
   code from Claude

 */

// Demo function
void demo() {
  geo_tree_t *tree;
  FILE *in;
  BC_CPREC dlon,dlat;
  unsigned int code;
  result_array_t *radius_results;
  printf("Geographic Tree Search Demo\n");
  printf("===========================\n");
  
  // Create tree
  tree = geo_tree_create(2000000);
  if (!tree) {
    printf("Failed to create tree\n");
    return;
  }
  in = fopen("cheng25.aki","r");
  code=0;
  while(fscanf(in,"%lf %lf %*f %*f %*f %*f %*f %*f %*f %*s",&dlon,&dlat)==2 ){
    geo_tree_add_point(tree, BC_D2R(dlat),BC_D2R(dlon),code);
    code++;
  }
  fclose(in);
  printf("Building tree with %d events...\n", tree->num_points);
  clock_t start = clock();
  geo_tree_build(tree);
  clock_t end = clock();
  BC_CPREC build_time = ((BC_CPREC)(end - start)) / CLOCKS_PER_SEC;
  printf("Tree built in %.3f seconds\n", build_time);
  
  // Test 1: Radius query
  dlon = -119;dlat = 36;
  printf("\n1. Events within 20 km of (%g, %g):",dlat,dlon);
  radius_results = geo_tree_query_radius(tree,  BC_D2R(dlat),  cos(BC_D2R(dlat)),BC_D2R(dlon), 20, 0);
  print_results(radius_results, "Events ");
  result_array_destroy(radius_results);
  
  // Test 2: K-nearest neighbors
  //printf("\n2. 3 nearest events to San Francisco (37.7749, -122.4194):");
  //result_array_t* knn_results = geo_tree_query_k_nearest(tree, 37.7749, -122.4194, 3);
  //print_results(knn_results, "3 nearest events to San Francisco");
  //result_array_destroy(knn_results);
    
 
  
  

  geo_tree_destroy(tree);
}

int main() {
    demo();
    return 0;
}
