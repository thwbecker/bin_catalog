#include "catalog.h"


/* 
   
   KD tree type implementation for geographic data, partially based on
   code from Claude code AI

 */
// Main implementation
/* 
   
   get initial tree structure, will return NULL if memory error

 */
geo_tree_t *geo_tree_create(int initial_capacity) {
  geo_tree_t *tree;
  tree = (geo_tree_t*)malloc(sizeof(geo_tree_t));
  if (!tree) {
    return NULL;
  }
  tree->points = (kd_tree_point *)malloc(sizeof(kd_tree_point) * initial_capacity);
  if (!tree->points) {
    free(tree);
    return NULL;
  }
  
  // Allocate nodes array (approximately same size as points for balanced tree)
  tree->nodes = (kd_node_t *)malloc(sizeof(kd_node_t) * initial_capacity);
  if (!tree->nodes) {
    free(tree->points);
    free(tree);
    return NULL;
  }
    
  tree->num_points = 0;
  tree->capacity = initial_capacity;
  tree->num_nodes = 0;
  tree->node_capacity = initial_capacity;
  tree->root_idx = -1;
  
  return tree;
}

void geo_tree_destroy(geo_tree_t *tree) {
  if (!tree) return;
  
  if (tree->points) {
    free(tree->points);
  }
  
  if (tree->nodes) {
    free(tree->nodes);
  }
  
  free(tree);
}

int allocate_node(geo_tree_t *tree) {
  int idx,new_capacity;
  kd_node_t *new_nodes;
  if (!tree) return -1;
  
  if (tree->num_nodes >= tree->node_capacity) {
    // Need to expand nodes array
    new_capacity = tree->node_capacity * 2;
    new_nodes = (kd_node_t *)realloc(tree->nodes,
				    sizeof(kd_node_t) * new_capacity);
    if (!new_nodes) return -1;
    
    tree->nodes = new_nodes;
    tree->node_capacity = new_capacity;
  }
  
  idx = tree->num_nodes;
  tree->num_nodes++;
  
  // Initialize node
  tree->nodes[idx].point_idx = -1;
  tree->nodes[idx].left_idx = -1;
  tree->nodes[idx].right_idx = -1;
  tree->nodes[idx].split_dim = 0;
  
  return idx;
}
/* create point, coordinates in radian */
int geo_tree_add_point(geo_tree_t *tree, BC_CPREC lat, BC_CPREC lon, int code) {
  int new_capacity;
  kd_tree_point *new_points;
  kd_node_t *new_nodes;
  kd_tree_point *point;
  if (!tree) return 0;
  
  // Resize points array if needed
  if (tree->num_points >= tree->capacity) {
    /* increase capacity by some factor to save time 
       this was 2 x before
    */
    new_capacity = (int)((float)tree->capacity * 1.25);
    /*  */
    new_points = (kd_tree_point *)realloc(tree->points, 
					sizeof(kd_tree_point) * new_capacity);
    if (!new_points) return 0;
    
    tree->points = new_points;
    tree->capacity = new_capacity;
    
    // Also resize nodes array
    new_nodes = (kd_node_t *)realloc(tree->nodes,
				     sizeof(kd_node_t) * new_capacity);
    if (!new_nodes) return 0;
    
    tree->nodes = new_nodes;
    tree->node_capacity = new_capacity;
  }
  
  // Add the point
  point = &tree->points[tree->num_points];
  point->lat = lat;
  point->cos_lat = cos((double)lat);
  point->lon = lon;
  point->code = code;
  
  tree->num_points++;
  return 1;
}

void geo_tree_build(geo_tree_t *tree) {
  int *indices;
  if (!tree || tree->num_points == 0) return;
    
  // Reset tree
  tree->num_nodes = 0;
  tree->root_idx = -1;
  
  // Create array of indices
  indices = (int*)malloc(sizeof(int) * tree->num_points);
  if (!indices) return;
  
  for (int i = 0; i < tree->num_points; i++) {
    indices[i] = i;
  }
  
  // Build the tree
  tree->root_idx = build_tree_recursive(tree, indices, tree->num_points, 0);
  
  free(indices);
}

int build_tree_recursive(geo_tree_t *tree, int *indices, int count, int depth) {
  sort_item_t *items;
  int dim, node_idx, median,i;
  kd_node_t *node;

  if (count <= 0) return -1;
    
  // Sort by current dimension (0 = lat, 1 = lon)
  dim = depth % 2;
  
  
  items = (sort_item_t *)malloc(sizeof(sort_item_t) * count);
  if (!items) return -1;
  
  for (i = 0; i < count; i++) {
    items[i].idx = indices[i];
    items[i].value = (dim == 0) ? tree->points[indices[i]].lat : tree->points[indices[i]].lon;
  }
  
  // Sort by the current dimension
  qsort(items, count, sizeof(sort_item_t), compare_by_val);
  
  // Update indices array with sorted order
  for (i = 0; i < count; i++) {
    indices[i] = items[i].idx;
  }
  free(items);
  
  // Choose median
  median = count / 2;
  
  // Allocate new node
  node_idx = allocate_node(tree);
  if (node_idx == -1) return -1;
  
  node = &tree->nodes[node_idx];
  node->point_idx = indices[median];
  node->split_dim = dim;
  node->left_idx = build_tree_recursive(tree, indices, median, depth + 1);
  node->right_idx = build_tree_recursive(tree, indices + median + 1, count - median - 1, depth + 1);
  
  return node_idx;
}

result_array_t *geo_tree_query_radius(geo_tree_t *tree, BC_CPREC center_lat, BC_CPREC cos_center_lat,
				      BC_CPREC center_lon,BC_CPREC radius_km, int max_results) {
  result_array_t *results;
  if (!tree || tree->root_idx == -1) return NULL;
    
  results = result_array_create(100);
  if (!results) return NULL;
  
  search_radius_recursive(tree, tree->root_idx, center_lat, cos_center_lat, center_lon, 
			  radius_km, results, max_results);
  
  // Sort results by distance
  qsort(results->results, results->count, sizeof(query_result_t), compare_by_distance);
  
  // Truncate to max_results if specified
  if (max_results > 0 && results->count > max_results) {
    results->count = max_results;
  }
  
  return results;
}

void search_radius_recursive(geo_tree_t *tree, int node_idx, BC_CPREC center_lat, BC_CPREC cos_center_lat, BC_CPREC center_lon,
			     BC_CPREC radius_km, result_array_t *results, int max_results) {
  kd_node_t *node;
  kd_tree_point *point;
  BC_CPREC distance;
  BC_CPREC diff, coord_diff;
  int dim;
  if (node_idx == -1) return;
  if (max_results > 0 && results->count >= max_results) return;
  
  node = &tree->nodes[node_idx];
  point = &tree->points[node->point_idx];

  distance = distance_geo(center_lon, center_lat, point->lon, point->lat, cos_center_lat, point->cos_lat);
  
  
  // If point is within radius, add to results
  if (distance <= radius_km) {
    result_array_add(results, *point, distance);
  }
  
  // Determine which side to search
  dim = node->split_dim;
  
  if (dim == 0) { // latitude dimension
    diff = center_lat - point->lat;
    coord_diff = fabs(diff);
  } else { // longitude dimension
    diff = center_lon - point->lon;
    coord_diff = fabs(diff);
  }
  
  // Search the near side first
  if (diff < 0) {
    search_radius_recursive(tree, node->left_idx, center_lat, cos_center_lat, center_lon, 
			    radius_km, results, max_results);
    // Check if we need to search the far side (rough approximation: 1 degree ≈ 111 km)
    if (coord_diff * 111.0 <= radius_km) {
      search_radius_recursive(tree, node->right_idx, center_lat, cos_center_lat, center_lon, 
			      radius_km, results, max_results);
    }
  } else {
    search_radius_recursive(tree, node->right_idx, center_lat, cos_center_lat, center_lon, 
			    radius_km, results, max_results);
    // Check if we need to search the far side
    if (coord_diff * 111.0 <= radius_km) {
      search_radius_recursive(tree, node->left_idx, center_lat, cos_center_lat, center_lon, 
			      radius_km, results, max_results);
    }
  }
}

result_array_t *geo_tree_query_k_nearest(geo_tree_t *tree, BC_CPREC center_lat, BC_CPREC center_lon, int k) {
  int best_count = 0;
  query_result_t *best_results;
  result_array_t *results;
  
  if (!tree || tree->root_idx == -1 || k <= 0) return NULL;
  
  best_results = (query_result_t *)malloc(sizeof(query_result_t) * k);
  if (!best_results) return NULL;
  
  best_count = 0;
    
  search_knn_recursive(tree, tree->root_idx, center_lat, center_lon, best_results, &best_count, k);
  
  // Sort results by distance
  qsort(best_results, best_count, sizeof(query_result_t), compare_by_distance);
  
  // Create result array
  results = result_array_create(best_count);
  if (results) {
    for (int i = 0; i < best_count; i++) {
      result_array_add(results, best_results[i].point, best_results[i].distance_km);
    }
  }
  
  free(best_results);
  return results;
}

void search_knn_recursive(geo_tree_t *tree, int node_idx, BC_CPREC center_lat, BC_CPREC center_lon,
			  query_result_t *best_results, int *best_count, int k) {
  kd_node_t *node;
  kd_tree_point *point;
  BC_CPREC distance ;
  int dim ;
  BC_CPREC diff;
  BC_CPREC coord_diff;
    
  if (node_idx == -1) return;
    
  node = &tree->nodes[node_idx];
  point = &tree->points[node->point_idx];

  distance = distance_geo(center_lon, center_lat, point->lon, point->lat,cos(center_lat),point->cos_lat);
  
  // Add to best results if we have room or if closer than furthest point
  if (*best_count < k) {
    best_results[*best_count].point = *point;
    best_results[*best_count].distance_km = distance;
    (*best_count)++;
    
    // Sort to keep furthest point at the end
    if (*best_count == k) {
      qsort(best_results, *best_count, sizeof(query_result_t), compare_by_distance);
    }
  } else if (distance < best_results[k-1].distance_km) {
    best_results[k-1].point = *point;
    best_results[k-1].distance_km = distance;
    
    // Re-sort to maintain order
    qsort(best_results, k, sizeof(query_result_t), compare_by_distance);
  }
  
  // Determine search order
  dim = node->split_dim;
  diff = (dim == 0) ? (center_lat - point->lat) : (center_lon - point->lon);
  coord_diff = fabs(diff);
  
  // Search near side first
  if (diff < 0) {
    search_knn_recursive(tree, node->left_idx, center_lat, center_lon, best_results, best_count, k);
    // Check if we need to search far side
    if (*best_count < k || coord_diff * 111.0 < best_results[(*best_count < k) ? *best_count - 1 : k - 1].distance_km) {
      search_knn_recursive(tree, node->right_idx, center_lat, center_lon, best_results, best_count, k);
    }
  } else {
    search_knn_recursive(tree, node->right_idx, center_lat, center_lon, best_results, best_count, k);
    // Check if we need to search far side
    if (*best_count < k || coord_diff * 111.0 < best_results[(*best_count < k) ? *best_count - 1 : k - 1].distance_km) {
      search_knn_recursive(tree, node->left_idx, center_lat, center_lon, best_results, best_count, k);
    }
  }
}



// Result array functions
result_array_t *result_array_create(int initial_capacity) {
  result_array_t *arr = (result_array_t *)malloc(sizeof(result_array_t));
  if (!arr) return NULL;
  
  arr->results = (query_result_t *)malloc(sizeof(query_result_t) * initial_capacity);
  if (!arr->results) {
    free(arr);
    return NULL;
  }
  
  arr->count = 0;
  arr->capacity = initial_capacity;
  return arr;
}

void result_array_add(result_array_t *arr, kd_tree_point point, BC_CPREC distance) {
  if (!arr) return;
  
  if (arr->count >= arr->capacity) {
    int new_capacity = arr->capacity * 2;
    query_result_t* new_results = (query_result_t*)realloc(arr->results, 
							   sizeof(query_result_t) * new_capacity);
    if (!new_results) return;
    
    arr->results = new_results;
    arr->capacity = new_capacity;
  }
  
  arr->results[arr->count].point = point;
  arr->results[arr->count].distance_km = distance;
  arr->count++;
}

void result_array_destroy(result_array_t *results) {
  if (!results) return;
  
  if (results->results) {
        free(results->results);
  }
  free(results);
}

// Comparison functions for qsort
int compare_by_val(const void* a, const void* b) {
  typedef struct { int idx; BC_CPREC value; } sort_item_t;
  const sort_item_t* ia = (const sort_item_t*)a;
  const sort_item_t* ib = (const sort_item_t*)b;
  
  if (ia->value < ib->value) return -1;
  if (ia->value > ib->value) return 1;
  return 0;
}


int compare_by_distance(const void* a, const void* b) {
  const query_result_t* ra = (const query_result_t*)a;
  const query_result_t* rb = (const query_result_t*)b;
  
  if (ra->distance_km < rb->distance_km) return -1;
  if (ra->distance_km > rb->distance_km) return 1;
  return 0;
}

void print_results(result_array_t* results, const char* title) {
  query_result_t *result;
  
  if (!results) {
    printf("No results\n");
    return;
  }
  
  printf("\n%s:\n", title);
  printf("Found %d results:\n", results->count);
  
  for (int i = 0; i < results->count; i++) {
    result = &results->results[i];
    printf("  %d.  %.1f km\n", i+1, result->distance_km);
  }
}
