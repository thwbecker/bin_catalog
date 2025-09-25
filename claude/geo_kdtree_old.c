#include "catalog.h"


/* 
   
   KD tree type implementation for geographic data, partially based on
   code from Claude code AI

   THIS NEEDS TO BE OPTIMIZED, SOME REALLY SILLY STUFF IN HERE


*/

geo_tree_t* geo_tree_create(int initial_capacity) {
    geo_tree_t* tree = (geo_tree_t*)malloc(sizeof(geo_tree_t));
    if (!tree) return NULL;
    
    tree->points = (geo_point_t*)malloc(sizeof(geo_point_t) * initial_capacity);
    if (!tree->points) {
        free(tree);
        return NULL;
    }
    
    tree->num_points = 0;
    tree->capacity = initial_capacity;
    tree->root = NULL;
    
    return tree;
}

void geo_tree_destroy(geo_tree_t* tree) {
    if (!tree) return;
    
    if (tree->points) {
        free(tree->points);
    }
    
    destroy_tree_recursive(tree->root);
    free(tree);
}

void destroy_tree_recursive(kd_node_t* node) {
    if (!node) return;
    
    destroy_tree_recursive(node->left);
    destroy_tree_recursive(node->right);
    free(node);
}

int geo_tree_add_point(geo_tree_t* tree, BC_CPREC lat, BC_CPREC lon, int code) {
    if (!tree) return 0;
    
    // Resize array if needed
    if (tree->num_points >= tree->capacity) {
        int new_capacity = tree->capacity * 2;
        geo_point_t* new_points = (geo_point_t*)realloc(tree->points, 
                                                       sizeof(geo_point_t) * new_capacity);
        if (!new_points) return 0;
        
        tree->points = new_points;
        tree->capacity = new_capacity;
    }
    
    // Add the point
    geo_point_t* point = &tree->points[tree->num_points];
    point->lat = lat;
    point->lon = lon;
    point->code = code;
    
   
    tree->num_points++;
    return 1;
}

void geo_tree_build(geo_tree_t* tree) {
    if (!tree || tree->num_points == 0) return;
    
    // Destroy existing tree
    destroy_tree_recursive(tree->root);
    tree->root = NULL;
    
    // Create array of indices
    int* indices = (int*)malloc(sizeof(int) * tree->num_points);
    if (!indices) return;
    
    for (int i = 0; i < tree->num_points; i++) {
        indices[i] = i;
    }
    
    // Build the tree
    tree->root = build_tree_recursive(tree, indices, tree->num_points, 0);
    
    free(indices);
}

kd_node_t* build_tree_recursive(geo_tree_t* tree, int* indices, int count, int depth) {
    if (count <= 0) return NULL;
    
    // Sort by current dimension (0 = lat, 1 = lon)
    int dim = depth % 2;
    
    // Create array of index/point pairs for sorting
    idx_point_t* items = (idx_point_t*)malloc(sizeof(idx_point_t) * count);
    if (!items) return NULL;
    
    for (int i = 0; i < count; i++) {
        items[i].idx = indices[i];
        items[i].point = &tree->points[indices[i]];
    }
    
    // Sort by the current dimension
    if (dim == 0) {
        qsort(items, count, sizeof(idx_point_t), compare_points_by_lat);
    } else {
        qsort(items, count, sizeof(idx_point_t), compare_points_by_lon);
    }
    
    // Update indices array with sorted order
    for (int i = 0; i < count; i++) {
        indices[i] = items[i].idx;
    }
    free(items);
    
    // Choose median
    int median = count / 2;
    
    // Create node
    kd_node_t* node = (kd_node_t*)malloc(sizeof(kd_node_t));
    if (!node) return NULL;
    
    node->point_idx = indices[median];
    node->left = build_tree_recursive(tree, indices, median, depth + 1);
    node->right = build_tree_recursive(tree, indices + median + 1, count - median - 1, depth + 1);
    
    return node;
}

result_array_t* geo_tree_query_radius(geo_tree_t* tree, BC_CPREC center_lat, BC_CPREC center_lon,
                                     BC_CPREC radius_km, int max_results) {
    if (!tree || !tree->root) return NULL;
    
    result_array_t* results = result_array_create(100);
    if (!results) return NULL;
    
    search_radius_recursive(tree, tree->root, center_lat, center_lon, 
                           radius_km, results, 0, max_results);
    
    // Sort results by distance
    qsort(results->results, results->count, sizeof(query_result_t), compare_by_distance);
    
    // Truncate to max_results if specified
    if (max_results > 0 && results->count > max_results) {
        results->count = max_results;
    }
    
    return results;
}

void search_radius_recursive(geo_tree_t* tree, kd_node_t* node, BC_CPREC center_lat, BC_CPREC center_lon,
                            BC_CPREC radius_km, result_array_t* results, int depth, int max_results) {
    if (!node) return;
    if (max_results > 0 && results->count >= max_results) return;
    
    geo_point_t* point = &tree->points[node->point_idx];
    
    // Fast bounding box pre-check to avoid expensive haversine when obviously outside
    BC_CPREC lat_diff = fabs(center_lat - point->lat);
    BC_CPREC lon_diff = fabs(center_lon - point->lon);
    
    // Handle longitude wraparound
    if (lon_diff > 180.0) {
        lon_diff = 360.0 - lon_diff;
    }
    
    // Quick rejection: if either coordinate difference is much larger than radius could allow
    // Use very conservative estimates: 1 degree ≈ 60 km (underestimate to be safe)
    BC_CPREC max_lat_km = lat_diff * 60.0;  // Conservative underestimate
    BC_CPREC max_lon_km = lon_diff * 60.0;  // Conservative underestimate
    
    int definitely_outside = 0;
    
    // Only skip haversine if we're absolutely certain the point is outside
    // Use very conservative bounds (much smaller than actual) to avoid false rejections
    if (max_lat_km > radius_km * 3.0 && max_lon_km > radius_km * 3.0) {
        // Point is very likely outside radius, but still not 100% certain due to spherical geometry
        // So we'll still compute haversine, just note that it's likely to be rejected
        definitely_outside = 1;
    }
    
    // Always compute exact distance - never skip this
    BC_CPREC distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
    
    // If point is within radius, add to results
    if (distance <= radius_km) {
        result_array_add(results, *point, distance);
    }
    
    // ALWAYS search both subtrees - no pruning at all
    // This ensures we never miss any points
    search_radius_recursive(tree, node->left, center_lat, center_lon, 
                           radius_km, results, depth + 1, max_results);
    search_radius_recursive(tree, node->right, center_lat, center_lon, 
                           radius_km, results, depth + 1, max_results);
}

result_array_t* geo_tree_query_k_nearest(geo_tree_t* tree, BC_CPREC center_lat, BC_CPREC center_lon, int k) {
    if (!tree || !tree->root || k <= 0) return NULL;
    
    query_result_t* best_results = (query_result_t*)calloc(k, sizeof(query_result_t));
    if (!best_results) return NULL;
    
    // Initialize with very large distances
    for (int i = 0; i < k; i++) {
        best_results[i].distance_km = INFINITY;
    }
    
    int best_count = 0;
    
    search_knn_recursive(tree, tree->root, center_lat, center_lon, best_results, &best_count, k, 0);
    
    // Sort results by distance and remove any unfilled slots
    qsort(best_results, best_count, sizeof(query_result_t), compare_by_distance);
    
    // Create result array
    result_array_t* results = result_array_create(best_count);
    if (results) {
        for (int i = 0; i < best_count; i++) {
            result_array_add(results, best_results[i].point, best_results[i].distance_km);
        }
    }
    
    free(best_results);
    return results;
}

void search_knn_recursive(geo_tree_t* tree, kd_node_t* node, BC_CPREC center_lat, BC_CPREC center_lon,
                         query_result_t* best_results, int* best_count, int k, int depth) {
    if (!node) return;
    
    geo_point_t* point = &tree->points[node->point_idx];
    BC_CPREC distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
    
    // Add to best results if we have room or if closer than furthest point
    if (*best_count < k) {
        best_results[*best_count].point = *point;
        best_results[*best_count].distance_km = distance;
        (*best_count)++;
        
        // Keep sorted
        if (*best_count > 1) {
            qsort(best_results, *best_count, sizeof(query_result_t), compare_by_distance);
        }
    } else if (distance < best_results[k-1].distance_km) {
        // Replace the furthest point
        best_results[k-1].point = *point;
        best_results[k-1].distance_km = distance;
        
        // Re-sort to maintain order
        qsort(best_results, k, sizeof(query_result_t), compare_by_distance);
    }
    
    // Conservative pruning for k-NN
    int dim = depth % 2;
    BC_CPREC split_value = (dim == 0) ? point->lat : point->lon;
    BC_CPREC query_coord = (dim == 0) ? center_lat : center_lon;
    BC_CPREC coord_diff = fabs(query_coord - split_value);
    
    // Get current k-th best distance
    BC_CPREC kth_best_distance = (*best_count >= k) ? best_results[k-1].distance_km : INFINITY;
    
    // Very conservative distance estimation with large safety margin
    BC_CPREC max_coord_distance = coord_diff * BC_DEG_SCALE; // Worst case for both lat and lon
    BC_CPREC safety_margin = 1.5; // 50% safety margin
    
    // Search both sides unless coordinate distance is much larger than current k-th best
    int search_left = 1;
    int search_right = 1;
    
    if (query_coord <= split_value) {
        // Search left (near) side first
        search_knn_recursive(tree, node->left, center_lat, center_lon, best_results, best_count, k, depth + 1);
        
        // Update k-th best distance after near side search
        kth_best_distance = (*best_count >= k) ? best_results[k-1].distance_km : INFINITY;
        
        // Only skip right side if coordinate distance is much larger than k-th best
        if (max_coord_distance <= kth_best_distance * safety_margin) {
            search_knn_recursive(tree, node->right, center_lat, center_lon, best_results, best_count, k, depth + 1);
        }
    } else {
        // Search right (near) side first
        search_knn_recursive(tree, node->right, center_lat, center_lon, best_results, best_count, k, depth + 1);
        
        // Update k-th best distance after near side search
        kth_best_distance = (*best_count >= k) ? best_results[k-1].distance_km : INFINITY;
        
        // Only skip left side if coordinate distance is much larger than k-th best
        if (max_coord_distance <= kth_best_distance * safety_margin) {
            search_knn_recursive(tree, node->left, center_lat, center_lon, best_results, best_count, k, depth + 1);
        }
    }
}

BC_CPREC haversine_distance(BC_CPREC lat1, BC_CPREC lon1, BC_CPREC lat2, BC_CPREC lon2) {
    BC_CPREC rlat1 = deg_to_rad(lat1);
    BC_CPREC rlon1 = deg_to_rad(lon1);
    BC_CPREC rlat2 = deg_to_rad(lat2);
    BC_CPREC rlon2 = deg_to_rad(lon2);
    
    BC_CPREC dlat = rlat2 - rlat1;
    BC_CPREC dlon = rlon2 - rlon1;
    
    BC_CPREC a = sin(dlat/2) * sin(dlat/2) + 
               cos(rlat1) * cos(rlat2) * sin(dlon/2) * sin(dlon/2);
    BC_CPREC c = 2 * asin(sqrt(a));
    
    return BC_RADIUS * c;
}

BC_CPREC deg_to_rad(BC_CPREC deg) {
    return deg * M_PI / 180.0;
}

// Fast bounding box check for spherical coordinates
int point_in_bounding_box(BC_CPREC center_lat, BC_CPREC center_lon, BC_CPREC radius_km,
                         BC_CPREC point_lat, BC_CPREC point_lon) {
    // Convert radius to approximate degrees
    BC_CPREC lat_delta = radius_km / BC_DEG_SCALE;  // 1 degree lat ≈ 111 km
    
    // Check latitude bounds (simple)
    if (point_lat < center_lat - lat_delta || point_lat > center_lat + lat_delta) {
        return 0;
    }
    
    // For longitude, account for latitude compression
    BC_CPREC lat_rad = deg_to_rad(center_lat);
    BC_CPREC lon_delta = radius_km / (BC_DEG_SCALE * cos(lat_rad));
    
    // Handle longitude wraparound
    BC_CPREC lon_diff = fabs(point_lon - center_lon);
    if (lon_diff > 180.0) {
        lon_diff = 360.0 - lon_diff;
    }
    
    return lon_diff <= lon_delta;
}

// Calculate minimum possible distance from a point to a splitting plane
BC_CPREC min_distance_to_plane(BC_CPREC center_lat, BC_CPREC center_lon, 
                            BC_CPREC split_value, int split_dim) {
    if (split_dim == 0) { // latitude split
        // Minimum distance is along a line of constant longitude
        BC_CPREC lat_diff = fabs(center_lat - split_value);
        return lat_diff * BC_DEG_SCALE; // 1 degree lat ≈ 111 km
    } else { // longitude split
        // Minimum distance depends on latitude
        BC_CPREC lat_rad = deg_to_rad(center_lat);
        BC_CPREC lon_diff = fabs(center_lon - split_value);
        
        // Handle wraparound
        if (lon_diff > 180.0) {
            lon_diff = 360.0 - lon_diff;
        }
        
        return lon_diff * BC_DEG_SCALE * cos(lat_rad);
    }
}

// Result array functions
result_array_t* result_array_create(int initial_capacity) {
    result_array_t* arr = (result_array_t*)malloc(sizeof(result_array_t));
    if (!arr) return NULL;
    
    arr->results = (query_result_t*)malloc(sizeof(query_result_t) * initial_capacity);
    if (!arr->results) {
        free(arr);
        return NULL;
    }
    
    arr->count = 0;
    arr->capacity = initial_capacity;
    return arr;
}

void result_array_add(result_array_t* arr, geo_point_t point, BC_CPREC distance) {
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

void result_array_destroy(result_array_t* results) {
    if (!results) return;
    
    if (results->results) {
        free(results->results);
    }
    free(results);
}

// Comparison functions for qsort
int compare_points_by_lat(const void* a, const void* b) {
    const idx_point_t* ia = (const idx_point_t*)a;
    const idx_point_t* ib = (const idx_point_t*)b;
    
    if (ia->point->lat < ib->point->lat) return -1;
    if (ia->point->lat > ib->point->lat) return 1;
    return 0;
}

int compare_points_by_lon(const void* a, const void* b) {
    const idx_point_t* ia = (const idx_point_t*)a;
    const idx_point_t* ib = (const idx_point_t*)b;
    
    if (ia->point->lon < ib->point->lon) return -1;
    if (ia->point->lon > ib->point->lon) return 1;
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
    if (!results) {
        printf("No results\n");
        return;
    }
    
    printf("\n%s:\n", title);
    printf("Found %d results:\n", results->count);
    
    for (int i = 0; i < results->count; i++) {
        query_result_t* result = &results->results[i];
        printf("  %d. %i: %.1f km\n", i+1, result->point.code, result->distance_km);
    }
}




void print_results_file(result_array_t* results, FILE *out) {
  query_result_t *result;
  
  if (!results) {
    printf("No results\n");
    return;
  }
  
  
  for (int i = 0; i < results->count; i++) {
    result = &results->results[i];
    fprintf(out,"%g %g %d %g\n", BC_R2D(result->point.lon),BC_R2D(result->point.lat),i+1, result->distance_km);
  }
}
