#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define EARTH_RADIUS_KM 6371.0
#define MAX_NAME_LEN 64
#define PI 3.14159265358979323846

// Structure to represent a geographic point
typedef struct {
    double lat;
    double lon;
    char name[MAX_NAME_LEN];
    void* data;  // Optional user data pointer
} geo_point_t;

// Structure for query results
typedef struct {
    geo_point_t point;
    double distance_km;
} query_result_t;

// KD-tree node structure
typedef struct kd_node {
    int point_idx;              // Index into the points array
    struct kd_node* left;
    struct kd_node* right;
} kd_node_t;

// Main tree structure
typedef struct {
    geo_point_t* points;        // Array of all points
    int num_points;
    int capacity;
    kd_node_t* root;
} geo_tree_t;

// Dynamic array for storing query results
typedef struct {
    query_result_t* results;
    int count;
    int capacity;
} result_array_t;

// Function declarations
geo_tree_t* geo_tree_create(int initial_capacity);
void geo_tree_destroy(geo_tree_t* tree);
int geo_tree_add_point(geo_tree_t* tree, double lat, double lon, const char* name, void* data);
void geo_tree_build(geo_tree_t* tree);
result_array_t* geo_tree_query_radius(geo_tree_t* tree, double center_lat, double center_lon, 
                                     double radius_km, int max_results);
result_array_t* geo_tree_query_k_nearest(geo_tree_t* tree, double center_lat, double center_lon, int k);
void result_array_destroy(result_array_t* results);

// Utility functions
double haversine_distance(double lat1, double lon1, double lat2, double lon2);
double deg_to_rad(double deg);
void print_results(result_array_t* results, const char* title);

// Internal helper functions
kd_node_t* build_tree_recursive(geo_tree_t* tree, int* indices, int count, int depth);
void search_radius_recursive(geo_tree_t* tree, kd_node_t* node, double center_lat, double center_lon,
                            double radius_km, result_array_t* results, int depth, int max_results);
void search_knn_recursive(geo_tree_t* tree, kd_node_t* node, double center_lat, double center_lon,
                         query_result_t* best_results, int* best_count, int k, int depth);
void destroy_tree_recursive(kd_node_t* node);
int compare_points_by_lat(const void* a, const void* b);
int compare_points_by_lon(const void* a, const void* b);
int compare_by_distance(const void* a, const void* b);
result_array_t* result_array_create(int initial_capacity);
void result_array_add(result_array_t* results, geo_point_t point, double distance);

// Index/point pair for sorting
typedef struct {
    int idx;
    geo_point_t* point;
} idx_point_t;

// Main implementation

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

int geo_tree_add_point(geo_tree_t* tree, double lat, double lon, const char* name, void* data) {
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
    point->data = data;
    
    if (name) {
        strncpy(point->name, name, MAX_NAME_LEN - 1);
        point->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        point->name[0] = '\0';
    }
    
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

result_array_t* geo_tree_query_radius(geo_tree_t* tree, double center_lat, double center_lon,
                                     double radius_km, int max_results) {
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

void search_radius_recursive(geo_tree_t* tree, kd_node_t* node, double center_lat, double center_lon,
                            double radius_km, result_array_t* results, int depth, int max_results) {
    if (!node) return;
    if (max_results > 0 && results->count >= max_results) return;
    
    geo_point_t* point = &tree->points[node->point_idx];
    
    // Fast bounding box pre-check to avoid expensive haversine when obviously outside
    double lat_diff = fabs(center_lat - point->lat);
    double lon_diff = fabs(center_lon - point->lon);
    
    // Handle longitude wraparound
    if (lon_diff > 180.0) {
        lon_diff = 360.0 - lon_diff;
    }
    
    // Quick rejection: if either coordinate difference is much larger than radius could allow
    // Use very conservative estimates: 1 degree ≈ 60 km (underestimate to be safe)
    double max_lat_km = lat_diff * 60.0;  // Conservative underestimate
    double max_lon_km = lon_diff * 60.0;  // Conservative underestimate
    
    int definitely_outside = 0;
    
    // Only skip haversine if we're absolutely certain the point is outside
    // Use very conservative bounds (much smaller than actual) to avoid false rejections
    if (max_lat_km > radius_km * 3.0 && max_lon_km > radius_km * 3.0) {
        // Point is very likely outside radius, but still not 100% certain due to spherical geometry
        // So we'll still compute haversine, just note that it's likely to be rejected
        definitely_outside = 1;
    }
    
    // Always compute exact distance - never skip this
    double distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
    
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

result_array_t* geo_tree_query_k_nearest(geo_tree_t* tree, double center_lat, double center_lon, int k) {
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

void search_knn_recursive(geo_tree_t* tree, kd_node_t* node, double center_lat, double center_lon,
                         query_result_t* best_results, int* best_count, int k, int depth) {
    if (!node) return;
    
    geo_point_t* point = &tree->points[node->point_idx];
    double distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
    
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
    double split_value = (dim == 0) ? point->lat : point->lon;
    double query_coord = (dim == 0) ? center_lat : center_lon;
    double coord_diff = fabs(query_coord - split_value);
    
    // Get current k-th best distance
    double kth_best_distance = (*best_count >= k) ? best_results[k-1].distance_km : INFINITY;
    
    // Very conservative distance estimation with large safety margin
    double max_coord_distance = coord_diff * 111.0; // Worst case for both lat and lon
    double safety_margin = 1.5; // 50% safety margin
    
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

double haversine_distance(double lat1, double lon1, double lat2, double lon2) {
    double rlat1 = deg_to_rad(lat1);
    double rlon1 = deg_to_rad(lon1);
    double rlat2 = deg_to_rad(lat2);
    double rlon2 = deg_to_rad(lon2);
    
    double dlat = rlat2 - rlat1;
    double dlon = rlon2 - rlon1;
    
    double a = sin(dlat/2) * sin(dlat/2) + 
               cos(rlat1) * cos(rlat2) * sin(dlon/2) * sin(dlon/2);
    double c = 2 * asin(sqrt(a));
    
    return EARTH_RADIUS_KM * c;
}

double deg_to_rad(double deg) {
    return deg * PI / 180.0;
}

// Fast bounding box check for spherical coordinates
int point_in_bounding_box(double center_lat, double center_lon, double radius_km,
                         double point_lat, double point_lon) {
    // Convert radius to approximate degrees
    double lat_delta = radius_km / 111.0;  // 1 degree lat ≈ 111 km
    
    // Check latitude bounds (simple)
    if (point_lat < center_lat - lat_delta || point_lat > center_lat + lat_delta) {
        return 0;
    }
    
    // For longitude, account for latitude compression
    double lat_rad = deg_to_rad(center_lat);
    double lon_delta = radius_km / (111.0 * cos(lat_rad));
    
    // Handle longitude wraparound
    double lon_diff = fabs(point_lon - center_lon);
    if (lon_diff > 180.0) {
        lon_diff = 360.0 - lon_diff;
    }
    
    return lon_diff <= lon_delta;
}

// Calculate minimum possible distance from a point to a splitting plane
double min_distance_to_plane(double center_lat, double center_lon, 
                            double split_value, int split_dim) {
    if (split_dim == 0) { // latitude split
        // Minimum distance is along a line of constant longitude
        double lat_diff = fabs(center_lat - split_value);
        return lat_diff * 111.0; // 1 degree lat ≈ 111 km
    } else { // longitude split
        // Minimum distance depends on latitude
        double lat_rad = deg_to_rad(center_lat);
        double lon_diff = fabs(center_lon - split_value);
        
        // Handle wraparound
        if (lon_diff > 180.0) {
            lon_diff = 360.0 - lon_diff;
        }
        
        return lon_diff * 111.0 * cos(lat_rad);
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

void result_array_add(result_array_t* arr, geo_point_t point, double distance) {
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
        printf("  %d. %s: %.1f km\n", i+1, result->point.name, result->distance_km);
    }
}

// Demo function
void demo() {
    printf("Geographic Tree Search Demo\n");
    printf("===========================\n");
    
    // Create tree
    geo_tree_t* tree = geo_tree_create(10);
    if (!tree) {
        printf("Failed to create tree\n");
        return;
    }
    
    // Add sample cities
    geo_tree_add_point(tree, 40.7128, -74.0060, "New York", NULL);
    geo_tree_add_point(tree, 51.5074, -0.1278, "London", NULL);
    geo_tree_add_point(tree, 48.8566, 2.3522, "Paris", NULL);
    geo_tree_add_point(tree, 35.6762, 139.6503, "Tokyo", NULL);
    geo_tree_add_point(tree, -33.8688, 151.2093, "Sydney", NULL);
    geo_tree_add_point(tree, 37.7749, -122.4194, "San Francisco", NULL);
    geo_tree_add_point(tree, 52.5200, 13.4050, "Berlin", NULL);
    geo_tree_add_point(tree, 55.7558, 37.6176, "Moscow", NULL);
    geo_tree_add_point(tree, 39.9042, 116.4074, "Beijing", NULL);
    geo_tree_add_point(tree, -23.5558, -46.6396, "Sao Paulo", NULL);
    
    printf("Building tree with %d cities...\n", tree->num_points);
    geo_tree_build(tree);
    
    // Test 1: Radius query
    printf("\n1. Cities within 2000 km of London (51.5074, -0.1278):");
    result_array_t* radius_results = geo_tree_query_radius(tree, 51.5074, -0.1278, 2000, 0);
    print_results(radius_results, "Cities within 2000 km of London");
    result_array_destroy(radius_results);
    
    // Test 2: K-nearest neighbors
    printf("\n2. 3 nearest cities to San Francisco (37.7749, -122.4194):");
    result_array_t* knn_results = geo_tree_query_k_nearest(tree, 37.7749, -122.4194, 3);
    print_results(knn_results, "3 nearest cities to San Francisco");
    result_array_destroy(knn_results);
    
    // Test 3: Verify correctness with brute force comparison
    printf("\n3. Verification test - comparing tree search vs brute force:");
    
    // Brute force search within 5000 km of New York
    printf("\nBrute force search within 5000 km of New York:\n");
    for (int i = 0; i < tree->num_points; i++) {
        geo_point_t* p = &tree->points[i];
        double dist = haversine_distance(40.7128, -74.0060, p->lat, p->lon);
        if (dist <= 5000.0) {
            printf("  %s: %.1f km\n", p->name, dist);
        }
    }
    
    // Tree search within 5000 km of New York
    result_array_t* tree_results = geo_tree_query_radius(tree, 40.7128, -74.0060, 5000, 0);
    print_results(tree_results, "Tree search within 5000 km of New York");
    result_array_destroy(tree_results);
    
    // Test 4: Performance test with random points
    printf("\n4. Performance test with random data...\n");
    
    geo_tree_t* big_tree = geo_tree_create(1000);
    
    // Add random points
    srand((unsigned int)time(NULL));
    int num_random = 10000;
    
    printf("Adding %d random points...\n", num_random);
    for (int i = 0; i < num_random; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;  // -90 to 90
        double lon = ((double)rand() / RAND_MAX) * 360.0 - 180.0; // -180 to 180
        
        char name[32];
        snprintf(name, sizeof(name), "Point_%d", i);
        geo_tree_add_point(big_tree, lat, lon, name, NULL);
    }
    
    printf("Building tree...\n");
    clock_t start = clock();
    geo_tree_build(big_tree);
    clock_t end = clock();
    double build_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tree built in %.3f seconds\n", build_time);
    
    // Test query performance
    printf("Querying points within 1000 km of New York...\n");
    start = clock();
    result_array_t* big_results = geo_tree_query_radius(big_tree, 40.7128, -74.0060, 1000, 0);
    end = clock();
    double query_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Found %d points within 1000 km in %.4f seconds\n", 
           big_results ? big_results->count : 0, query_time);
    
    result_array_destroy(big_results);
    geo_tree_destroy(big_tree);
    geo_tree_destroy(tree);
}


// Demo function
void demo2() {
  geo_tree_t *tree;
  FILE *in,*out;
  double dlon,dlat,dist;
  unsigned int code,mode;
  int ndist;
  char fname[300];
  result_array_t *results;
  // Create tree
  tree = geo_tree_create(2000000);
  if (!tree) {
    fprintf(stderr,"Failed to create tree\n");
    return;
  }
  in = fopen("cheng25.aki","r");
  code=0;
  while(fscanf(in,"%lf %lf %*f %*f %*f %*f %*f %*f %*f %*s",&dlon,&dlat)==2 ){
    geo_tree_add_point(tree, dlat,dlon,"test",NULL);
    code++;
  }
  fclose(in);
  fprintf(stderr,"Building tree with %d events...\n", tree->num_points);
  clock_t start = clock();
  geo_tree_build(tree);
  clock_t end = clock();
  double build_time = ((double)(end - start)) / CLOCKS_PER_SEC;
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
	results = geo_tree_query_k_nearest(tree, dlat,  dlon, ndist);
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
    demo2();
    return 0;
}
