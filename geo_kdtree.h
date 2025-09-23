/* 
   
   headers and defines for geo_kdtree, based on Claude generated code
   

 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


typedef struct {
  BC_CPREC lat,lon,cos_lat;
  int code;
  //void *dpointer;
} kd_tree_point;


//Create a structure for sorting with indices
typedef struct {
  int idx;
  BC_CPREC value;
} sort_item_t;
    
// Structure for query results
typedef struct {
  kd_tree_point point;
  BC_CPREC distance_km;
} query_result_t;

// KD-tree node structure - optimized for cache performance
typedef struct kd_node {
    int point_idx;              // Index into the points array
    int left_idx;               // Index of left child in node array (-1 if none)
    int right_idx;              // Index of right child in node array (-1 if none)
    char split_dim;             // Split dimension (0=lat, 1=lon)
} kd_node_t;

// Main tree structure - optimized for large datasets
typedef struct {
    kd_tree_point* points;        // Array of all points
    int num_points;
    int capacity;
    kd_node_t* nodes;           // Array-based tree nodes for better cache performance
    int num_nodes;
    int node_capacity;
    int root_idx;               // Index of root node (-1 if none)
} geo_tree_t;

// Dynamic array for storing query results
typedef struct {
    query_result_t* results;
    int count;
    int capacity;
} result_array_t;

// Function declarations
geo_tree_t *geo_tree_create(int );
void geo_tree_destroy(geo_tree_t *);
int geo_tree_add_point(geo_tree_t *, BC_CPREC, BC_CPREC, int);
void geo_tree_build(geo_tree_t *);
result_array_t *geo_tree_query_radius(geo_tree_t *, BC_CPREC , BC_CPREC , BC_CPREC,
				      BC_CPREC , int );
result_array_t *geo_tree_query_k_nearest(geo_tree_t *, BC_CPREC , BC_CPREC , int );
void result_array_destroy(result_array_t *);

// Utility functions

void print_results(result_array_t * , const char * );

// Internal helper functions
int build_tree_recursive(geo_tree_t * , int * , int , int );
void search_radius_recursive(geo_tree_t * , int , BC_CPREC , BC_CPREC ,BC_CPREC,
                            BC_CPREC , result_array_t * , int _);
void search_knn_recursive(geo_tree_t * , int , BC_CPREC , BC_CPREC ,
			  query_result_t * , int * , int );
int allocate_node(geo_tree_t * );
int compare_by_val(const void * , const void * );
int compare_by_distance(const void * , const void * );
result_array_t *result_array_create(int );
void result_array_add(result_array_t * , kd_tree_point , BC_CPREC );
