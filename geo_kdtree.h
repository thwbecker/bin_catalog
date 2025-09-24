/* 
   
   headers and defines for geo_kdtree, based on Claude generated code
   

 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Structure to represent a geographic point
typedef struct {
    double lat;
    double lon;
  int code;  // Optional user data pointer
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

// Index/point pair for sorting
typedef struct {
    int idx;
    geo_point_t* point;
} idx_point_t;

