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
  int code;
} geo_point_t;

// Structure for query results
typedef struct {
    geo_point_t point;
    double distance_km;
} query_result_t;

// Simple geographic search structure (no complex tree for now)
typedef struct {
    geo_point_t* points;        // Array of all points
    int num_points;
    int capacity;
} geo_search_t;

// Dynamic array for storing query results
typedef struct {
    query_result_t* results;
    int count;
    int capacity;
} result_array_t;
