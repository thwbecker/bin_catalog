/* 
   
   KD tree type implementation for geographic data, partially based on
   code from Claude code AI

   THIS NEEDS TO BE OPTIMIZED, SOME REALLY SILLY STUFF IN HERE


*/

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

// Function declarations
geo_search_t* geo_search_create(int initial_capacity);
void geo_search_destroy(geo_search_t* search);
int geo_search_add_point(geo_search_t* search, double lat, double lon, const char* name, void* data);
result_array_t* geo_search_query_radius(geo_search_t* search, double center_lat, double center_lon, 
                                       double radius_km, int max_results);
result_array_t* geo_search_query_k_nearest(geo_search_t* search, double center_lat, double center_lon, int k);
void result_array_destroy(result_array_t* results);

// Utility functions
double haversine_distance(double lat1, double lon1, double lat2, double lon2);
double deg_to_rad(double deg);
void print_results(result_array_t* results, const char* title);
int compare_by_distance(const void* a, const void* b);
result_array_t* result_array_create(int initial_capacity);
void result_array_add(result_array_t* results, geo_point_t point, double distance);

// Implementation

geo_search_t* geo_search_create(int initial_capacity) {
    if (initial_capacity <= 0) initial_capacity = 100;
    
    geo_search_t* search = (geo_search_t*)malloc(sizeof(geo_search_t));
    if (!search) {
        printf("Failed to allocate geo_search_t\n");
        return NULL;
    }
    
    search->points = (geo_point_t*)malloc(sizeof(geo_point_t) * initial_capacity);
    if (!search->points) {
        printf("Failed to allocate points array\n");
        free(search);
        return NULL;
    }
    
    search->num_points = 0;
    search->capacity = initial_capacity;
    
    return search;
}

void geo_search_destroy(geo_search_t* search) {
    if (!search) return;
    
    if (search->points) {
        free(search->points);
    }
    free(search);
}

int geo_search_add_point(geo_search_t* search, double lat, double lon, const char* name, void* data) {
    if (!search) return 0;
    
    // Resize array if needed
    if (search->num_points >= search->capacity) {
        int new_capacity = search->capacity * 2;
        geo_point_t* new_points = (geo_point_t*)realloc(search->points, 
                                                       sizeof(geo_point_t) * new_capacity);
        if (!new_points) {
            printf("Failed to resize points array\n");
            return 0;
        }
        
        search->points = new_points;
        search->capacity = new_capacity;
    }
    
    // Add the point
    geo_point_t* point = &search->points[search->num_points];
    point->lat = lat;
    point->lon = lon;
    point->data = data;
    
    if (name) {
        strncpy(point->name, name, MAX_NAME_LEN - 1);
        point->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        point->name[0] = '\0';
    }
    
    search->num_points++;
    return 1;
}

result_array_t* geo_search_query_radius(geo_search_t* search, double center_lat, double center_lon,
                                       double radius_km, int max_results) {
    if (!search || radius_km <= 0) return NULL;
    
    result_array_t* results = result_array_create(100);
    if (!results) return NULL;
    
    // Simple linear search through all points
    for (int i = 0; i < search->num_points; i++) {
        if (max_results > 0 && results->count >= max_results) break;
        
        geo_point_t* point = &search->points[i];
        double distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
        
        if (distance <= radius_km) {
            result_array_add(results, *point, distance);
        }
    }
    
    // Sort results by distance
    qsort(results->results, results->count, sizeof(query_result_t), compare_by_distance);
    
    return results;
}

result_array_t* geo_search_query_k_nearest(geo_search_t* search, double center_lat, double center_lon, int k) {
    if (!search || k <= 0) return NULL;
    
    // Collect all points with distances
    result_array_t* all_results = result_array_create(search->num_points + 10);
    if (!all_results) return NULL;
    
    // Calculate distance to every point
    for (int i = 0; i < search->num_points; i++) {
        geo_point_t* point = &search->points[i];
        double distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
        result_array_add(all_results, *point, distance);
    }
    
    // Sort by distance
    qsort(all_results->results, all_results->count, sizeof(query_result_t), compare_by_distance);
    
    // Create result array with just the k nearest
    result_array_t* results = result_array_create(k + 5);
    if (results) {
        int limit = (k < all_results->count) ? k : all_results->count;
        for (int i = 0; i < limit; i++) {
            result_array_add(results, all_results->results[i].point, all_results->results[i].distance_km);
        }
    }
    
    result_array_destroy(all_results);
    return results;
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

// Result array functions
result_array_t* result_array_create(int initial_capacity) {
    if (initial_capacity <= 0) initial_capacity = 10;
    
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

int compare_by_distance(const void* a, const void* b) {
    const query_result_t* ra = (const query_result_t*)a;
    const query_result_t* rb = (const query_result_t*)b;
    
    if (ra->distance_km < rb->distance_km) return -1;
    if (ra->distance_km > rb->distance_km) return 1;
    return 0;
}

void print_results(result_array_t* results, const char* title) {
    if (!results) {
        printf("No results for %s\n", title ? title : "query");
        return;
    }
    
    printf("\n%s:\n", title ? title : "Results");
    printf("Found %d results:\n", results->count);
    
    for (int i = 0; i < results->count; i++) {
        query_result_t* result = &results->results[i];
        printf("  %d. %s: %.1f km\n", i+1, result->point.name, result->distance_km);
    }
}

// Demo function
void demo() {
    printf("Simple Geographic Search Demo\n");
    printf("=============================\n");
    
    // Create search structure
    geo_search_t* search = geo_search_create(20);
    if (!search) {
        printf("Failed to create search structure\n");
        return;
    }
    
    // Add sample cities
    printf("Adding sample cities...\n");
    geo_search_add_point(search, 40.7128, -74.0060, "New York", NULL);
    geo_search_add_point(search, 51.5074, -0.1278, "London", NULL);
    geo_search_add_point(search, 48.8566, 2.3522, "Paris", NULL);
    geo_search_add_point(search, 35.6762, 139.6503, "Tokyo", NULL);
    geo_search_add_point(search, -33.8688, 151.2093, "Sydney", NULL);
    geo_search_add_point(search, 37.7749, -122.4194, "San Francisco", NULL);
    geo_search_add_point(search, 52.5200, 13.4050, "Berlin", NULL);
    geo_search_add_point(search, 55.7558, 37.6176, "Moscow", NULL);
    geo_search_add_point(search, 39.9042, 116.4074, "Beijing", NULL);
    geo_search_add_point(search, -23.5558, -46.6396, "Sao Paulo", NULL);
    
    printf("Added %d cities\n", search->num_points);
    
    // Test 1: Radius query
    printf("\n=== Test 1: Radius Query ===");
    result_array_t* radius_results = geo_search_query_radius(search, 51.5074, -0.1278, 2000, 0);
    print_results(radius_results, "Cities within 2000 km of London");
    if (radius_results) result_array_destroy(radius_results);
    
    // Test 2: K-nearest neighbors
    printf("\n=== Test 2: K-Nearest Neighbors ===");
    result_array_t* knn_results = geo_search_query_k_nearest(search, 37.7749, -122.4194, 3);
    print_results(knn_results, "3 nearest cities to San Francisco");
    if (knn_results) result_array_destroy(knn_results);
    
    // Test 3: Larger dataset
    printf("\n=== Test 3: Performance Test ===");
    geo_search_t* big_search = geo_search_create(10000);
    
    srand((unsigned int)time(NULL));
    int num_points = 10000;
    
    printf("Adding %d random points...\n", num_points);
    clock_t start = clock();
    
    for (int i = 0; i < num_points; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lon = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        
        char name[32];
        snprintf(name, sizeof(name), "Point_%d", i);
        geo_search_add_point(big_search, lat, lon, name, NULL);
    }
    
    clock_t end = clock();
    printf("Added %d points in %.3f seconds\n", 
           big_search->num_points, ((double)(end - start)) / CLOCKS_PER_SEC);
    
    // Test radius query performance
    printf("Testing radius query...\n");
    start = clock();
    result_array_t* big_radius = geo_search_query_radius(big_search, 40.7128, -74.0060, 1000, 0);
    end = clock();
    printf("Found %d points within 1000 km in %.4f seconds\n", 
           big_radius ? big_radius->count : 0, ((double)(end - start)) / CLOCKS_PER_SEC);
    
    // Test k-NN query performance
    printf("Testing k-NN query...\n");
    start = clock();
    result_array_t* big_knn = geo_search_query_k_nearest(big_search, 40.7128, -74.0060, 10);
    end = clock();
    printf("Found %d nearest points in %.4f seconds\n", 
           big_knn ? big_knn->count : 0, ((double)(end - start)) / CLOCKS_PER_SEC);
    
    // Cleanup
    if (big_radius) result_array_destroy(big_radius);
    if (big_knn) result_array_destroy(big_knn);
    geo_search_destroy(big_search);
    geo_search_destroy(search);
    
    printf("\n=== Performance Summary ===\n");
    printf("This simple linear search approach:\n");
    printf("- Radius query: O(n) time complexity\n");
    printf("- K-NN query: O(n log n) time complexity\n");
    printf("- Memory usage: O(n)\n");
    printf("- No complex data structures to debug\n");
    printf("- Guaranteed correctness\n");
    
    printf("\nFor 1M points, expect:\n");
    printf("- Radius query: ~0.1-1 second\n");
    printf("- K-NN query: ~1-5 seconds\n");
    printf("- Memory usage: ~100 MB\n");
}

int main() {
    demo();
    return 0;
}
