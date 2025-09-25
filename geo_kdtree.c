#include "catalog.h"
/* 
   
   KD tree type implementation for geographic data, partially based on
   code from Claude code AI

   THIS NEEDS TO BE OPTIMIZED, SOME REALLY SILLY STUFF IN HERE


*/




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

int geo_search_add_point(geo_search_t* search, BC_CPREC lat, BC_CPREC lon, int code) {
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
    point->code = code;
    
    
    search->num_points++;
    return 1;
}

result_array_t* geo_search_query_radius(geo_search_t* search, BC_CPREC center_lat, BC_CPREC center_lon,
                                       BC_CPREC radius_km, int max_results) {
    if (!search || radius_km <= 0) return NULL;
    
    result_array_t* results = result_array_create(100);
    if (!results) return NULL;
    
    // Simple linear search through all points
    for (int i = 0; i < search->num_points; i++) {
        if (max_results > 0 && results->count >= max_results) break;
        
        geo_point_t* point = &search->points[i];
        BC_CPREC distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
        
        if (distance <= radius_km) {
            result_array_add(results, *point, distance);
        }
    }
    
    // Sort results by distance
    qsort(results->results, results->count, sizeof(query_result_t), compare_by_distance);
    
    return results;
}

result_array_t* geo_search_query_k_nearest(geo_search_t* search, BC_CPREC center_lat, BC_CPREC center_lon, int k) {
    if (!search || k <= 0) return NULL;
    
    // Collect all points with distances
    result_array_t* all_results = result_array_create(search->num_points + 10);
    if (!all_results) return NULL;
    
    // Calculate distance to every point
    for (int i = 0; i < search->num_points; i++) {
        geo_point_t* point = &search->points[i];
        BC_CPREC distance = haversine_distance(center_lat, center_lon, point->lat, point->lon);
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
        printf("  %d. %i: %.1f km\n", i+1, result->point.code, result->distance_km);
    }
}
