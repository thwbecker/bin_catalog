#include "catalog.h"
/* 
   
   some code from Claude code AI

   this is a very basic, slow search implementation - the R tree and
   KDtree versions were flaky - so this is NOT a KDtree 

   we leave as is for backward compatibility with an attempted KDtree
   search which is in claude/*.c

*/



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
  search->sorted = 0;
  
  return search;
}

void geo_search_destroy(geo_search_t* search) {
  if (!search) return;
  
  if (search->points) {
    free(search->points);
  }
  free(search);
}
/* add as degree */
int geo_search_add_point(geo_search_t* search, BC_CPREC dlat, BC_CPREC dlon, int code) {
  geo_point_t *point,*new_points;
  int new_capacity;
  if (!search) return 0;
  
  // Resize array if needed
  if (search->num_points >= search->capacity) {
    new_capacity = (int)((float)search->capacity * 1.5);
    new_points = (geo_point_t*)realloc(search->points, sizeof(geo_point_t) * new_capacity);
    if (!new_points) {
      printf("Failed to resize points array\n");
      return 0;
    }
    search->points = new_points;
    search->capacity = new_capacity;
  }
  
  // Add the point
  point = &search->points[search->num_points];
  point->lat = BC_D2R(dlat);
  point->cos_lat = cos(point->lat);
  point->lon = BC_D2R(dlon);
  point->code = code;
  
  search->sorted = 0;
  search->num_points++;
  return 1;
}

/* order by latitude, then by code so that the order is reproducible */
static int compare_by_lat(const void *a, const void *b)
{
  const geo_point_t *pa = (const geo_point_t *)a, *pb = (const geo_point_t *)b;
  if (pa->lat < pb->lat) return -1;
  if (pa->lat > pb->lat) return  1;
  return (pa->code > pb->code) - (pa->code < pb->code);
}
/* sort points by latitude once, so that queries only visit a latitude band */
static void geo_search_sort(geo_search_t *search)
{
  if (!search->sorted) {
    qsort(search->points, search->num_points, sizeof(geo_point_t), compare_by_lat);
    search->sorted = 1;
  }
}
/* first index with lat >= lat0 in the sorted points */
static int geo_search_lower_bound(const geo_search_t *search, BC_CPREC lat0)
{
  int lo = 0, hi = search->num_points, mid;
  while (lo < hi) {
    mid = (lo + hi) / 2;
    if (search->points[mid].lat < lat0) lo = mid + 1; else hi = mid;
  }
  return lo;
}

result_array_t* geo_search_query_radius(geo_search_t* search, BC_CPREC center_dlat,
					BC_CPREC center_dlon,BC_CPREC radius_km) {
  BC_CPREC center_lat, center_lon, center_cos_lat;
  result_array_t *results;
  geo_point_t *point;
  BC_CPREC distance,dlat_max;
  
  center_lat = BC_D2R(center_dlat);
  center_lon = BC_D2R(center_dlon);
  center_cos_lat = cos(center_lat);
  
  if (!search || radius_km <= 0)
    return NULL;
  
  results = result_array_create(100);
  if (!results) return NULL;
  
  /* the great circle distance is at least BC_RADIUS * |dlat|, so
     only points within this latitude band need to be checked. the
     points are sorted by latitude, the band is found by bisection */
  dlat_max = radius_km / BC_RADIUS;
  geo_search_sort(search);
  for (int i = geo_search_lower_bound(search, center_lat - dlat_max);
       i < search->num_points; i++) {
    point = &search->points[i];
    if(point->lat > center_lat + dlat_max)
      break;
    distance =  distance_geo(center_lon, center_lat,point->lon, point->lat,center_cos_lat,point->cos_lat);
    //fprintf(stderr,"%g %g %g %g - %g %g\n",center_lon, center_lat,point->lon, point->lat,distance,radius_km);
    if (distance <= radius_km) {
      result_array_add(results, *point, distance);
    }
  }
  
  // Sort results by distance
  qsort(results->results, results->count, sizeof(query_result_t), compare_by_distance);
 
  return results;
}

/* 
   partial selection (Hoare quickselect): on return, r[0..k-1] hold the
   k smallest distances in arbitrary order
*/
static void select_k_smallest(query_result_t *r, int n, int k)
{
  int left = 0, right = n - 1, i, j;
  BC_CPREC pivot;
  query_result_t tmp;
  if((k <= 0)||(k >= n))
    return;
  while(left < right){
    pivot = r[(left + right)/2].distance_km;
    i = left; j = right;
    while(i <= j){
      while(r[i].distance_km < pivot) i++;
      while(r[j].distance_km > pivot) j--;
      if(i <= j){
	tmp = r[i]; r[i] = r[j]; r[j] = tmp;
	i++; j--;
      }
    }
    /* now r[left..j] <= pivot <= r[i..right] */
    if(k - 1 <= j)
      right = j;
    else if(k - 1 >= i)
      left = i;
    else
      break;
  }
}

/*
   the k nearest points within radius_km, sorted by distance. fewer
   than k are returned if there are fewer within the radius. same set as
   geo_search_query_k_nearest when its k-th distance is <= radius_km,
   but only the latitude band of the radius is searched
*/
result_array_t* geo_search_query_k_nearest_within(geo_search_t* search, BC_CPREC center_dlat,
						   BC_CPREC center_dlon, int k, BC_CPREC radius_km)
{
  result_array_t *results;
  if (!search || k <= 0) return NULL;
  results = geo_search_query_radius(search, center_dlat, center_dlon, radius_km);
  if (results && results->count > k)
    results->count = k;	/* sorted by distance */
  return results;
}

result_array_t* geo_search_query_k_nearest(geo_search_t* search,BC_CPREC center_dlat, BC_CPREC center_dlon, int k)
{
  BC_CPREC center_lat, center_lon, center_cos_lat,distance;
  result_array_t *all_results, *results;
  geo_point_t *point;
  
  center_lat = BC_D2R(center_dlat);
  center_lon = BC_D2R(center_dlon);
  center_cos_lat = cos(center_lat);
  
  if (!search || k <= 0) return NULL;
  
  // Collect all points with distances
  all_results = result_array_create(search->num_points + 10);
  if (!all_results) return NULL;
  
  // Calculate distance to every point
  for (int i = 0; i < search->num_points; i++) {
    point = &search->points[i];
    distance =  distance_geo( center_lon,center_lat,point->lon, point->lat,center_cos_lat,  point->cos_lat);
    
    result_array_add(all_results, *point, distance);
  }
  
  /* move the k smallest distances to the front (average linear
     time), then sort only those */
  if(k < all_results->count){
    select_k_smallest(all_results->results, all_results->count, k);
    qsort(all_results->results, k, sizeof(query_result_t), compare_by_distance);
  }else{
    qsort(all_results->results, all_results->count, sizeof(query_result_t), compare_by_distance);
  }
  
  // Create result array with just the k nearest
  results = result_array_create(k + 5);
  if (results) {
    int limit = (k < all_results->count) ? k : all_results->count;
    for (int i = 0; i < limit; i++) {
      result_array_add(results, all_results->results[i].point, all_results->results[i].distance_km);
    }
  }
  
  result_array_destroy(all_results);
  return results;
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
    int new_capacity = (int)((float)arr->capacity * 1.5);
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

void print_results(result_array_t* results) {
  
  for (int i = 0; i < results->count; i++) {
    query_result_t* result = &results->results[i];
    printf("%11g %11g %11d %05i %8.1f km\n", BC_R2D(result->point.lon),BC_R2D(result->point.lat),i+1, result->point.code, result->distance_km);
  }
}
