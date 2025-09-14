/* 

   from Claude
   
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {
    double x, y;
    int id;
} Point;

typedef struct {
    Point **points;
    int count;
    int capacity;
} GridCell;

typedef struct {
    GridCell *cells;
    int grid_width, grid_height;
    double min_x, min_y, max_x, max_y;
    double cell_size;
    int total_points;
} SpatialGrid;

typedef struct {
    Point *points;
    int count;
    int capacity;
} PointList;

// Initialize a point list for query results
PointList* point_list_create(int initial_capacity) {
    PointList *list = malloc(sizeof(PointList));
    list->points = malloc(sizeof(Point) * initial_capacity);
    list->count = 0;
    list->capacity = initial_capacity;
    return list;
}

void point_list_add(PointList *list, Point *point) {
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        list->points = realloc(list->points, sizeof(Point) * list->capacity);
    }
    list->points[list->count++] = *point;
}

void point_list_free(PointList *list) {
    free(list->points);
    free(list);
}

// Create spatial grid
SpatialGrid* spatial_grid_create(double min_x, double min_y, double max_x, double max_y, double cell_size) {
    SpatialGrid *grid = malloc(sizeof(SpatialGrid));
    
    grid->min_x = min_x;
    grid->min_y = min_y;
    grid->max_x = max_x;
    grid->max_y = max_y;
    grid->cell_size = cell_size;
    grid->total_points = 0;
    
    // Calculate grid dimensions
    grid->grid_width = (int)ceil((max_x - min_x) / cell_size);
    grid->grid_height = (int)ceil((max_y - min_y) / cell_size);
    
    // Allocate grid cells
    int total_cells = grid->grid_width * grid->grid_height;
    grid->cells = calloc(total_cells, sizeof(GridCell));
    
    // Initialize each cell
    for (int i = 0; i < total_cells; i++) {
        grid->cells[i].points = NULL;
        grid->cells[i].count = 0;
        grid->cells[i].capacity = 0;
    }
    
    return grid;
}

// Get grid cell index from coordinates
int get_cell_index(SpatialGrid *grid, double x, double y) {
    int grid_x = (int)floor((x - grid->min_x) / grid->cell_size);
    int grid_y = (int)floor((y - grid->min_y) / grid->cell_size);
    
    // Clamp to grid bounds
    if (grid_x < 0) grid_x = 0;
    if (grid_x >= grid->grid_width) grid_x = grid->grid_width - 1;
    if (grid_y < 0) grid_y = 0;
    if (grid_y >= grid->grid_height) grid_y = grid->grid_height - 1;
    
    return grid_y * grid->grid_width + grid_x;
}

// Add point to grid cell
void grid_cell_add_point(GridCell *cell, Point *point) {
    if (cell->count >= cell->capacity) {
        cell->capacity = cell->capacity == 0 ? 4 : cell->capacity * 2;
        cell->points = realloc(cell->points, sizeof(Point*) * cell->capacity);
    }
    cell->points[cell->count++] = point;
}

// Insert point into spatial grid
void spatial_grid_insert(SpatialGrid *grid, Point *point) {
    int cell_index = get_cell_index(grid, point->x, point->y);
    grid_cell_add_point(&grid->cells[cell_index], point);
    grid->total_points++;
}

// Calculate distance between two points
double point_distance(Point *a, Point *b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    return sqrt(dx * dx + dy * dy);
}

// Query points within distance
PointList* spatial_grid_query_range(SpatialGrid *grid, double x, double y, double max_distance) {
    PointList *results = point_list_create(32);
    
    // Calculate grid cell range to check
    int cells_to_check = (int)ceil(max_distance / grid->cell_size) + 1;
    int center_x = (int)floor((x - grid->min_x) / grid->cell_size);
    int center_y = (int)floor((y - grid->min_y) / grid->cell_size);
    
    Point query_point = {x, y, -1};
    
    // Check all cells within the range
    for (int dy = -cells_to_check; dy <= cells_to_check; dy++) {
        for (int dx = -cells_to_check; dx <= cells_to_check; dx++) {
            int check_x = center_x + dx;
            int check_y = center_y + dy;
            
            // Skip cells outside grid bounds
            if (check_x < 0 || check_x >= grid->grid_width ||
                check_y < 0 || check_y >= grid->grid_height) {
                continue;
            }
            
            int cell_index = check_y * grid->grid_width + check_x;
            GridCell *cell = &grid->cells[cell_index];
            
            // Check all points in this cell
            for (int i = 0; i < cell->count; i++) {
                Point *point = cell->points[i];
                double distance = point_distance(&query_point, point);
                
                if (distance <= max_distance) {
                    point_list_add(results, point);
                }
            }
        }
    }
    
    return results;
}

// Free spatial grid memory
void spatial_grid_free(SpatialGrid *grid) {
    int total_cells = grid->grid_width * grid->grid_height;
    for (int i = 0; i < total_cells; i++) {
        free(grid->cells[i].points);
    }
    free(grid->cells);
    free(grid);
}

// Example usage
int main() {
    // Create a grid covering area from (0,0) to (1000,1000) with 10-unit cells
    SpatialGrid *grid = spatial_grid_create(0.0, 0.0, 1000.0, 1000.0, 10.0);
    
    // Create some sample points
    Point points[10];
    for (int i = 0; i < 10; i++) {
        points[i].x = (double)(rand() % 1000);
        points[i].y = (double)(rand() % 1000);
        points[i].id = i;
        spatial_grid_insert(grid, &points[i]);
    }
    
    printf("Inserted %d points into grid\n", grid->total_points);
    printf("Grid dimensions: %dx%d cells\n", grid->grid_width, grid->grid_height);
    
    // Query for points near (500, 500) within distance 50
    PointList *nearby = spatial_grid_query_range(grid, 500.0, 500.0, 50.0);
    
    printf("Found %d points within distance 50 of (500, 500):\n", nearby->count);
    for (int i = 0; i < nearby->count; i++) {
        Point *p = &nearby->points[i];
        double dist = point_distance(&(Point){500.0, 500.0, -1}, p);
        printf("  Point %d: (%.1f, %.1f) distance: %.2f\n", 
               p->id, p->x, p->y, dist);
    }
    
    // Cleanup
    point_list_free(nearby);
    spatial_grid_free(grid);
    
    return 0;
}
