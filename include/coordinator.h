#include <stdlib.h>

void initialize_cube(float *cube, size_t size_x, size_t size_y, size_t size_z);
void partition_cube(float **workers, size_t *total_size_x, size_t *total_size_y, size_t *total_size_z, unsigned int num_workers, float *cube, size_t size_x, size_t size_y, size_t size_z);
