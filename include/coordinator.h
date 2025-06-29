#include <stdlib.h>

void initialize_cube(float *cube, size_t size_x, size_t size_y, size_t size_z);
void partition_cube(float **workers, unsigned int num_workers, float *cube, size_t size_x, size_t size_y, size_t size_z);
