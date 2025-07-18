#pragma once

#include <stddef.h>

unsigned int greatest_common_divisor(unsigned int a, unsigned int b);
float*** unflatten_cube(float *flattened_cube, size_t size_x, size_t size_y, size_t size_z);
void free_unflattened_cube(float ***unflattened_cube, size_t size_x, size_t size_y, size_t size_z);
