#include <stdlib.h>
#include "../include/utils.h"

unsigned int greatest_common_divisor(unsigned int a, unsigned int b) {
  unsigned int d = 0;
    while ((a & 1) == 0 && (b & 1) == 0) {
      a >>= 1;
      b >>= 1;
      d++;
    }
    while ((a & 1) == 0) {
      a >>= 1;
    }
    while ((b & 1) == 0) {
      b >>= 1;
    }
    while(a != b) {
      if (a > b) {
        a -= b;
        while((a & 1) == 0) {
          a >>= 1;
        }
      } else if(b > a) {
        b -= a;
        while((b & 1) == 0) {
          b >>= 1;
        }
      } else {
        exit(EXIT_FAILURE);
      }
    }
    return a << d;
}

float*** unflatten_cube(float *flattened_cube, size_t size_x, size_t size_y, size_t size_z) {
  float ***result = malloc(sizeof(float**) * size_x);
  for(size_t x = 0; x < size_x; x++) {
    result[x] = malloc(sizeof(float*) * size_y);
    for(size_t y = 0; y < size_y; y++) {
      result[x][y] = malloc(sizeof(float) * size_z);
      for(size_t z = 0; z < size_z; z++) {
        result[x][y][z] = flattened_cube[x + y * size_x + z * size_x * size_y];
      }
    }
  }
  return result;
}

void free_unflattened_cube(float ***unflattened_cube, size_t size_x, size_t size_y, size_t size_z) {
  for(size_t x = 0; x < size_x; x++) {
    for(size_t y = 0; y < size_y; y++) {
      free(unflattened_cube[x][y]);
    }
    free(unflattened_cube[x]);
  }
  free(unflattened_cube);
}
