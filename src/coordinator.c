#include "../include/coordinator.h"
#include "../include/log.h"
#include <math.h>
#include <string.h>

void initialize_cube(float *cube, size_t size_x, size_t size_y, size_t size_z) {
  for (float *current = cube; current < cube + size_x * size_y * size_z; current++) {
    *current = current - cube;
  }
}

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
      a = (a - b) >> 1;
    } else if(b > a) {
      b = (b - a) >> 1;
    }
  }
  return a << d;
}

unsigned int min(unsigned int a, unsigned int b) {
  return (a < b) ? a : b;
}

void partition_cube(float **workers, unsigned int num_workers, float *cube, size_t size_x, size_t size_y, size_t size_z) {
  unsigned int partition_size_x = greatest_common_divisor(size_x, num_workers);
  unsigned int partition_size_y = greatest_common_divisor(size_y, num_workers);
  unsigned int partition_size_z = greatest_common_divisor(size_z, num_workers);
  size_t total_size_x = ceil((float) size_x / partition_size_x);
  size_t total_size_y = ceil((float) size_y / partition_size_y);
  size_t total_size_z = ceil((float) size_z / partition_size_z);
  for(unsigned int worker = 0; worker < num_workers; worker++) {
    workers[worker] = malloc(total_size_x * total_size_y * total_size_z * sizeof(float));
    if (workers[worker] == NULL) {
      log_error(0, "Out of memory! Exitting...", worker);
      exit(EXIT_FAILURE);
    }
  }
  for(size_t i = 0; i < size_x * size_y * size_z; i++) {
    size_t position_x = i % size_x;
    size_t position_y = (i / size_x) % size_y;
    size_t position_z = i / (size_x * size_y);
    unsigned char reversed = position_z / size_z % 2;
    unsigned int worker = reversed ? num_workers - (position_x / size_x + position_y / size_y) % num_workers : (position_x / size_x + position_y / size_y) % num_workers;
    workers[worker][i] = cube[i];
  }
}
