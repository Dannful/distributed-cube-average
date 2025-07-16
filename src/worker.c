#include <stdlib.h>

#include "../include/worker.h"
#include "../include/utils.h"

void extract_coordinates(size_t *position_x, size_t *position_y, size_t *position_z, size_t size_x, size_t size_y, size_t size_z, int index) {
  *position_x = index % size_x;
  *position_y = (index / size_x) % size_y;
  *position_z = index / (size_x * size_y);
}

unsigned int get_index_for_coordinates(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z) {
  return position_x + position_y * size_x + position_z * size_x * size_y;
}

unsigned int get_assigned_worker(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z, size_t worker_count) {
  unsigned int partition_count_x = greatest_common_divisor(size_x, worker_count);
  unsigned int partition_count_y = greatest_common_divisor(size_y, worker_count);
  unsigned int partition_count_z = greatest_common_divisor(size_z, worker_count);
  unsigned char reversed = position_z / size_z % 2;
  unsigned int worker = reversed ? worker_count - (position_x / partition_count_x + position_y / partition_count_y) % worker_count : (position_x / partition_count_x + position_y / partition_count_y) % worker_count;
  return worker;
}

void worker_process(unsigned int iterations, unsigned int stencil_size, float *data, size_t count, size_t *indices) {
  size_t *partition_sizes = calloc(sizeof(size_t), count);
  for(int i = 0; i < iterations; i++) {
    for(size_t j = 0; j < count; j++) {
      size_t position_x, position_y, position_z;
      extract_coordinates(&position_x, &position_y, &position_z, indices[0], indices[1], indices[2], j);
      unsigned int worker = get_assigned_worker(position_x, position_y, position_z, indices[0], indices[1], indices[2], indices[3]);
    }
  }
  free(partition_sizes);
}
