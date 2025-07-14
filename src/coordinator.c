#include "../include/coordinator.h"
#include "../include/utils.h"
#include "../include/worker.h"
#include <math.h>
#include <string.h>

void initialize_cube(float *cube, size_t size_x, size_t size_y, size_t size_z) {
  for (float *current = cube; current < cube + size_x * size_y * size_z; current++) {
    *current = current - cube;
  }
}

unsigned int min(unsigned int a, unsigned int b) {
  return (a < b) ? a : b;
}

void partition_cube(float **workers, size_t **worker_indices, size_t *worker_count, size_t *partition_size_x, size_t *partition_size_y, size_t *partition_size_z, unsigned int num_workers, float *cube, size_t size_x, size_t size_y, size_t size_z) {
  unsigned int partition_count_x = greatest_common_divisor(size_x, num_workers);
  unsigned int partition_count_y = greatest_common_divisor(size_y, num_workers);
  unsigned int partition_count_z = greatest_common_divisor(size_z, num_workers);
  *partition_size_x = ceil((float) size_x / partition_count_x);
  *partition_size_y = ceil((float) size_y / partition_count_y);
  *partition_size_z = ceil((float) size_z / partition_count_z);
  size_t total_partition_size = *partition_size_x * *partition_size_y * *partition_size_z;
  for(size_t i = 0; i < size_x * size_y * size_z; i++) {
    size_t position_x, position_y, position_z;
    extract_coordinates(&position_x, &position_y, &position_z, size_x, size_y, size_z, i);
    unsigned char reversed = position_z / size_z % 2;
    unsigned int worker = get_assigned_worker(position_x, position_y, position_z, size_x, size_y, size_y, num_workers);
    if(worker_count[worker] == 0) {
      workers[worker] = malloc(total_partition_size * sizeof(float));
      worker_indices[worker] = malloc(total_partition_size * sizeof(size_t));
    } else {
      if(worker_count[worker] % total_partition_size == 0) {
        size_t new_size = (size_t) (((float) worker_count[worker] / total_partition_size + 1) * total_partition_size * sizeof(float));
        size_t new_size_indices = (size_t) (((float) worker_count[worker] / total_partition_size + 1) * total_partition_size * sizeof(size_t));
        workers[worker] = realloc(workers[worker], new_size);
        worker_indices[worker] = realloc(worker_indices[worker], new_size_indices);
      }
    }
    workers[worker][worker_count[worker]] = cube[i];
    worker_indices[worker][worker_count[worker]] = i;
    worker_count[worker] += 1;
    if(worker_count[worker] % total_partition_size == total_partition_size - 1 && position_x == size_x - 1) {
      workers[worker] = realloc(workers[worker], worker_count[worker] * sizeof(float));
      worker_indices[worker] = realloc(worker_indices[worker], worker_count[worker] * sizeof(size_t));
    }
  }
}
