#pragma once

#include "definitions.h"
#include <stddef.h>

static inline HOST_DEVICE void
dc_extract_coordinates(size_t *position_x, size_t *position_y,
                       size_t *position_z, size_t size_x, size_t size_y,
                       size_t size_z, int index) {
  *position_x = index % size_x;
  *position_y = (index / size_x) % size_y;
  *position_z = index / (size_x * size_y);
}

static inline HOST_DEVICE unsigned int
dc_get_index_for_coordinates(size_t position_x, size_t position_y,
                             size_t position_z, size_t size_x, size_t size_y,
                             size_t size_z) {
  return position_x + position_y * size_x + position_z * size_x * size_y;
}

static inline HOST_DEVICE unsigned int
dc_get_global_coordinates(const int worker_coordinates[DIMENSIONS],
                          const size_t worker_sizes[DIMENSIONS],
                          const size_t global_sizes[DIMENSIONS],
                          const size_t local_coordinates[DIMENSIONS],
                          const int topology[DIMENSIONS]) {
  size_t local_x = local_coordinates[0] - STENCIL;
  size_t local_y = local_coordinates[1] - STENCIL;
  size_t local_z = local_coordinates[2] - STENCIL;
  size_t size_x = (global_sizes[0] - 2 * STENCIL) / topology[0];
  size_t size_y = (global_sizes[1] - 2 * STENCIL) / topology[1];
  size_t size_z = (global_sizes[2] - 2 * STENCIL) / topology[2];
  size_t global_index = dc_get_index_for_coordinates(
      STENCIL + worker_coordinates[0] * size_x + local_x,
      STENCIL + worker_coordinates[1] * size_y + local_y,
      STENCIL + worker_coordinates[2] * size_z + local_z, global_sizes[0],
      global_sizes[1], global_sizes[2]);
  return global_index;
}

static inline size_t
dc_compute_count_from_sizes(const size_t sizes[DIMENSIONS]) {
  size_t count = 1;
  for (int i = 0; i < DIMENSIONS; i++) {
    count *= sizes[i];
  }
  return count;
}