#pragma once

#include "device_data.h"
#include "mpi.h"
#include "setup.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__CUDACC__)
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

#define STENCIL 4
#define PP_TAG 3
#define QP_TAG 6
#define PC_TAG 9
#define QC_TAG 12

typedef struct {
  size_t count;
  MPI_Request *requests;
  void **buffers_to_free;
} worker_requests_t;

typedef struct {
  worker_requests_t requests;
  size_t halo_count;
  size_t *halo_sizes;
  float **halo_data;
  int (*halo_dirs)[DIMENSIONS];
} worker_halos_t;

static inline HOST_DEVICE void dc_extract_coordinates(
    size_t *position_x, size_t *position_y, size_t *position_z, size_t size_x,
    size_t size_y, size_t size_z, int index) {
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

void dc_worker_receive_data(dc_process_t *process);
void dc_worker_process(dc_process_t *process);
void dc_worker_free(dc_process_t process);

void dc_send_halo_to_neighbours(dc_process_t process, int tag,
                                dc_device_data *data, float *from,
                                worker_requests_t *requests);
worker_halos_t dc_receive_halos(dc_process_t process, int tag);
void dc_send_data_to_coordinator(dc_process_t process);

void dc_compute_boundaries(const dc_process_t *process, dc_device_data *data);
void dc_compute_interior(const dc_process_t *process, dc_device_data *data);

void dc_free_worker_halos(worker_halos_t *halos);
void dc_free_worker_requests(worker_requests_t *requests);
void dc_concatenate_worker_requests(int rank, worker_requests_t *target,
                                    worker_requests_t *source);

size_t dc_compute_count_from_sizes(size_t sizes[DIMENSIONS]);

void dc_worker_swap_arrays(dc_process_t *process);

void dc_worker_insert_halos(const dc_process_t *process,
                            const worker_halos_t *halos, dc_device_data *data,
                            float *to_array);

#ifdef __cplusplus
}
#endif