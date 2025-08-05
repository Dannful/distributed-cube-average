#pragma once

#include "mpi.h"
#include "setup.h"
#include <stddef.h>

typedef struct {
  size_t count;
  MPI_Request *requests;
} worker_requests_t;

typedef struct {
  worker_requests_t requests;
  size_t halo_count;
  size_t *halo_sizes;
  size_t **halo_indices;
  size_t **halo_data;
} worker_halos_t;

void dc_extract_coordinates(size_t *position_x, size_t *position_y, size_t *position_z, size_t size_x, size_t size_y, size_t size_z, int index);
unsigned int dc_get_index_for_coordinates(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z);
void dc_worker_receive_data(dc_process_t *process);
void dc_worker_process(dc_process_t process);
void dc_worker_free(dc_process_t process);
worker_halos_t dc_receive_halos(dc_process_t process);
worker_requests_t dc_send_halo_to_neighbours(dc_process_t process, float *from);
void dc_free_worker_halos(worker_halos_t *halos);
void dc_free_worker_requests(worker_requests_t *requests);
void dc_concatenate_worker_requests(worker_requests_t *target, worker_requests_t *source);
