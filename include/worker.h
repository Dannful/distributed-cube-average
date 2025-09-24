#pragma once

#include "mpi.h"
#include "setup.h"
#include <stddef.h>

#define STENCIL 4
#define PP_TAG 3
#define QP_TAG 6

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
} worker_halos_t;

void dc_extract_coordinates(size_t *position_x, size_t *position_y,
                            size_t *position_z, size_t size_x, size_t size_y,
                            size_t size_z, int index);
unsigned int dc_get_index_for_coordinates(size_t position_x, size_t position_y,
                                          size_t position_z, size_t size_x,
                                          size_t size_y, size_t size_z);

void dc_worker_receive_data(dc_process_t *process);
void dc_worker_process(dc_process_t *process);
void dc_worker_free(dc_process_t process);

void dc_send_halo_to_neighbours(dc_process_t process, int tag, float *from,
                                worker_requests_t *requests);
worker_halos_t dc_receive_halos(dc_process_t process, int tag);
void dc_send_data_to_coordinator(dc_process_t process);

void dc_compute_boundaries(const dc_process_t *process, float *pp_copy, float *qp_copy);
void dc_compute_interior(const dc_process_t *process);

void dc_free_worker_halos(worker_halos_t *halos);
void dc_free_worker_requests(worker_requests_t *requests);
void dc_concatenate_worker_requests(worker_requests_t *target,
                                    worker_requests_t *source);

size_t dc_compute_count_from_sizes(size_t sizes[DIMENSIONS]);

void dc_worker_swap_arrays(dc_process_t *process);

void dc_worker_insert_halos(const dc_process_t *process,
                            const worker_halos_t *halos, float *data);
