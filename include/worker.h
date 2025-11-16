#pragma once

#include "dc_process.h"
#include "device_data.h"
#include "mpi.h"
#include "setup.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

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

void dc_worker_receive_data(dc_process_t *process, MPI_Comm comm);
void dc_worker_process(dc_process_t *process, MPI_Comm comm);
void dc_worker_free(dc_process_t process);

void dc_send_halo_to_neighbours(dc_process_t process, MPI_Comm comm, int tag,
                                dc_device_data *data, float *from,
                                worker_requests_t *requests);
worker_halos_t dc_receive_halos(dc_process_t process, MPI_Comm comm, int tag);
void dc_send_data_to_coordinator(dc_process_t process, MPI_Comm comm);

void dc_compute_boundaries(const dc_process_t *process, dc_device_data *data);
void dc_compute_interior(const dc_process_t *process, dc_device_data *data);

void dc_free_worker_halos(worker_halos_t *halos);
void dc_free_worker_requests(worker_requests_t *requests);
void dc_concatenate_worker_requests(int rank, worker_requests_t *target,
                                    worker_requests_t *source);

void dc_worker_swap_arrays(dc_process_t *process);

void dc_worker_insert_halos(const dc_process_t *process,
                            const worker_halos_t *halos, dc_device_data *data,
                            float *to_array);

#ifdef __cplusplus
}
#endif