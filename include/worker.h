#pragma once

#include "dc_process.h"
#include "device_data.h"
#include "mpi.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PP_TAG 3
#define QP_TAG 6
#define PC_TAG 9
#define QC_TAG 12

// Persistent MPI derived datatypes describing the halo exchange for every
// neighbour. The send type selects the owned boundary cells, the recv type
// selects the outer halo layer. Both address the field arrays in place, so no
// extraction/reinsertion copies are needed. Unused faces hold MPI_DATATYPE_NULL.
typedef struct {
  MPI_Datatype send_types[NEIGHBOURHOOD];
  MPI_Datatype recv_types[NEIGHBOURHOOD];
  size_t neighbour_count;
} dc_halo_exchange_t;

void dc_worker_init_from_partition_info(dc_process_t *process, MPI_Comm comm);
double dc_worker_process(dc_process_t *process, MPI_Comm comm);
void dc_worker_free(dc_process_t process);

void dc_halo_exchange_init(const dc_process_t *process,
                           dc_halo_exchange_t *exchange);
void dc_halo_exchange_free(dc_halo_exchange_t *exchange);
void dc_post_halo_recvs(const dc_process_t *process,
                        const dc_halo_exchange_t *exchange, MPI_Comm comm,
                        int tag, float *array, MPI_Request *requests,
                        size_t *count);
void dc_post_halo_sends(const dc_process_t *process,
                        const dc_halo_exchange_t *exchange, MPI_Comm comm,
                        int tag, float *array, MPI_Request *requests,
                        size_t *count);

void dc_send_data_to_coordinator(dc_process_t process, MPI_Comm comm);

void dc_compute_boundaries(const dc_process_t *process, dc_device_data *data);
void dc_compute_interior(const dc_process_t *process, dc_device_data *data);

void dc_worker_swap_arrays(dc_process_t *process);

#ifdef __cplusplus
}
#endif
