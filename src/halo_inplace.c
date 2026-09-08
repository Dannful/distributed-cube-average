#include "dc_process.h"
#include "definitions.h"
#include "log.h"
#include "worker.h"
#include <mpi.h>
#include <stddef.h>


#ifdef CUDA_AWARE
#include <mpi-ext.h>
#endif

static void dc_face_range(int direction, size_t size, size_t radius,
                          int is_send, int *start, int *subsize) {
  if (direction < 0) {
    *start = is_send ? (int)radius : 0;
    *subsize = (int)radius;
  } else if (direction > 0) {
    *start = is_send ? (int)(size - 2 * radius) : (int)(size - radius);
    *subsize = (int)radius;
  } else {
    *start = (int)radius;
    *subsize = (int)(size - 2 * radius);
  }
}

static void dc_require_cuda_aware(const dc_process_t *process) {
#ifdef CUDA_AWARE
  if (!MPIX_Query_cuda_support()) {
    dc_log_error(process->rank,
                 "CUDA-aware MPI is required by this build but the linked MPI "
                 "library does not support it; rebuild with the plain cuda "
                 "backend for host-staged halos");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#endif
}

void dc_halo_exchange_init(const dc_process_t *process,
                           dc_halo_exchange_t *exchange) {
  dc_require_cuda_aware(process);

  const size_t radius = STENCIL;
  // x fastest, z slowest (indexing.h), so MPI_ORDER_C runs over {z, y, x}.
  const int array_sizes[DIMENSIONS] = {
      (int)process->sizes[2], (int)process->sizes[1], (int)process->sizes[0]};
  exchange->neighbour_count = 0;

  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    exchange->send_types[face_index] = MPI_DATATYPE_NULL;
    exchange->recv_types[face_index] = MPI_DATATYPE_NULL;

    if (process->neighbours[face_index] == MPI_PROC_NULL) {
      continue;
    }

    const int displacement[DIMENSIONS] = {
        (int)(face_index % 3) - 1,
        (int)((face_index % 9) / 3) - 1,
        (int)(face_index / 9) - 1,
    };

    int send_start[DIMENSIONS], send_subsize[DIMENSIONS];
    int recv_start[DIMENSIONS], recv_subsize[DIMENSIONS];
    for (unsigned int d = 0; d < DIMENSIONS; d++) {
      dc_face_range(displacement[d], process->sizes[d], radius, 1,
                    &send_start[d], &send_subsize[d]);
      dc_face_range(displacement[d], process->sizes[d], radius, 0,
                    &recv_start[d], &recv_subsize[d]);
    }

    // Reorder to the {z, y, x} axis order expected by MPI_ORDER_C.
    const int send_subsize_c[DIMENSIONS] = {send_subsize[2], send_subsize[1],
                                            send_subsize[0]};
    const int send_start_c[DIMENSIONS] = {send_start[2], send_start[1],
                                          send_start[0]};
    const int recv_subsize_c[DIMENSIONS] = {recv_subsize[2], recv_subsize[1],
                                            recv_subsize[0]};
    const int recv_start_c[DIMENSIONS] = {recv_start[2], recv_start[1],
                                          recv_start[0]};

    MPI_Type_create_subarray(DIMENSIONS, array_sizes, send_subsize_c,
                             send_start_c, MPI_ORDER_C, MPI_FLOAT,
                             &exchange->send_types[face_index]);
    MPI_Type_commit(&exchange->send_types[face_index]);

    MPI_Type_create_subarray(DIMENSIONS, array_sizes, recv_subsize_c,
                             recv_start_c, MPI_ORDER_C, MPI_FLOAT,
                             &exchange->recv_types[face_index]);
    MPI_Type_commit(&exchange->recv_types[face_index]);

    exchange->neighbour_count++;
  }
}

void dc_halo_exchange_free(dc_halo_exchange_t *exchange) {
  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    if (exchange->send_types[face_index] != MPI_DATATYPE_NULL) {
      MPI_Type_free(&exchange->send_types[face_index]);
    }
    if (exchange->recv_types[face_index] != MPI_DATATYPE_NULL) {
      MPI_Type_free(&exchange->recv_types[face_index]);
    }
  }
  exchange->neighbour_count = 0;
}

void dc_post_halo_recvs(const dc_process_t *process,
                        const dc_halo_exchange_t *exchange, MPI_Comm comm,
                        int tag, int field, float *array,
                        MPI_Request *requests, size_t *count) {
  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    const int neighbour_rank = process->neighbours[face_index];
    if (neighbour_rank == MPI_PROC_NULL) {
      continue;
    }
    MPI_Irecv(array, 1, exchange->recv_types[face_index], neighbour_rank, tag,
              comm, &requests[(*count)++]);
  }
}

void dc_post_halo_sends(const dc_process_t *process,
                        const dc_halo_exchange_t *exchange, MPI_Comm comm,
                        int tag, int field, dc_device_data *data, float *array,
                        MPI_Request *requests, size_t *count) {
  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    const int neighbour_rank = process->neighbours[face_index];
    if (neighbour_rank == MPI_PROC_NULL) {
      continue;
    }
    MPI_Isend(array, 1, exchange->send_types[face_index], neighbour_rank, tag,
              comm, &requests[(*count)++]);
  }
}

void dc_finish_halo_recvs(const dc_process_t *process,
                          const dc_halo_exchange_t *exchange, int field,
                          dc_device_data *data, float *array) {
  // Receives land directly in the halo layer, so nothing remains to do.
}
