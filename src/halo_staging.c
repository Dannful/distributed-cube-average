#include "dc_process.h"
#include "definitions.h"
#include "device_data.h"
#include "log.h"
#include "worker.h"
#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>


static void dc_face_range(int direction, size_t size, size_t radius,
                          int is_send, size_t *start, size_t *end) {
  if (direction < 0) {
    *start = is_send ? radius : 0;
    *end = is_send ? 2 * radius : radius;
  } else if (direction > 0) {
    *start = is_send ? size - 2 * radius : size - radius;
    *end = is_send ? size - radius : size;
  } else {
    *start = radius;
    *end = size - radius;
  }
}

static void dc_face_box(size_t face_index, const size_t sizes[DIMENSIONS],
                        size_t radius, int is_send, size_t start[DIMENSIONS],
                        size_t end[DIMENSIONS]) {
  const int displacement[DIMENSIONS] = {
      (int)(face_index % 3) - 1,
      (int)((face_index % 9) / 3) - 1,
      (int)(face_index / 9) - 1,
  };
  for (unsigned int d = 0; d < DIMENSIONS; d++) {
    dc_face_range(displacement[d], sizes[d], radius, is_send, &start[d],
                  &end[d]);
  }
}

static float *dc_alloc_face_buffer(int rank, size_t count) {
  float *buffer = malloc(count * sizeof(float));
  if (buffer == NULL) {
    dc_log_error(rank, "OOM: could not allocate halo staging buffer");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return buffer;
}

void dc_halo_exchange_init(const dc_process_t *process,
                           dc_halo_exchange_t *exchange) {
  const size_t radius = STENCIL;
  exchange->neighbour_count = 0;

  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    exchange->face_counts[face_index] = 0;
    for (int field = 0; field < 2; field++) {
      exchange->send_buffers[field][face_index] = NULL;
      exchange->recv_buffers[field][face_index] = NULL;
    }

    if (process->neighbours[face_index] == MPI_PROC_NULL) {
      continue;
    }

    size_t start[DIMENSIONS], end[DIMENSIONS];
    dc_face_box(face_index, process->sizes, radius, 1, start, end);
    size_t count = (end[0] - start[0]) * (end[1] - start[1]) *
                   (end[2] - start[2]);
    exchange->face_counts[face_index] = count;

    for (int field = 0; field < 2; field++) {
      exchange->send_buffers[field][face_index] =
          dc_alloc_face_buffer(process->rank, count);
      exchange->recv_buffers[field][face_index] =
          dc_alloc_face_buffer(process->rank, count);
    }

    exchange->neighbour_count++;
  }
}

void dc_halo_exchange_free(dc_halo_exchange_t *exchange) {
  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    for (int field = 0; field < 2; field++) {
      free(exchange->send_buffers[field][face_index]);
      free(exchange->recv_buffers[field][face_index]);
      exchange->send_buffers[field][face_index] = NULL;
      exchange->recv_buffers[field][face_index] = NULL;
    }
    exchange->face_counts[face_index] = 0;
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
    MPI_Irecv(exchange->recv_buffers[field][face_index],
              exchange->face_counts[face_index], MPI_FLOAT, neighbour_rank, tag,
              comm, &requests[(*count)++]);
  }
}

void dc_post_halo_sends(const dc_process_t *process,
                        const dc_halo_exchange_t *exchange, MPI_Comm comm,
                        int tag, int field, dc_device_data *data, float *array,
                        MPI_Request *requests, size_t *count) {
  const size_t radius = STENCIL;
  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    const int neighbour_rank = process->neighbours[face_index];
    if (neighbour_rank == MPI_PROC_NULL) {
      continue;
    }
    size_t start[DIMENSIONS], end[DIMENSIONS];
    dc_face_box(face_index, process->sizes, radius, 1, start, end);
    dc_device_extract_halo_face(data, exchange->send_buffers[field][face_index],
                                start, end, process->sizes, array);
    MPI_Isend(exchange->send_buffers[field][face_index],
              exchange->face_counts[face_index], MPI_FLOAT, neighbour_rank, tag,
              comm, &requests[(*count)++]);
  }
}

void dc_finish_halo_recvs(const dc_process_t *process,
                          const dc_halo_exchange_t *exchange, int field,
                          dc_device_data *data, float *array) {
  const size_t radius = STENCIL;
  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    if (process->neighbours[face_index] == MPI_PROC_NULL) {
      continue;
    }
    size_t start[DIMENSIONS], end[DIMENSIONS];
    dc_face_box(face_index, process->sizes, radius, 0, start, end);
    dc_device_insert_halo_face(data, exchange->recv_buffers[field][face_index],
                               start, end, process->sizes, array);
  }
}
