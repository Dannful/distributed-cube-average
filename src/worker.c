#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "calculate_source.h"
#include "coordinator.h"
#include "device_data.h"
#include "indexing.h"
#include "log.h"
#include "propagate.h"
#include "worker.h"

void dc_worker_receive_data(dc_process_t *process, MPI_Comm comm) {
  MPI_Recv(&process->source_index, 1, MPI_INT, COORDINATOR, MPI_ANY_TAG, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(&process->iterations, 1, MPI_UINT32_T, COORDINATOR, MPI_ANY_TAG,
           comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

  size_t count = dc_compute_count_from_sizes(process->sizes);
  process->pp = (float *)malloc(count * sizeof(float));
  if (process->pp == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for pp in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->pc = (float *)malloc(count * sizeof(float));
  if (process->pc == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for pc in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->qp = (float *)malloc(count * sizeof(float));
  if (process->qp == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for qp in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->qc = (float *)malloc(count * sizeof(float));
  if (process->qc == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for qc in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }

  MPI_Recv(process->pp, count, MPI_FLOAT, COORDINATOR, MPI_ANY_TAG, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->pc, count, MPI_FLOAT, COORDINATOR, MPI_ANY_TAG, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->qp, count, MPI_FLOAT, COORDINATOR, MPI_ANY_TAG, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->qc, count, MPI_FLOAT, COORDINATOR, MPI_ANY_TAG, comm,
           MPI_STATUS_IGNORE);

  process->precomp_vars.ch1dxx = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dxx == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dxx in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dyy = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dyy == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dyy in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dzz = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dzz == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dzz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dxy = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dxy == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dxy in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dyz = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dyz == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dyz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dxz = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dxz == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dxz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2px = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.v2px == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2px in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2pz = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.v2pz == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2pz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2sz = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.v2sz == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2sz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2pn = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.v2pn == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2pn in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }

  MPI_Recv(process->precomp_vars.ch1dxx, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dyy, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dzz, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dxy, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dyz, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dxz, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2px, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2pz, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2sz, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2pn, count, MPI_FLOAT, COORDINATOR,
           MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
}

void dc_send_halo_to_neighbours(dc_process_t process, MPI_Comm comm, int tag,
                                dc_device_data *data, float *from,
                                worker_requests_t *requests) {
  worker_requests_t reqs;
  size_t radius = STENCIL;

  reqs.count = 0;
  reqs.requests = malloc(NEIGHBOURHOOD * sizeof(MPI_Request));
  if (reqs.requests == NULL) {
    dc_log_error(process.rank,
                 "OOM: could not allocate memory for reqs.requests in "
                 "dc_send_halo_to_neighbours");
    MPI_Finalize();
    exit(1);
  }
  reqs.buffers_to_free = malloc(NEIGHBOURHOOD * sizeof(void *));
  if (reqs.buffers_to_free == NULL) {
    dc_log_error(process.rank,
                 "OOM: could not allocate memory for reqs.buffers_to_free in "
                 "dc_send_halo_to_neighbours");
    MPI_Finalize();
    exit(1);
  }

  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    int neighbour_rank = process.neighbours[face_index];
    if (neighbour_rank == MPI_PROC_NULL)
      continue;
    int dz = face_index / 9;
    int dy = (face_index % 9) / 3;
    int dx = face_index % 3;

    dx -= 1;
    dy -= 1;
    dz -= 1;

    size_t send_starts[3], send_ends[3];
    if (dx == -1) {
      send_starts[0] = radius;
      send_ends[0] = 2 * radius;
    } else if (dx == 1) {
      send_starts[0] = process.sizes[0] - 2 * radius;
      send_ends[0] = process.sizes[0] - radius;
    } else {
      send_starts[0] = radius;
      send_ends[0] = process.sizes[0] - radius;
    }

    if (dy == -1) {
      send_starts[1] = radius;
      send_ends[1] = 2 * radius;
    } else if (dy == 1) {
      send_starts[1] = process.sizes[1] - 2 * radius;
      send_ends[1] = process.sizes[1] - radius;
    } else {
      send_starts[1] = radius;
      send_ends[1] = process.sizes[1] - radius;
    }

    if (dz == -1) {
      send_starts[2] = radius;
      send_ends[2] = 2 * radius;
    } else if (dz == 1) {
      send_starts[2] = process.sizes[2] - 2 * radius;
      send_ends[2] = process.sizes[2] - radius;
    } else {
      send_starts[2] = radius;
      send_ends[2] = process.sizes[2] - radius;
    }

    size_t data_size = (send_ends[0] - send_starts[0]) *
                       (send_ends[1] - send_starts[1]) *
                       (send_ends[2] - send_starts[2]);
    float *send_buffer = malloc(data_size * sizeof(float));
    if (send_buffer == NULL) {
      dc_log_error(process.rank,
                   "OOM: could not allocate memory for send_buffer in "
                   "dc_send_halo_to_neighbours");
      MPI_Finalize();
      exit(1);
    }

    dc_device_extract_halo_face(data, send_buffer, send_starts, send_ends,
                                process.sizes, from);

    reqs.buffers_to_free[reqs.count] = send_buffer;
    MPI_Isend(send_buffer, data_size, MPI_FLOAT, neighbour_rank, tag, comm,
              &reqs.requests[reqs.count]);
    reqs.count++;
  }
  dc_concatenate_worker_requests(process.rank, requests, &reqs);
}

worker_halos_t dc_receive_halos(dc_process_t process, MPI_Comm comm, int tag) {
  worker_halos_t result;
  size_t radius = STENCIL;
  result.halo_count = 0;

  result.requests.count = 0;
  result.requests.requests = malloc(NEIGHBOURHOOD * sizeof(MPI_Request));
  if (result.requests.requests == NULL) {
    dc_log_error(
        process.rank,
        "OOM: could not allocate memory for requests in dc_receive_halos");
    MPI_Finalize();
    exit(1);
  }
  result.requests.buffers_to_free = NULL;

  result.halo_sizes = calloc(NEIGHBOURHOOD, sizeof(size_t));
  if (result.halo_sizes == NULL) {
    dc_log_error(
        process.rank,
        "OOM: could not allocate memory for halo_sizes in dc_receive_halos");
    MPI_Finalize();
    exit(1);
  }
  result.halo_data = calloc(NEIGHBOURHOOD, sizeof(float *));
  if (result.halo_data == NULL) {
    dc_log_error(
        process.rank,
        "OOM: could not allocate memory for halo_data in dc_receive_halos");
    MPI_Finalize();
    exit(1);
  }

  for (size_t face_index = 0; face_index < NEIGHBOURHOOD; face_index++) {
    int neighbour_rank = process.neighbours[face_index];
    if (neighbour_rank == MPI_PROC_NULL)
      continue;
    int dz = face_index / 9;
    int dy = (face_index % 9) / 3;
    int dx = face_index % 3;
    int displacement[DIMENSIONS] = {dx - 1, dy - 1, dz - 1};

    size_t recv_data_size = 1;
    for (unsigned int i = 0; i < DIMENSIONS; i++) {
      if (displacement[i] == 0) {
        recv_data_size *= process.sizes[i] - 2 * radius;
      } else {
        recv_data_size *= radius;
      }
    }
    result.halo_sizes[face_index] = recv_data_size;
    result.halo_data[face_index] = malloc(recv_data_size * sizeof(float));
    if (result.halo_data[face_index] == NULL) {
      dc_log_error(process.rank, "OOM: could not allocate memory for "
                                 "halo_data[face_index] in dc_receive_halos");
      MPI_Finalize();
      exit(1);
    }

    MPI_Irecv(result.halo_data[face_index], recv_data_size, MPI_FLOAT,
              neighbour_rank, tag, comm,
              &result.requests.requests[result.requests.count]);

    result.halo_count++;
    result.requests.count++;
  }
  return result;
}

void dc_compute_boundaries(const dc_process_t *process, dc_device_data *data) {
  const int radius = STENCIL;

  for (unsigned int dimension = 0; dimension < DIMENSIONS; dimension++) {
    for (int direction = -1; direction <= 1; direction += 2) {
      int displacement[DIMENSIONS] = {0};
      displacement[dimension] = direction;
      size_t start_coords[DIMENSIONS], end_coords[DIMENSIONS];
      for (unsigned int i = 0; i < DIMENSIONS; i++) {
        start_coords[i] =
            (displacement[i] > 0) ? (process->sizes[i] - 2 * radius) : radius;
        end_coords[i] =
            (displacement[i] < 0) ? 2 * radius : process->sizes[i] - radius;
      }
      dc_propagate(start_coords, end_coords, process->sizes,
                   process->coordinates, process->topology, data, process->dx,
                   process->dy, process->dz, process->dt);
    }
  }
}

void dc_compute_interior(const dc_process_t *process, dc_device_data *data) {
  const int radius = STENCIL;

  size_t start_coords[DIMENSIONS] = {2 * radius, 2 * radius, 2 * radius};
  size_t end_coords[DIMENSIONS] = {process->sizes[0] - 2 * radius,
                                   process->sizes[1] - 2 * radius,
                                   process->sizes[2] - 2 * radius};

  dc_propagate(start_coords, end_coords, process->sizes, process->coordinates,
               process->topology, data, process->dx, process->dy, process->dz,
               process->dt);
}

void dc_send_data_to_coordinator(dc_process_t process, MPI_Comm comm) {
  if (process.rank == COORDINATOR)
    return;
  MPI_Send(process.sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR, 0, comm);
  MPI_Send(process.pc, dc_compute_count_from_sizes(process.sizes), MPI_FLOAT,
           COORDINATOR, 0, comm);
  MPI_Send(process.qc, dc_compute_count_from_sizes(process.sizes), MPI_FLOAT,
           COORDINATOR, 0, comm);
}

void dc_worker_process(dc_process_t *process, MPI_Comm comm) {
  worker_requests_t all_send_requests;
  all_send_requests.buffers_to_free = NULL;
  all_send_requests.requests = NULL;
  all_send_requests.count = 0;

  dc_log_info(process->rank, "Starting %u iterations with sizes %d %d %d",
              process->iterations, process->sizes[0], process->sizes[1],
              process->sizes[2]);

  dc_device_data *data = dc_device_data_init(process);

  for (unsigned int i = 0; i < process->iterations; i++) {
    if (process->source_index != -1) {
      float source = dc_calculate_source(process->dt, i);
      dc_log_info(process->rank, "Inserting source %f at %d", source,
                  process->source_index);
      dc_device_add_source(data, process->source_index, source);
    }

    dc_device_data_copy_to_device_copies(data, process->sizes);

    worker_halos_t new_pp_halos = dc_receive_halos(*process, comm, PP_TAG);
    worker_halos_t new_qp_halos = dc_receive_halos(*process, comm, QP_TAG);

    dc_compute_boundaries(process, data);
    dc_send_halo_to_neighbours(*process, comm, PP_TAG, data, data->pp,
                               &all_send_requests);
    dc_send_halo_to_neighbours(*process, comm, QP_TAG, data, data->qp,
                               &all_send_requests);
    dc_compute_interior(process, data);

    dc_concatenate_worker_requests(process->rank, &new_pp_halos.requests,
                                   &new_qp_halos.requests);
    MPI_Waitall(new_pp_halos.requests.count, new_pp_halos.requests.requests,
                MPI_STATUSES_IGNORE);
    dc_worker_insert_halos(process, &new_pp_halos, data, data->pp);
    dc_worker_insert_halos(process, &new_qp_halos, data, data->qp);
    dc_free_worker_halos(&new_pp_halos);
    dc_free_worker_halos(&new_qp_halos);
    dc_device_swap_arrays(data);
    MPI_Waitall(all_send_requests.count, all_send_requests.requests,
                MPI_STATUSES_IGNORE);
    dc_free_worker_requests(&all_send_requests);
  }

  process->pp = data->pp;
  process->pc = data->pc;
  process->qp = data->qp;
  process->qc = data->qc;

  dc_device_data_get_results(process, data);
  dc_device_data_free(data);
  dc_log_info(process->rank, "Processing complete.");
}

void dc_free_worker_requests(worker_requests_t *requests) {
  if (requests->buffers_to_free != NULL) {
    for (size_t i = 0; i < requests->count; i++) {
      if (requests->buffers_to_free[i] != NULL) {
        free(requests->buffers_to_free[i]);
      }
    }
    free(requests->buffers_to_free);
  }
  if (requests->requests != NULL) {
    free(requests->requests);
  }
  requests->requests = NULL;
  requests->buffers_to_free = NULL;
  requests->count = 0;
}

void dc_free_worker_halos(worker_halos_t *halos) {
  if (halos->halo_data != NULL) {
    for (size_t i = 0; i < NEIGHBOURHOOD; i++) {
      if (halos->halo_data[i] != NULL) {
        free(halos->halo_data[i]);
      }
    }
    free(halos->halo_data);
  }
  if (halos->halo_sizes != NULL) {
    free(halos->halo_sizes);
  }
  dc_free_worker_requests(&halos->requests);
  halos->halo_data = NULL;
  halos->halo_sizes = NULL;
  halos->halo_count = 0;
}

void dc_worker_free(dc_process_t process) {
  free(process.pp);
  free(process.pc);
  free(process.qp);
  free(process.qc);
}

void dc_concatenate_worker_requests(int rank, worker_requests_t *target,
                                    worker_requests_t *source) {
  if (source == NULL || source->count == 0)
    return;
  size_t original_target_count = target->count;
  size_t new_count = original_target_count + source->count;
  target->requests = realloc(target->requests, new_count * sizeof(MPI_Request));
  if (target->requests == NULL) {
    dc_log_error(rank, "OOM: could not allocate memory for target->requests in "
                       "dc_concatenate_worker_requests");
    MPI_Finalize();
    exit(1);
  }
  memcpy(target->requests + original_target_count, source->requests,
         source->count * sizeof(MPI_Request));
  if (source->buffers_to_free != NULL) {
    if (target->buffers_to_free == NULL) {
      target->buffers_to_free = malloc(new_count * sizeof(void *));
      if (target->buffers_to_free == NULL) {
        dc_log_error(
            rank, "OOM: could not allocate memory for target->buffers_to_free "
                  "in dc_concatenate_worker_requests");
        MPI_Finalize();
        exit(1);
      }
      memset(target->buffers_to_free, 0,
             original_target_count * sizeof(void *));
    } else {
      target->buffers_to_free =
          realloc(target->buffers_to_free, new_count * sizeof(void *));
      if (target->buffers_to_free == NULL) {
        dc_log_error(rank, "OOM: could not re-allocate memory for "
                           "target->buffers_to_free in "
                           "dc_concatenate_worker_requests");
        MPI_Finalize();
        exit(1);
      }
    }
    memcpy(target->buffers_to_free + original_target_count,
           source->buffers_to_free, source->count * sizeof(void *));
  } else if (target->buffers_to_free != NULL) {
    target->buffers_to_free =
        realloc(target->buffers_to_free, new_count * sizeof(void *));
    if (target->buffers_to_free == NULL) {
      dc_log_error(rank, "OOM: could not re-allocate memory for "
                         "target->buffers_to_free in "
                         "dc_concatenate_worker_requests");
      MPI_Finalize();
      exit(1);
    }
    memset(target->buffers_to_free + original_target_count, 0,
           source->count * sizeof(void *));
  }
  target->count = new_count;

  free(source->requests);
  free(source->buffers_to_free);
  source->requests = NULL;
  source->buffers_to_free = NULL;
  source->count = 0;
}

void dc_worker_swap_arrays(dc_process_t *process) {
  float *temp;

  temp = process->pp;
  process->pp = process->pc;
  process->pc = temp;

  temp = process->qp;
  process->qp = process->qc;
  process->qc = temp;
}

void dc_worker_insert_halos(const dc_process_t *process,
                            const worker_halos_t *halos, dc_device_data *data,
                            float *to_array) {
  const size_t radius = STENCIL;

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        if (dx == 0 && dy == 0 && dz == 0) {
          continue;
        }

        size_t face_index = 9 * (dz + 1) + 3 * (dy + 1) + dx + 1;
        if (process->neighbours[face_index] == MPI_PROC_NULL) {
          continue;
        }
        float *halo_buffer = halos->halo_data[face_index];

        if (halo_buffer == NULL) {
          continue;
        }

        size_t recv_starts[3], recv_ends[3];
        if (dx == -1) {
          recv_starts[0] = 0;
          recv_ends[0] = radius;
        } else if (dx == 1) {
          recv_starts[0] = process->sizes[0] - radius;
          recv_ends[0] = process->sizes[0];
        } else {
          recv_starts[0] = radius;
          recv_ends[0] = process->sizes[0] - radius;
        }

        if (dy == -1) {
          recv_starts[1] = 0;
          recv_ends[1] = radius;
        } else if (dy == 1) {
          recv_starts[1] = process->sizes[1] - radius;
          recv_ends[1] = process->sizes[1];
        } else {
          recv_starts[1] = radius;
          recv_ends[1] = process->sizes[1] - radius;
        }

        if (dz == -1) {
          recv_starts[2] = 0;
          recv_ends[2] = radius;
        } else if (dz == 1) {
          recv_starts[2] = process->sizes[2] - radius;
          recv_ends[2] = process->sizes[2];
        } else {
          recv_starts[2] = radius;
          recv_ends[2] = process->sizes[2] - radius;
        }

        dc_device_insert_halo_face(data, halo_buffer, recv_starts, recv_ends,
                                   process->sizes, to_array);
      }
    }
  }
}
