#include "bits/types/struct_timeval.h"
#include "boundary.h"
#include "dc_process.h"
#include "definitions.h"
#include "precomp.h"
#include <math.h>
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
#include "sys/time.h"
#include "worker.h"
#include <smpi/smpi.h>

#ifdef SIMGRID
#include <smpi/smpi.h>
#endif

void dc_worker_init_from_partition_info(dc_process_t *process, MPI_Comm comm) {
  dc_partition_info_t info;
  MPI_Recv(&info, sizeof(dc_partition_info_t), MPI_BYTE, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);

  process->sizes[0] = info.local_sizes[0];
  process->sizes[1] = info.local_sizes[1];
  process->sizes[2] = info.local_sizes[2];
  process->iterations = info.iterations;
  process->source_index = info.source_index;

  size_t count = dc_compute_count_from_sizes(process->sizes);
  dc_log_info(process->rank,
              "Received partition info: local %zux%zux%zu, global %zux%zux%zu",
              info.local_sizes[0], info.local_sizes[1], info.local_sizes[2],
              info.global_sizes[0], info.global_sizes[1], info.global_sizes[2]);

  process->pp = (float *)calloc(count, sizeof(float));
  process->pc = (float *)calloc(count, sizeof(float));
  process->qp = (float *)calloc(count, sizeof(float));
  process->qc = (float *)calloc(count, sizeof(float));
  if (process->pp == NULL || process->pc == NULL || process->qp == NULL ||
      process->qc == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate field arrays");
    MPI_Finalize();
    exit(1);
  }

  size_t sx = info.local_sizes[0];
  size_t sy = info.local_sizes[1];
  size_t sz = info.local_sizes[2];

  process->anisotropy_vars.vpz = (float *)malloc(count * sizeof(float));
  process->anisotropy_vars.vsv = (float *)malloc(count * sizeof(float));
  process->anisotropy_vars.epsilon = (float *)malloc(count * sizeof(float));
  process->anisotropy_vars.delta = (float *)malloc(count * sizeof(float));
  process->anisotropy_vars.phi = (float *)malloc(count * sizeof(float));
  process->anisotropy_vars.theta = (float *)malloc(count * sizeof(float));
  if (process->anisotropy_vars.vpz == NULL ||
      process->anisotropy_vars.vsv == NULL ||
      process->anisotropy_vars.epsilon == NULL ||
      process->anisotropy_vars.delta == NULL ||
      process->anisotropy_vars.phi == NULL ||
      process->anisotropy_vars.theta == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate anisotropy arrays");
    MPI_Finalize();
    exit(1);
  }

  // Initialize anisotropy with default values
  for (size_t i = 0; i < count; i++) {
    process->anisotropy_vars.vpz[i] = 3000.0f;
    process->anisotropy_vars.epsilon[i] = 0.24f;
    process->anisotropy_vars.delta[i] = 0.1f;
    process->anisotropy_vars.phi[i] = 1.0f;
    process->anisotropy_vars.theta[i] = atanf(1.0);
    if (SIGMA > MAX_SIGMA) {
      process->anisotropy_vars.vsv[i] = 0.0f;
    } else {
      process->anisotropy_vars.vsv[i] =
          process->anisotropy_vars.vpz[i] *
          sqrtf(fabsf(process->anisotropy_vars.epsilon[i] -
                      process->anisotropy_vars.delta[i]) /
                SIGMA);
    }
  }

  unsigned int seed = 0;
  randomVelocityBoundaryPartition(sx, sy, sz, // Local sizes
                                  info.global_sizes[0], info.global_sizes[1],
                                  info.global_sizes[2], // Global sizes
                                  info.start_coords[0], info.start_coords[1],
                                  info.start_coords[2], // Start coords
                                  info.problem_sizes[0], info.problem_sizes[1],
                                  info.problem_sizes[2], // Problem sizes
                                  STENCIL, info.absorption_size,
                                  process->anisotropy_vars.vpz,
                                  process->anisotropy_vars.vsv, &seed);

  process->precomp_vars.ch1dxx = (float *)malloc(count * sizeof(float));
  process->precomp_vars.ch1dyy = (float *)malloc(count * sizeof(float));
  process->precomp_vars.ch1dzz = (float *)malloc(count * sizeof(float));
  process->precomp_vars.ch1dxy = (float *)malloc(count * sizeof(float));
  process->precomp_vars.ch1dyz = (float *)malloc(count * sizeof(float));
  process->precomp_vars.ch1dxz = (float *)malloc(count * sizeof(float));
  process->precomp_vars.v2px = (float *)malloc(count * sizeof(float));
  process->precomp_vars.v2pz = (float *)malloc(count * sizeof(float));
  process->precomp_vars.v2sz = (float *)malloc(count * sizeof(float));
  process->precomp_vars.v2pn = (float *)malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dxx == NULL ||
      process->precomp_vars.ch1dyy == NULL ||
      process->precomp_vars.ch1dzz == NULL ||
      process->precomp_vars.ch1dxy == NULL ||
      process->precomp_vars.ch1dyz == NULL ||
      process->precomp_vars.ch1dxz == NULL ||
      process->precomp_vars.v2px == NULL ||
      process->precomp_vars.v2pz == NULL ||
      process->precomp_vars.v2sz == NULL ||
      process->precomp_vars.v2pn == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate precomp_vars");
    MPI_Finalize();
    exit(1);
  }

  for (size_t i = 0; i < count; i++) {
    float sinTheta = sin(process->anisotropy_vars.theta[i]);
    float cosTheta = cos(process->anisotropy_vars.theta[i]);
    float sin2Theta = sin(2.0 * process->anisotropy_vars.theta[i]);
    float sinPhi = sin(process->anisotropy_vars.phi[i]);
    float cosPhi = cos(process->anisotropy_vars.phi[i]);
    float sin2Phi = sin(2.0 * process->anisotropy_vars.phi[i]);

    process->precomp_vars.ch1dxx[i] = sinTheta * sinTheta * cosPhi * cosPhi;
    process->precomp_vars.ch1dyy[i] = sinTheta * sinTheta * sinPhi * sinPhi;
    process->precomp_vars.ch1dzz[i] = cosTheta * cosTheta;
    process->precomp_vars.ch1dxy[i] = sinTheta * sinTheta * sin2Phi;
    process->precomp_vars.ch1dyz[i] = sin2Theta * sinPhi;
    process->precomp_vars.ch1dxz[i] = sin2Theta * cosPhi;

    process->precomp_vars.v2sz[i] =
        process->anisotropy_vars.vsv[i] * process->anisotropy_vars.vsv[i];
    process->precomp_vars.v2pz[i] =
        process->anisotropy_vars.vpz[i] * process->anisotropy_vars.vpz[i];
    process->precomp_vars.v2px[i] =
        process->precomp_vars.v2pz[i] *
        (1.0 + 2.0 * process->anisotropy_vars.epsilon[i]);
    process->precomp_vars.v2pn[i] =
        process->precomp_vars.v2pz[i] *
        (1.0 + 2.0 * process->anisotropy_vars.delta[i]);
  }

  dc_log_info(process->rank, "Local initialization complete");
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
  const size_t radius = STENCIL;
  const size_t *sizes = process->sizes;

  int has_interior = (sizes[0] >= 4 * radius && sizes[1] >= 4 * radius &&
                      sizes[2] >= 4 * radius);

  if (!has_interior) {
    size_t start[DIMENSIONS] = {radius, radius, radius};
    size_t end[DIMENSIONS] = {sizes[0] - radius, sizes[1] - radius,
                              sizes[2] - radius};
    if (start[0] < end[0] && start[1] < end[1] && start[2] < end[2]) {
      dc_propagate(start, end, process->sizes, process->coordinates,
                   process->topology, data, process->dx, process->dy,
                   process->dz, process->dt);
    }
    return;
  }

  size_t start[DIMENSIONS], end[DIMENSIONS];

  for (int side = 0; side < 2; side++) {
    start[0] = (side == 0) ? radius : sizes[0] - 2 * radius;
    end[0] = (side == 0) ? 2 * radius : sizes[0] - radius;
    start[1] = radius;
    end[1] = sizes[1] - radius;
    start[2] = radius;
    end[2] = sizes[2] - radius;
    if (start[0] < end[0] && start[1] < end[1] && start[2] < end[2]) {
      dc_propagate(start, end, process->sizes, process->coordinates,
                   process->topology, data, process->dx, process->dy,
                   process->dz, process->dt);
    }
  }

  for (int side = 0; side < 2; side++) {
    start[0] = 2 * radius;
    end[0] = sizes[0] - 2 * radius;
    start[1] = (side == 0) ? radius : sizes[1] - 2 * radius;
    end[1] = (side == 0) ? 2 * radius : sizes[1] - radius;
    start[2] = radius;
    end[2] = sizes[2] - radius;
    if (start[0] < end[0] && start[1] < end[1] && start[2] < end[2]) {
      dc_propagate(start, end, process->sizes, process->coordinates,
                   process->topology, data, process->dx, process->dy,
                   process->dz, process->dt);
    }
  }

  for (int side = 0; side < 2; side++) {
    start[0] = 2 * radius;
    end[0] = sizes[0] - 2 * radius;
    start[1] = 2 * radius;
    end[1] = sizes[1] - 2 * radius;
    start[2] = (side == 0) ? radius : sizes[2] - 2 * radius;
    end[2] = (side == 0) ? 2 * radius : sizes[2] - radius;
    if (start[0] < end[0] && start[1] < end[1] && start[2] < end[2]) {
      dc_propagate(start, end, process->sizes, process->coordinates,
                   process->topology, data, process->dx, process->dy,
                   process->dz, process->dt);
    }
  }
}

void dc_compute_interior(const dc_process_t *process, dc_device_data *data) {
  const size_t radius = STENCIL;
  const size_t *sizes = process->sizes;

  if (sizes[0] < 4 * radius || sizes[1] < 4 * radius || sizes[2] < 4 * radius) {
    return;
  }

  size_t start[DIMENSIONS] = {2 * radius, 2 * radius, 2 * radius};
  size_t end[DIMENSIONS] = {sizes[0] - 2 * radius, sizes[1] - 2 * radius,
                            sizes[2] - 2 * radius};

  if (start[0] < end[0] && start[1] < end[1] && start[2] < end[2]) {
    dc_propagate(start, end, process->sizes, process->coordinates,
                 process->topology, data, process->dx, process->dy, process->dz,
                 process->dt);
  }
}

void dc_send_data_to_coordinator(dc_process_t process, MPI_Comm comm) {
  if (process.rank == COORDINATOR)
    return;
#ifdef SIMGRID
  return;
#endif
  MPI_Send(process.sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR, 0, comm);
  MPI_Send(process.pc, dc_compute_count_from_sizes(process.sizes), MPI_FLOAT,
           COORDINATOR, 0, comm);
  MPI_Send(process.qc, dc_compute_count_from_sizes(process.sizes), MPI_FLOAT,
           COORDINATOR, 0, comm);
}

double get_time_micros() {
  struct timeval time;
  gettimeofday(&time, NULL);
  return ((double)time.tv_sec * 1e6) + (double)time.tv_usec;
}

#ifdef SIMGRID
void sampled_computation(double *average, int *count, int *stopped,
                         dc_process_t *process, dc_device_data *device_data,
                         void (*computation)(const dc_process_t *,
                                             dc_device_data *)) {
  const double threshold = 0.05;
  const unsigned short min_samples = 10;

  if (*stopped) {
    smpi_execute_benched(*average / 1e6);
    return;
  }

  double start = get_time_micros();
  computation(process, device_data);
  double end = get_time_micros();
  double value = end - start;
  double new_average =
      *average == -1 ? value : (*average * *count + value) / (*count + 1);
  (*count)++;

  *stopped = *average != -1 && *count >= min_samples &&
             (fabs(*average - value) / *average) < threshold;

  *average = new_average;
}
#endif

double dc_worker_process(dc_process_t *process, MPI_Comm comm) {
  dc_log_info(process->rank, "Starting %u iterations with sizes %d %d %d",
              process->iterations, process->sizes[0], process->sizes[1],
              process->sizes[2]);

  dc_device_data *data = dc_device_data_init(process);

  worker_requests_t all_send_requests = {0};

  double start_time = MPI_Wtime();

  int count = 0;
  int stopped = 0;
  double average = -1;

  for (unsigned int i = 0; i < process->iterations; i++) {
    if (process->source_index != -1) {
      float source = dc_calculate_source(process->dt, i);
      dc_device_add_source(data, process->source_index, source);
    }

    worker_halos_t new_pp_halos = dc_receive_halos(*process, comm, PP_TAG);
    worker_halos_t new_qp_halos = dc_receive_halos(*process, comm, QP_TAG);

#ifdef SIMGRID
    sampled_computation(&average, &count, &stopped, process, data,
                        dc_compute_boundaries);
#else
    dc_compute_boundaries(process, data);
#endif

    dc_send_halo_to_neighbours(*process, comm, PP_TAG, data, data->pp,
                               &all_send_requests);
    dc_send_halo_to_neighbours(*process, comm, QP_TAG, data, data->qp,
                               &all_send_requests);

#ifdef SIMGRID
    sampled_computation(&average, &count, &stopped, process, data,
                        dc_compute_interior);
#else
    dc_compute_interior(process, data);
#endif

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

  dc_device_data_get_results(process, data);
  dc_device_data_free(data);

  double end_time = MPI_Wtime();
  double elapsed = end_time - start_time;
  size_t compute_size_x = process->sizes[0] - 2 * STENCIL;
  size_t compute_size_y = process->sizes[1] - 2 * STENCIL;
  size_t compute_size_z = process->sizes[2] - 2 * STENCIL;
  double msamples = ((double)compute_size_x * compute_size_y * compute_size_z *
                     process->iterations) /
                    1000000.0;
  return msamples / elapsed;
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

  free(process.hostnames);
  process.hostnames = NULL;
}

void dc_concatenate_worker_requests(int rank, worker_requests_t *target,
                                    worker_requests_t *source) {
  if (source == NULL)
    return;
  if (source->count == 0) {
    if (source->requests != NULL) {
      free(source->requests);
      source->requests = NULL;
    }
    if (source->buffers_to_free != NULL) {
      free(source->buffers_to_free);
      source->buffers_to_free = NULL;
    }
    return;
  }
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
