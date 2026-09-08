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

  // Init the exchange first: the CUDA-aware backend aborts here, before any
  // device memory is touched, if MPI lacks CUDA support.
  dc_halo_exchange_t exchange;
  dc_halo_exchange_init(process, &exchange);

  dc_device_data *data = dc_device_data_init(process);

  // pp and qp are exchanged independently, so up to 2 * NEIGHBOURHOOD
  // outstanding requests per direction.
  MPI_Request recv_requests[2 * NEIGHBOURHOOD];
  MPI_Request send_requests[2 * NEIGHBOURHOOD];

  double start_time = MPI_Wtime();

  int count = 0;
  int stopped = 0;
  double average = -1;

  for (unsigned int i = 0; i < process->iterations; i++) {
    if (process->source_index >= 0) {
      float source = dc_calculate_source(process->dt, i);
      dc_device_add_source(data, process->source_index, source);
    }

    size_t recv_count = 0;
    dc_post_halo_recvs(process, &exchange, comm, PP_TAG, 0, data->pp,
                       recv_requests, &recv_count);
    dc_post_halo_recvs(process, &exchange, comm, QP_TAG, 1, data->qp,
                       recv_requests, &recv_count);

#ifdef SIMGRID
    sampled_computation(&average, &count, &stopped, process, data,
                        dc_compute_boundaries);
#else
    dc_compute_boundaries(process, data);
#endif

    size_t send_count = 0;
    dc_post_halo_sends(process, &exchange, comm, PP_TAG, 0, data, data->pp,
                       send_requests, &send_count);
    dc_post_halo_sends(process, &exchange, comm, QP_TAG, 1, data, data->qp,
                       send_requests, &send_count);

#ifdef SIMGRID
    sampled_computation(&average, &count, &stopped, process, data,
                        dc_compute_interior);
#else
    dc_compute_interior(process, data);
#endif

    MPI_Waitall(recv_count, recv_requests, MPI_STATUSES_IGNORE);

    dc_finish_halo_recvs(process, &exchange, 0, data, data->pp);
    dc_finish_halo_recvs(process, &exchange, 1, data, data->qp);

    dc_device_swap_arrays(data);

    MPI_Waitall(send_count, send_requests, MPI_STATUSES_IGNORE);
  }

  dc_device_data_get_results(process, data);
  dc_device_data_free(data);
  dc_halo_exchange_free(&exchange);

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

void dc_worker_free(dc_process_t process) {
  free(process.pp);
  free(process.pc);
  free(process.qp);
  free(process.qc);

  free(process.hostnames);
  process.hostnames = NULL;
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
