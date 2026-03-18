#include "definitions.h"
#include "memory.h"
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef SIMGRID
#include <smpi/smpi.h>

// =============================================================================
// SMPI Calibration Model
// =============================================================================
// Corrects for systematic differences between SMPI simulation and real MPI.
// Uses physically-motivated scaling: overhead ∝ N² (surface), compute ∝ N³ (volume)

#define BASE_FLOPS_PER_SAMPLE 35000ULL
#define BASE_FLOPS_PER_ITERATION 5000000000ULL
#define FLOPS_PER_HALO_ELEMENT 0ULL
#define SMPI_OVERHEAD_FLOPS_PER_SURFACE_ELEMENT 850000.0

#define SIMGRID_MAX_HALO_SIZE (4 * 1024 * 1024)
#define SIMGRID_NUM_BUFFERS (NEIGHBOURHOOD * 2 * 2)

static float *simgrid_buffer_pool[SIMGRID_NUM_BUFFERS];
static int simgrid_buffer_next = 0;
static int simgrid_buffers_initialized = 0;

// Approximate log10(x) without libm
static inline double approx_log10(size_t x) {
  double result = 0.0;
  while (x >= 10) { result += 1.0; x /= 10; }
  result += (double)(x - 1) / 9.0;
  return result;
}

// Approximate exp(x) using Taylor series
static inline double approx_exp(double x) {
  if (x > 10) return 22026.0;
  if (x < -10) return 0.000045;
  double result = 1.0, term = 1.0;
  for (int i = 1; i <= 10; i++) {
    term *= x / i;
    result += term;
  }
  return result;
}

// FLOPS scaling factor based on problem size (piecewise linear in log-space)
static inline double get_size_scaling_factor(size_t compute_count) {
  if (compute_count == 0) return 1.0;
  double log_count = approx_log10(compute_count);

  // Calibration points: log10(count) -> correction
  // 4.14 (size 32) -> 0.966, 5.24 (size 64) -> 0.809
  // 6.24 (size 128) -> 1.65, 7.18 (size 256) -> 0.25
  if (log_count < 4.14) {
    return 0.966 + (4.14 - log_count) * 0.01;
  } else if (log_count < 5.24) {
    double t = (log_count - 4.14) / 1.1;
    return 0.966 * (1.0 - t) + 0.809 * t;
  } else if (log_count < 6.24) {
    double t = (log_count - 5.24) / 1.0;
    return 0.809 * (1.0 - t) + 1.55 * t;
  } else if (log_count < 7.18) {
    double t = (log_count - 6.24) / 0.94;
    return 1.65 * (1.0 - t) + 0.25 * t;
  } else {
    double correction = 0.35 * (1.0 + 0.3 * (log_count - 7.18));
    return (correction < 0.2) ? 0.2 : (correction > 1.5) ? 1.5 : correction;
  }
}

// Overhead compensation (negative FLOPS) based on surface area with logistic scaling
static inline long long calc_overhead_compensation(size_t *sizes) {
  size_t face_size = (sizes[0] - 2*STENCIL) * (sizes[1] - 2*STENCIL);
  size_t surface_elements = 6 * face_size * STENCIL;

  // Logistic scaling: minimal at small sizes, full at large sizes
  double log_face = approx_log10(face_size);
  double scale = 1.0 / (1.0 + approx_exp(-5.0 * (log_face - 4.3)));

  return (long long)(-1.0 * surface_elements * SMPI_OVERHEAD_FLOPS_PER_SURFACE_ELEMENT * scale);
}

// Buffer pool management
static void simgrid_init_buffer_pool(void) {
  if (simgrid_buffers_initialized) return;
  for (int i = 0; i < SIMGRID_NUM_BUFFERS; i++) {
    simgrid_buffer_pool[i] = (float *)SMPI_SHARED_MALLOC(SIMGRID_MAX_HALO_SIZE * sizeof(float));
  }
  simgrid_buffers_initialized = 1;
}

static float *simgrid_get_buffer(void) {
  if (!simgrid_buffers_initialized) simgrid_init_buffer_pool();
  float *buf = simgrid_buffer_pool[simgrid_buffer_next];
  simgrid_buffer_next = (simgrid_buffer_next + 1) % SIMGRID_NUM_BUFFERS;
  return buf;
}

static void simgrid_reset_buffer_pool(void) {
  simgrid_buffer_next = 0;
}

#endif

#include "calculate_source.h"
#include "coordinator.h"
#include "device_data.h"
#include "indexing.h"
#include "log.h"
#include "propagate.h"
#include "worker.h"

void dc_worker_receive_data(dc_process_t *process, MPI_Comm comm) {
  MPI_Recv(&process->source_index, 1, MPI_INT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(&process->iterations, 1, MPI_UINT32_T, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);

  size_t count = dc_compute_count_from_sizes(process->sizes);
  process->pp = (float *)shared_malloc(count * sizeof(float));
  if (process->pp == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for pp in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->pc = (float *)shared_malloc(count * sizeof(float));
  if (process->pc == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for pc in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->qp = (float *)shared_malloc(count * sizeof(float));
  if (process->qp == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for qp in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->qc = (float *)shared_malloc(count * sizeof(float));
  if (process->qc == NULL) {
    dc_log_error(
        process->rank,
        "OOM: could not allocate memory for qc in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }

  MPI_Recv(process->pp, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->pc, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->qp, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->qc, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);

  process->precomp_vars.ch1dxx = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dxx == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dxx in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dyy = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dyy == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dyy in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dzz = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dzz == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dzz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dxy = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dxy == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dxy in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dyz = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dyz == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dyz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.ch1dxz = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.ch1dxz == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for "
                 "precomp_vars.ch1dxz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2px = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.v2px == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2px in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2pz = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.v2pz == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2pz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2sz = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.v2sz == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2sz in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }
  process->precomp_vars.v2pn = (float *)shared_malloc(count * sizeof(float));
  if (process->precomp_vars.v2pn == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for "
                                "precomp_vars.v2pn in dc_worker_receive_data");
    MPI_Finalize();
    exit(1);
  }

  MPI_Recv(process->precomp_vars.ch1dxx, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dyy, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dzz, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dxy, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dyz, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.ch1dxz, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2px, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2pz, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2sz, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
  MPI_Recv(process->precomp_vars.v2pn, count, MPI_FLOAT, COORDINATOR, 0, comm,
           MPI_STATUS_IGNORE);
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

#ifdef SIMGRID
    // Use pre-allocated buffer pool to avoid allocation overhead
    float *send_buffer = simgrid_get_buffer();
    // Model halo extraction as FLOPS proportional to data size
    SMPI_SAMPLE_FLOPS(data_size * FLOPS_PER_HALO_ELEMENT);
    reqs.buffers_to_free[reqs.count] = NULL;  // Don't free pooled buffer
#else
    float *send_buffer = shared_malloc(data_size * sizeof(float));
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
#endif
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
#ifdef SIMGRID
    // Use pre-allocated buffer pool to avoid allocation overhead
    result.halo_data[face_index] = simgrid_get_buffer();
#else
    result.halo_data[face_index] =
        shared_malloc(recv_data_size * sizeof(float));
    if (result.halo_data[face_index] == NULL) {
      dc_log_error(process.rank, "OOM: could not allocate memory for "
                                 "halo_data[face_index] in dc_receive_halos");
      MPI_Finalize();
      exit(1);
    }
#endif

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
        if (i < dimension) {
          start_coords[i] = 2 * radius;
          end_coords[i] = process->sizes[i] - 2 * radius;
        } else if (i == dimension) {
          start_coords[i] =
              (displacement[i] > 0) ? (process->sizes[i] - 2 * radius) : radius;
          end_coords[i] =
              (displacement[i] < 0) ? 2 * radius : process->sizes[i] - radius;
        } else {
          start_coords[i] = radius;
          end_coords[i] = process->sizes[i] - radius;
        }
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

double dc_worker_process(dc_process_t *process, MPI_Comm comm) {
  worker_requests_t all_send_requests;
  all_send_requests.buffers_to_free = NULL;
  all_send_requests.requests = NULL;
  all_send_requests.count = 0;

  dc_log_info(process->rank, "Starting %u iterations with sizes %d %d %d",
              process->iterations, process->sizes[0], process->sizes[1],
              process->sizes[2]);

  dc_device_data *data = dc_device_data_init(process);

#ifdef SIMGRID
  size_t total_compute_count = 1;
  size_t interior_compute_count = 1;
  for (int i = 0; i < 3; i++) {
    total_compute_count *= (process->sizes[i] - 2 * STENCIL);
    if (process->sizes[i] > 4 * STENCIL) {
      interior_compute_count *= (process->sizes[i] - 4 * STENCIL);
    } else {
      interior_compute_count = 0;
    }
  }
  size_t boundary_compute_count = total_compute_count - interior_compute_count;

  // Apply size-dependent scaling to match real execution times
  double scale = get_size_scaling_factor(total_compute_count);
  unsigned long long FLOPS_PER_SAMPLE = (unsigned long long)(BASE_FLOPS_PER_SAMPLE * scale);
  unsigned long long FLOPS_PER_ITERATION = (unsigned long long)(BASE_FLOPS_PER_ITERATION * scale);

  // Calculate overhead compensation (negative FLOPS to subtract SMPI overhead)
  // Distributed across 2 SMPI_SAMPLE_FLOPS calls per iteration
  long long overhead_compensation = calc_overhead_compensation(process->sizes) / 2;
#endif

  double start_time = MPI_Wtime();

  for (unsigned int i = 0; i < process->iterations; i++) {
#ifdef SIMGRID
    simgrid_reset_buffer_pool();  // Reuse pre-allocated buffers each iteration
#endif
    if (process->source_index != -1) {
      float source = dc_calculate_source(process->dt, i);
      dc_device_add_source(data, process->source_index, source);
    }

    worker_halos_t new_pp_halos = dc_receive_halos(*process, comm, PP_TAG);
    worker_halos_t new_qp_halos = dc_receive_halos(*process, comm, QP_TAG);

#ifdef SIMGRID
    // Add overhead compensation to subtract SMPI tracking overhead
    long long boundary_flops = (long long)(boundary_compute_count * FLOPS_PER_SAMPLE +
                                           FLOPS_PER_ITERATION / 2);
    long long compensated_flops = boundary_flops + overhead_compensation;
    if (compensated_flops < 0) compensated_flops = 0;  // Clamp to non-negative
    SMPI_SAMPLE_FLOPS((double)compensated_flops);
#else
    dc_compute_boundaries(process, data);
#endif
    dc_send_halo_to_neighbours(*process, comm, PP_TAG, data, data->pp,
                               &all_send_requests);
    dc_send_halo_to_neighbours(*process, comm, QP_TAG, data, data->qp,
                               &all_send_requests);
#ifdef SIMGRID
    // Add overhead compensation to subtract SMPI tracking overhead
    long long interior_flops = (long long)(interior_compute_count * FLOPS_PER_SAMPLE +
                                           FLOPS_PER_ITERATION / 2);
    long long compensated_interior = interior_flops + overhead_compensation;
    if (compensated_interior < 0) compensated_interior = 0;  // Clamp to non-negative
    SMPI_SAMPLE_FLOPS((double)compensated_interior);
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
  double msamples_per_s = msamples / elapsed;
  return msamples_per_s;
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
#ifndef SIMGRID
    // Only free in non-SIMGRID mode; pooled buffers are reused
    for (size_t i = 0; i < NEIGHBOURHOOD; i++) {
      if (halos->halo_data[i] != NULL) {
        free(halos->halo_data[i]);
      }
    }
#endif
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

#ifdef SIMGRID
        // Model halo insertion as FLOPS proportional to data size
        size_t halo_data_size = (recv_ends[0] - recv_starts[0]) *
                                (recv_ends[1] - recv_starts[1]) *
                                (recv_ends[2] - recv_starts[2]);
        SMPI_SAMPLE_FLOPS(halo_data_size * FLOPS_PER_HALO_ELEMENT);
#else
        dc_device_insert_halo_face(data, halo_buffer, recv_starts, recv_ends,
                                   process->sizes, to_array);
#endif
      }
    }
  }
}
