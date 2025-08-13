#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coordinator.h"
#include "log.h"
#include "setup.h"
#include "worker.h"

void dc_extract_coordinates(size_t *position_x, size_t *position_y,
                            size_t *position_z, size_t size_x, size_t size_y,
                            size_t size_z, int index) {
  *position_x = index % size_x;
  *position_y = (index / size_x) % size_y;
  *position_z = index / (size_x * size_y);
}

unsigned int dc_get_index_for_coordinates(size_t position_x, size_t position_y,
                                          size_t position_z, size_t size_x,
                                          size_t size_y, size_t size_z) {
  return position_x + position_y * size_x + position_z * size_x * size_y;
}

void dc_worker_receive_data(dc_process_t *process) {
  MPI_Recv(&process->stencil_size, 1, MPI_UINT32_T, COORDINATOR, MPI_ANY_TAG,
           process->communicator, MPI_STATUS_IGNORE);
  MPI_Recv(&process->iterations, 1, MPI_UINT32_T, COORDINATOR, MPI_ANY_TAG,
           process->communicator, MPI_STATUS_IGNORE);
  MPI_Recv(&process->count, 1, MPI_UNSIGNED_LONG, COORDINATOR, MPI_ANY_TAG,
           process->communicator, MPI_STATUS_IGNORE);
  MPI_Recv(process->sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR,
           MPI_ANY_TAG, process->communicator, MPI_STATUS_IGNORE);

  process->indices = malloc(process->count * sizeof(size_t));
  process->data = malloc(process->count * sizeof(float));

  MPI_Recv(process->indices, process->count, MPI_UNSIGNED_LONG, COORDINATOR,
           MPI_ANY_TAG, process->communicator, MPI_STATUS_IGNORE);
  MPI_Recv(process->data, process->count, MPI_FLOAT, COORDINATOR, MPI_ANY_TAG,
           process->communicator, MPI_STATUS_IGNORE);
}

worker_requests_t dc_send_halo_to_neighbours(dc_process_t process,
                                             float *from) {
  worker_requests_t reqs;
  size_t radius = (process.stencil_size - 1) / 2;

  reqs.count = 0;
  reqs.requests = malloc(6 * sizeof(MPI_Request));
  reqs.buffers_to_free = malloc(6 * sizeof(void *));

  for (int dim = 0; dim < DIMENSIONS; dim++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int neighbour_rank;
      int coordinates[DIMENSIONS];
      memcpy(coordinates, process.coordinates, sizeof(int) * DIMENSIONS);
      coordinates[dim] += dir;
      if(coordinates[dim] < 0 || coordinates[dim] >= process.topology[dim])
        continue;
      MPI_Cart_rank(process.communicator, coordinates, &neighbour_rank);

      int displacement[DIMENSIONS] = {0};
      displacement[dim] = dir;
      size_t data_size = 1;
      for (unsigned int i = 0; i < DIMENSIONS; i++) {
        data_size *= (displacement[i] == 0) ? process.sizes[i] : radius;
      }
      float *send_buffer = malloc(sizeof(float) * data_size);

      size_t start_coords[DIMENSIONS], end_coords[DIMENSIONS];
      for (unsigned int i = 0; i < DIMENSIONS; i++) {
        start_coords[i] =
            (displacement[i] > 0) ? (process.sizes[i] - radius) : 0;
        end_coords[i] = (displacement[i] < 0) ? radius : process.sizes[i];
      }
      size_t data_index = 0;
      for (size_t z = start_coords[2]; z < end_coords[2]; z++)
        for (size_t y = start_coords[1]; y < end_coords[1]; y++)
          for (size_t x = start_coords[0]; x < end_coords[0]; x++)
            send_buffer[data_index++] = from[dc_get_index_for_coordinates(
                x, y, z, process.sizes[0], process.sizes[1], process.sizes[2])];

      dc_log_info(process.rank,
                  "Sending halo of %zu element(s) to neighbour %d", data_size,
                  neighbour_rank);

      reqs.buffers_to_free[reqs.count] = send_buffer;
      MPI_Isend(send_buffer, data_size, MPI_FLOAT, neighbour_rank, 0,
                process.communicator, &reqs.requests[reqs.count]);
      reqs.count++;
    }
  }
  return reqs;
}

worker_halos_t dc_receive_halos(dc_process_t process) {
  worker_halos_t result;
  size_t radius = (process.stencil_size - 1) / 2;
  result.halo_count = 0;

  result.requests.count = 0;
  result.requests.requests = malloc(6 * sizeof(MPI_Request));
  result.requests.buffers_to_free = NULL;

  result.halo_sizes = malloc(6 * sizeof(size_t));
  result.halo_data = malloc(6 * sizeof(float *));

  for (int dim = 0; dim < DIMENSIONS; dim++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int source_rank;
      int coordinates[DIMENSIONS];
      memcpy(coordinates, process.coordinates, sizeof(int) * DIMENSIONS);
      coordinates[dim] += dir;
      if(coordinates[dim] < 0 || coordinates[dim] >= process.topology[dim])
        continue;
      MPI_Cart_rank(process.communicator, coordinates, &source_rank);

      size_t recv_data_size = 1;
      for (unsigned int i = 0; i < DIMENSIONS; i++) {
        recv_data_size *= (i == dim) ? radius : process.sizes[i];
      }

      result.halo_sizes[result.halo_count] = recv_data_size;
      result.halo_data[result.halo_count] =
          malloc(recv_data_size * sizeof(float));

      dc_log_info(process.rank,
                  "Receiving halo from worker %d...", source_rank);
      MPI_Irecv(result.halo_data[result.halo_count], recv_data_size, MPI_FLOAT,
                source_rank, 0, process.communicator,
                &result.requests.requests[result.requests.count]);

      result.halo_count++;
      result.requests.count++;
    }
  }
  return result;
}

float *dc_populate_padded_grid(const dc_process_t *process,
                               const worker_halos_t *halos,
                               const float *input_data) {
  const size_t radius = (process->stencil_size - 1) / 2;
  const size_t padded_size_x = process->sizes[0] + 2 * radius;
  const size_t padded_size_y = process->sizes[1] + 2 * radius;
  const size_t padded_size_z = process->sizes[2] + 2 * radius;
  const size_t padded_grid_count =
      padded_size_x * padded_size_y * padded_size_z;
  float *padded_grid = calloc(padded_grid_count, sizeof(float));

  for (size_t x = 0; x < process->sizes[0]; x++) {
    for (size_t y = 0; y < process->sizes[1]; y++) {
      for (size_t z = 0; z < process->sizes[2]; z++) {
        size_t local_idx = dc_get_index_for_coordinates(
            x, y, z, process->sizes[0], process->sizes[1], process->sizes[2]);
        size_t padded_idx = dc_get_index_for_coordinates(
            x + radius, y + radius, z + radius, padded_size_x, padded_size_y,
            padded_size_z);
        padded_grid[padded_idx] = input_data[local_idx];
      }
    }
  }

  int halo_idx_counter = 0;
  for (int dim = 0; dim < DIMENSIONS; ++dim) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int source_rank;
      int coordinates[DIMENSIONS];
      memcpy(coordinates, process->coordinates, sizeof(int) * DIMENSIONS);
      coordinates[dim] += dir;
      if(coordinates[dim] < 0 || coordinates[dim] >= process->topology[dim])
        continue;
      MPI_Cart_rank(process->communicator, coordinates, &source_rank);

      size_t start_x = (dir == -1) ? 0 : process->sizes[0] + radius;
      size_t start_y = (dir == -1) ? 0 : process->sizes[1] + radius;
      size_t start_z = (dir == -1) ? 0 : process->sizes[2] + radius;

      if (dim == 0) {
        start_y = radius;
        start_z = radius;
      }
      if (dim == 1) {
        start_x = radius;
        start_z = radius;
      }
      if (dim == 2) {
        start_x = radius;
        start_y = radius;
      }

      size_t plane_size_x = (dim == 0) ? radius : process->sizes[0];
      size_t plane_size_y = (dim == 1) ? radius : process->sizes[1];
      size_t plane_size_z = (dim == 2) ? radius : process->sizes[2];

      size_t current_halo_data_idx = 0;
      for (size_t x = 0; x < plane_size_x; x++) {
        for (size_t y = 0; y < plane_size_y; y++) {
          for (size_t z = 0; z < plane_size_z; z++) {
            size_t padded_idx = dc_get_index_for_coordinates(
                start_x + x, start_y + y, start_z + z, padded_size_x,
                padded_size_y, padded_size_z);
            if (halo_idx_counter < halos->halo_count &&
                current_halo_data_idx < halos->halo_sizes[halo_idx_counter]) {
              padded_grid[padded_idx] =
                  halos->halo_data[halo_idx_counter][current_halo_data_idx++];
            }
          }
        }
      }
      halo_idx_counter++;
    }
  }
  return padded_grid;
}

void dc_compute_boundaries(const dc_process_t *process, float *output_data,
                           const float *input_data,
                           const worker_halos_t *halos) {
  float *padded_grid = dc_populate_padded_grid(process, halos, input_data);

  const size_t radius = (process->stencil_size - 1) / 2;
  const size_t size_x = process->sizes[0];
  const size_t size_y = process->sizes[1];
  const size_t size_z = process->sizes[2];
  const size_t padded_size_x = size_x + 2 * radius;
  const size_t padded_size_y = size_y + 2 * radius;
  const size_t padded_size_z = size_z + 2 * radius;

  for (size_t x = 0; x < size_x; x++) {
    for (size_t y = 0; y < size_y; y++) {
      for (size_t z = 0; z < size_z; z++) {
        if (x < radius || x >= size_x - radius || y < radius ||
            y >= size_y - radius || z < radius || z >= size_z - radius) {
          const size_t px = x + radius, py = y + radius, pz = z + radius;
          double sum = padded_grid[dc_get_index_for_coordinates(
              px, py, pz, padded_size_x, padded_size_y, padded_size_z)];
          sum += padded_grid[dc_get_index_for_coordinates(
              px - 1, py, pz, padded_size_x, padded_size_y, padded_size_z)];
          sum += padded_grid[dc_get_index_for_coordinates(
              px + 1, py, pz, padded_size_x, padded_size_y, padded_size_z)];
          sum += padded_grid[dc_get_index_for_coordinates(
              px, py - 1, pz, padded_size_x, padded_size_y, padded_size_z)];
          sum += padded_grid[dc_get_index_for_coordinates(
              px, py + 1, pz, padded_size_x, padded_size_y, padded_size_z)];
          sum += padded_grid[dc_get_index_for_coordinates(
              px, py, pz - 1, padded_size_x, padded_size_y, padded_size_z)];
          sum += padded_grid[dc_get_index_for_coordinates(
              px, py, pz + 1, padded_size_x, padded_size_y, padded_size_z)];
          output_data[dc_get_index_for_coordinates(x, y, z, size_x, size_y,
                                                   size_z)] = sum / ((process->stencil_size - 1) * DIMENSIONS + 1);
        }
      }
    }
  }
  free(padded_grid);
}

void dc_compute_interior(const dc_process_t *process, float *output_data,
                         const float *input_data) {
  const size_t radius = (process->stencil_size - 1) / 2;
  const size_t size_x = process->sizes[0];
  const size_t size_y = process->sizes[1];
  const size_t size_z = process->sizes[2];

  for (size_t z = radius; z < size_z - radius; ++z) {
    for (size_t y = radius; y < size_y - radius; ++y) {
      for (size_t x = radius; x < size_x - radius; ++x) {
        double sum = input_data[dc_get_index_for_coordinates(x, y, z, size_x,
                                                             size_y, size_z)];
        sum += input_data[dc_get_index_for_coordinates(x - 1, y, z, size_x,
                                                       size_y, size_z)];
        sum += input_data[dc_get_index_for_coordinates(x + 1, y, z, size_x,
                                                       size_y, size_z)];
        sum += input_data[dc_get_index_for_coordinates(x, y - 1, z, size_x,
                                                       size_y, size_z)];
        sum += input_data[dc_get_index_for_coordinates(x, y + 1, z, size_x,
                                                       size_y, size_z)];
        sum += input_data[dc_get_index_for_coordinates(x, y, z - 1, size_x,
                                                       size_y, size_z)];
        sum += input_data[dc_get_index_for_coordinates(x, y, z + 1, size_x,
                                                       size_y, size_z)];
        output_data[dc_get_index_for_coordinates(x, y, z, size_x, size_y,
                                                 size_z)] = sum / ((process->stencil_size - 1) * DIMENSIONS + 1);
      }
    }
  }
}

void dc_worker_process(dc_process_t process) {
  float *current_data = malloc(process.count * sizeof(float));
  memcpy(current_data, process.data, process.count * sizeof(float));
  float *next_data = malloc(process.count * sizeof(float));
  float *temp_ptr;

  worker_halos_t current_halos, future_halos;
  dc_log_info(process.rank,
              "Starting %u iterations of stencil computation...",
              process.iterations);
  {
    worker_requests_t initial_send_reqs =
        dc_send_halo_to_neighbours(process, current_data);
    current_halos = dc_receive_halos(process);
    dc_concatenate_worker_requests(&initial_send_reqs, &current_halos.requests);
    MPI_Waitall(initial_send_reqs.count, initial_send_reqs.requests,
                MPI_STATUS_IGNORE);
    dc_free_worker_requests(&initial_send_reqs);
  }

  for (unsigned int i = 0; i < process.iterations; i++) {
    future_halos = dc_receive_halos(process);
    dc_compute_boundaries(&process, next_data, current_data, &current_halos);
    worker_requests_t send_reqs =
        dc_send_halo_to_neighbours(process, next_data);
    dc_compute_interior(&process, next_data, current_data);

    dc_concatenate_worker_requests(&send_reqs, &future_halos.requests);
    MPI_Waitall(send_reqs.count, send_reqs.requests, MPI_STATUS_IGNORE);

    dc_free_worker_requests(&send_reqs);
    dc_free_worker_halos(&current_halos);
    current_halos = future_halos;

    temp_ptr = current_data;
    current_data = next_data;
    next_data = temp_ptr;
  }

  dc_free_worker_halos(&current_halos);
  memcpy(process.data, current_data, process.count * sizeof(float));
  free(current_data);
  free(next_data);
  dc_log_info(process.rank, "Processing complete.");
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
    for (size_t i = 0; i < halos->halo_count; i++) {
      free(halos->halo_data[i]);
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
  free(process.indices);
  free(process.data);
}

void dc_concatenate_worker_requests(worker_requests_t *target,
                                    worker_requests_t *source) {
  if (source == NULL || source->count == 0)
    return;

  size_t original_target_count = target->count;
  size_t new_count = original_target_count + source->count;

  target->requests = realloc(target->requests, new_count * sizeof(MPI_Request));
  memcpy(target->requests + original_target_count, source->requests,
         source->count * sizeof(MPI_Request));

  if (source->buffers_to_free != NULL) {
    if (target->buffers_to_free == NULL) {
      target->buffers_to_free = malloc(new_count * sizeof(void *));
      memset(target->buffers_to_free, 0,
             original_target_count * sizeof(void *));
    } else {
      target->buffers_to_free =
          realloc(target->buffers_to_free, new_count * sizeof(void *));
    }
    memcpy(target->buffers_to_free + original_target_count,
           source->buffers_to_free, source->count * sizeof(void *));
  } else if (target->buffers_to_free != NULL) {
    target->buffers_to_free =
        realloc(target->buffers_to_free, new_count * sizeof(void *));
    memset(target->buffers_to_free + original_target_count, 0,
           source->count * sizeof(void *));
  }

  target->count = new_count;
}
