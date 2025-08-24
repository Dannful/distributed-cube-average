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
  MPI_Recv(process->sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR,
           MPI_ANY_TAG, process->communicator, MPI_STATUS_IGNORE);

  process->data =
      malloc(dc_compute_count_from_sizes(process->sizes) * sizeof(float));

  MPI_Recv(process->data, dc_compute_count_from_sizes(process->sizes),
           MPI_FLOAT, COORDINATOR, MPI_ANY_TAG, process->communicator,
           MPI_STATUS_IGNORE);
}

worker_requests_t dc_send_halo_to_neighbours(dc_process_t process,
                                             float *from) {
  worker_requests_t reqs;
  size_t radius = process.stencil_size;

  reqs.count = 0;
  reqs.requests = malloc(6 * sizeof(MPI_Request));
  reqs.buffers_to_free = malloc(6 * sizeof(void *));

  for (int dim = 0; dim < DIMENSIONS; dim++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int neighbour_rank;
      int target_coords[DIMENSIONS];
      memcpy(target_coords, process.coordinates, sizeof(int) * DIMENSIONS);
      target_coords[dim] += dir;

      if (target_coords[dim] < 0 ||
          target_coords[dim] >= process.topology[dim]) {
        continue;
      }
      MPI_Cart_rank(process.communicator, target_coords, &neighbour_rank);

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
      for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
        for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
          for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
            send_buffer[data_index++] = from[dc_get_index_for_coordinates(
                x, y, z, process.sizes[0], process.sizes[1], process.sizes[2])];
          }
        }
      }

      dc_log_info(process.rank, "Sending halo of %zu elements to neighbor %d",
                  data_size, neighbour_rank);

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
  size_t radius = process.stencil_size;
  result.halo_count = 0;

  result.requests.count = 0;
  result.requests.requests = malloc(2 * DIMENSIONS * sizeof(MPI_Request));
  result.requests.buffers_to_free = NULL;

  result.halo_sizes = calloc(2 * DIMENSIONS, sizeof(size_t));
  result.halo_data = calloc(2 * DIMENSIONS, sizeof(float *));

  for (int dim = 0; dim < DIMENSIONS; dim++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int source_rank;
      int source_coords[DIMENSIONS];
      memcpy(source_coords, process.coordinates, sizeof(int) * DIMENSIONS);
      source_coords[dim] += dir;

      if (source_coords[dim] < 0 ||
          source_coords[dim] >= process.topology[dim]) {
        continue;
      }
      MPI_Cart_rank(process.communicator, source_coords, &source_rank);

      size_t recv_data_size = 1;
      for (unsigned int i = 0; i < DIMENSIONS; i++) {
        recv_data_size *= (i == dim) ? radius : process.sizes[i];
      }

      size_t face_index = 2 * dim + (dir + 1) / 2;
      result.halo_sizes[face_index] = recv_data_size;
      result.halo_data[face_index] = malloc(recv_data_size * sizeof(float));

      MPI_Irecv(result.halo_data[face_index], recv_data_size, MPI_FLOAT,
                source_rank, 0, process.communicator,
                &result.requests.requests[result.requests.count]);

      result.halo_count++;
      result.requests.count++;
    }
  }
  return result;
}

static int get_value_at_coord(int x, int y, int z, const dc_process_t *process,
                              const float *input_data,
                              const worker_halos_t *halos, float *out_value) {
  const size_t size_x = process->sizes[0];
  const size_t size_y = process->sizes[1];
  const size_t size_z = process->sizes[2];
  const int radius = process->stencil_size;

  if (x >= 0 && x < size_x && y >= 0 && y < size_y && z >= 0 && z < size_z) {
    *out_value = input_data[dc_get_index_for_coordinates(x, y, z, size_x,
                                                         size_y, size_z)];
    return 1;
  }

  int halo_idx = -1;
  size_t halo_x, halo_y, halo_z;
  size_t halo_size_x, halo_size_y, halo_size_z;

  if (x < 0) {
    halo_idx = 0;
    halo_x = x + radius;
    halo_y = y;
    halo_z = z;
    halo_size_x = radius;
    halo_size_y = size_y;
    halo_size_z = size_z;
  } else if (x >= size_x) {
    halo_idx = 1;
    halo_x = x - size_x;
    halo_y = y;
    halo_z = z;
    halo_size_x = radius;
    halo_size_y = size_y;
    halo_size_z = size_z;
  } else if (y < 0) {
    halo_idx = 2;
    halo_x = x;
    halo_y = y + radius;
    halo_z = z;
    halo_size_x = size_x;
    halo_size_y = radius;
    halo_size_z = size_z;
  } else if (y >= size_y) {
    halo_idx = 3;
    halo_x = x;
    halo_y = y - size_y;
    halo_z = z;
    halo_size_x = size_x;
    halo_size_y = radius;
    halo_size_z = size_z;
  } else if (z < 0) {
    halo_idx = 4;
    halo_x = x;
    halo_y = y;
    halo_z = z + radius;
    halo_size_x = size_x;
    halo_size_y = size_y;
    halo_size_z = radius;
  } else if (z >= size_z) {
    halo_idx = 5;
    halo_x = x;
    halo_y = y;
    halo_z = z - size_z;
    halo_size_x = size_x;
    halo_size_y = size_y;
    halo_size_z = radius;
  }

  if (halo_idx == -1 || halos->halo_sizes[halo_idx] == 0)
    return 0;

  size_t halo_internal_index = dc_get_index_for_coordinates(
      halo_x, halo_y, halo_z, halo_size_x, halo_size_y, halo_size_z);
  if (halo_internal_index < halos->halo_sizes[halo_idx]) {
    *out_value = halos->halo_data[halo_idx][halo_internal_index];
    return 1;
  }
  return 0;
}

static void dc_compute_stencil_for_point(size_t x, size_t y, size_t z,
                                         const dc_process_t *process,
                                         float *output_data,
                                         const float *input_data,
                                         const worker_halos_t *halos) {
  const int radius = process->stencil_size;
  double sum = 0.0;
  int count = 0;
  float neighbour_value;

  for (int displacement = 0; displacement <= radius; displacement++) {
    if (displacement == 0) {
      if (get_value_at_coord(x, y, z, process, input_data, halos,
                             &neighbour_value)) {
        sum += neighbour_value;
        count++;
      }
      continue;
    }
    if (get_value_at_coord(x + displacement, y, z, process, input_data, halos,
                           &neighbour_value)) {
      sum += neighbour_value;
      count++;
    }
    if (get_value_at_coord(x - displacement, y, z, process, input_data, halos,
                           &neighbour_value)) {
      sum += neighbour_value;
      count++;
    }
    if (get_value_at_coord(x, y + displacement, z, process, input_data, halos,
                           &neighbour_value)) {
      sum += neighbour_value;
      count++;
    }
    if (get_value_at_coord(x, y - displacement, z, process, input_data, halos,
                           &neighbour_value)) {
      sum += neighbour_value;
      count++;
    }
    if (get_value_at_coord(x, y, z + displacement, process, input_data, halos,
                           &neighbour_value)) {
      sum += neighbour_value;
      count++;
    }
    if (get_value_at_coord(x, y, z - displacement, process, input_data, halos,
                           &neighbour_value)) {
      sum += neighbour_value;
      count++;
    }
  }
  size_t local_idx = dc_get_index_for_coordinates(
      x, y, z, process->sizes[0], process->sizes[1], process->sizes[2]);
  if (count > 0) {
    output_data[local_idx] = sum / count;
  } else {
    float original_value;
    if (get_value_at_coord(x, y, z, process, input_data, halos,
                           &original_value)) {
      output_data[local_idx] = original_value;
    }
  }
}

void dc_compute_boundaries(const dc_process_t *process, float *output_data,
                           const float *input_data,
                           const worker_halos_t *halos) {
  const int radius = process->stencil_size;
  const size_t size_x = process->sizes[0];
  const size_t size_y = process->sizes[1];
  const size_t size_z = process->sizes[2];

  for (size_t x_layer = 0; x_layer < radius; x_layer++) {
    size_t x_left = x_layer;
    size_t x_right = size_x - 1 - x_layer;
    for (size_t y = 0; y < size_y; y++) {
      for (size_t z = 0; z < size_z; z++) {
        dc_compute_stencil_for_point(x_left, y, z, process, output_data,
                                     input_data, halos);
        if (x_left != x_right)
          dc_compute_stencil_for_point(x_right, y, z, process, output_data,
                                       input_data, halos);
      }
    }
  }

  for (size_t y_layer = 0; y_layer < radius; y_layer++) {
    size_t y_bottom = y_layer;
    size_t y_top = size_y - 1 - y_layer;
    for (size_t x = radius; x < size_x - radius; x++) {
      for (size_t z = 0; z < size_z; z++) {
        dc_compute_stencil_for_point(x, y_bottom, z, process, output_data,
                                     input_data, halos);
        if (y_bottom != y_top)
          dc_compute_stencil_for_point(x, y_top, z, process, output_data,
                                       input_data, halos);
      }
    }
  }

  for (size_t z_layer = 0; z_layer < radius; z_layer++) {
    size_t z_back = z_layer;
    size_t z_front = size_z - 1 - z_layer;
    for (size_t x = radius; x < size_x - radius; x++) {
      for (size_t y = radius; y < size_y - radius; y++) {
        dc_compute_stencil_for_point(x, y, z_back, process, output_data,
                                     input_data, halos);
        if (z_back != z_front)
          dc_compute_stencil_for_point(x, y, z_front, process, output_data,
                                       input_data, halos);
      }
    }
  }
}

void dc_compute_interior(const dc_process_t *process, float *output_data,
                         const float *input_data) {
  const int radius = process->stencil_size;
  const size_t size_x = process->sizes[0];
  const size_t size_y = process->sizes[1];
  const size_t size_z = process->sizes[2];

  for (size_t z = radius; z < size_z - radius; z++) {
    for (size_t y = radius; y < size_y - radius; y++) {
      for (size_t x = radius; x < size_x - radius; x++) {

        double sum = 0.0;
        int count = 0;

        for (int displacement = 0; displacement <= radius; displacement++) {
          if (displacement == 0) {
            sum += input_data[dc_get_index_for_coordinates(x, y, z, size_x,
                                                           size_y, size_z)];
            count++;
            continue;
          }
          sum += input_data[dc_get_index_for_coordinates(
              x + displacement, y, z, size_x, size_y, size_z)];
          sum += input_data[dc_get_index_for_coordinates(
              x - displacement, y, z, size_x, size_y, size_z)];
          sum += input_data[dc_get_index_for_coordinates(
              x, y + displacement, z, size_x, size_y, size_z)];
          sum += input_data[dc_get_index_for_coordinates(
              x, y - displacement, z, size_x, size_y, size_z)];
          sum += input_data[dc_get_index_for_coordinates(
              x, y, z + displacement, size_x, size_y, size_z)];
          sum += input_data[dc_get_index_for_coordinates(
              x, y, z - displacement, size_x, size_y, size_z)];
          count += 6;
        }

        size_t current_index =
            dc_get_index_for_coordinates(x, y, z, size_x, size_y, size_z);
        output_data[current_index] = sum / count;
      }
    }
  }
}

void dc_send_data_to_coordinator(dc_process_t process) {
  if (process.rank == COORDINATOR)
    return;
  MPI_Send(process.sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR, 0,
           process.communicator);
  MPI_Send(process.data, dc_compute_count_from_sizes(process.sizes), MPI_FLOAT,
           COORDINATOR, 0, process.communicator);
}

void dc_worker_process(dc_process_t process) {
  size_t count = dc_compute_count_from_sizes(process.sizes);
  float *current_data = malloc(count * sizeof(float));
  memcpy(current_data, process.data, count * sizeof(float));
  float *next_data = malloc(count * sizeof(float));
  float *temp_ptr;

  worker_requests_t all_send_requests;
  all_send_requests.buffers_to_free = NULL;
  all_send_requests.requests = NULL;
  all_send_requests.count = 0;
  worker_halos_t current_halos, future_halos;
  dc_log_info(process.rank, "Starting %u iterations", process.iterations);

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

    dc_concatenate_worker_requests(&all_send_requests, &send_reqs);

    MPI_Waitall(future_halos.requests.count, future_halos.requests.requests,
                MPI_STATUS_IGNORE);

    dc_free_worker_halos(&current_halos);
    current_halos = future_halos;
    temp_ptr = current_data;
    current_data = next_data;
    next_data = temp_ptr;
  }

  MPI_Waitall(all_send_requests.count, all_send_requests.requests,
              MPI_STATUS_IGNORE);
  dc_free_worker_requests(&all_send_requests);

  dc_free_worker_halos(&current_halos);
  memcpy(process.data, current_data, count * sizeof(float));
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
    for (size_t i = 0; i < DIMENSIONS * 2; i++) {
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

void dc_worker_free(dc_process_t process) { free(process.data); }

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

  free(source->requests);
  free(source->buffers_to_free);
  source->requests = NULL;
  source->buffers_to_free = NULL;
  source->count = 0;
}

size_t dc_compute_count_from_sizes(size_t sizes[DIMENSIONS]) {
  size_t count = 1;
  for (int i = 0; i < DIMENSIONS; i++) {
    count *= sizes[i];
  }
  return count;
}
