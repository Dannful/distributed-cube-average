#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "worker.h"
#include "setup.h"
#include "coordinator.h"
#include "log.h"

void dc_extract_coordinates(size_t *position_x, size_t *position_y, size_t *position_z, size_t size_x, size_t size_y, size_t size_z, int index) {
  *position_x = index % size_x;
  *position_y = (index / size_x) % size_y;
  *position_z = index / (size_x * size_y);
}

unsigned int dc_get_index_for_coordinates(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z) {
  return position_x + position_y * size_x + position_z * size_x * size_y;
}

void dc_worker_receive_data(dc_process_t *process) {
  MPI_Recv(&process->stencil_size, 1, MPI_UINT32_T, COORDINATOR, MPI_ANY_TAG, process->communicator, MPI_STATUSES_IGNORE);
  MPI_Recv(&process->iterations, 1, MPI_UINT32_T, COORDINATOR, MPI_ANY_TAG, process->communicator, MPI_STATUSES_IGNORE);
  MPI_Recv(&process->count, 1, MPI_UNSIGNED_LONG, COORDINATOR, MPI_ANY_TAG, process->communicator, MPI_STATUSES_IGNORE);
  MPI_Recv(process->sizes, DIMENSIONS, MPI_UNSIGNED_LONG, COORDINATOR, MPI_ANY_TAG, process->communicator, MPI_STATUSES_IGNORE);
  process->indices = malloc(process->count * sizeof(size_t));
  process->data = malloc(process->count * sizeof(float));
  MPI_Recv(process->indices, process->count, MPI_UNSIGNED_LONG, COORDINATOR, MPI_ANY_TAG, process->communicator, MPI_STATUSES_IGNORE);
  MPI_Recv(process->data, process->count, MPI_FLOAT, COORDINATOR, MPI_ANY_TAG, process->communicator, MPI_STATUSES_IGNORE);
}

worker_requests_t dc_send_halo_to_neighbours(dc_process_t process, float *from) {
  size_t cell_size = (process.stencil_size - 1) / 2;
  worker_requests_t worker_requests;
  worker_requests.count = 0;
  worker_requests.requests = NULL;

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        if (dx + dy + dz != 1 && dx + dy + dz != -1) continue;

        int displacement[DIMENSIONS] = {dx, dy, dz};
        int target_coords[DIMENSIONS];
        unsigned char is_displacement_valid = 1;
        for(unsigned int i = 0; i < DIMENSIONS; i++) {
          if(process.coordinates[i] + displacement[i] < 0 || process.coordinates[i] + displacement[i] >= process.topology[i]) {
            is_displacement_valid = 0;
            break;
          }
          target_coords[i] = process.coordinates[i] + displacement[i];
        }
        if(!is_displacement_valid) continue;

        int neighbour_rank;
        MPI_Cart_rank(process.communicator, target_coords, &neighbour_rank);

        if (neighbour_rank < 0) continue;

        size_t data_size = 1;
        size_t data_dims[DIMENSIONS];
        for(unsigned int i = 0; i < DIMENSIONS; i++) {
          data_dims[i] = (displacement[i] == 0) ? process.sizes[i] : cell_size;
          data_size *= data_dims[i];
        }

        float *send_data = malloc(sizeof(float) * data_size);
        size_t *send_indices = malloc(sizeof(size_t) * data_size);

        dc_log_info(process.rank, "Preparing %zu elements for neighbor %d at displacement (%d,%d,%d)",
                  data_size, neighbour_rank, dx, dy, dz);

        size_t start_coords[DIMENSIONS], end_coords[DIMENSIONS];
        for(unsigned int i = 0; i < DIMENSIONS; i++) {
            start_coords[i] = (displacement[i] > 0) ? (process.sizes[i] - cell_size) : 0;
            end_coords[i]   = (displacement[i] < 0) ? cell_size : process.sizes[i];
        }

        size_t data_index = 0;
        for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
          for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
            for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
              size_t index = dc_get_index_for_coordinates(x, y, z, process.sizes[0], process.sizes[1], process.sizes[2]);
              send_indices[data_index] = index;
              send_data[data_index++] = from[index];
            }
          }
        }

        MPI_Request request;
        MPI_Status status;

        MPI_Isend(&data_size, 1, MPI_UNSIGNED_LONG, neighbour_rank, 0, process.communicator, &request);
        MPI_Isend(send_indices, data_size, MPI_UNSIGNED_LONG, neighbour_rank, 0, process.communicator, &request);
        MPI_Isend(send_data, data_size, MPI_FLOAT, neighbour_rank, 0, process.communicator, &request);
        free(send_indices);
        free(send_data);

        if(worker_requests.requests == NULL) {
          worker_requests.requests = malloc(sizeof(MPI_Request));
        } else {
          worker_requests.requests = realloc(worker_requests.requests, (worker_requests.count + 1) * sizeof(MPI_Request));
        }
        worker_requests.requests[worker_requests.count++] = request;
      }
    }
  }
  return worker_requests;
}

worker_halos_t dc_receive_halos(dc_process_t process) {
  worker_halos_t result;
  result.halo_count = 0;
  result.requests.count = 0;
  result.requests.requests = NULL;
  result.halo_sizes = NULL;
  result.halo_indices = NULL;
  result.halo_data = NULL;
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        if (dx + dy + dz != 1 && dx + dy + dz != -1) continue;

        int displacement[DIMENSIONS] = {dx, dy, dz};
        int target_coords[DIMENSIONS];
        unsigned char is_displacement_valid = 1;
        for(unsigned int i = 0; i < DIMENSIONS; i++) {
          if(process.coordinates[i] + displacement[i] < 0 || process.coordinates[i] + displacement[i] >= process.topology[i]) {
            is_displacement_valid = 0;
            break;
          }
          target_coords[i] = process.coordinates[i] + displacement[i];
        }
        if(!is_displacement_valid) continue;

        int neighbour_rank;

        if(result.requests.requests == NULL) {
          result.requests.requests = malloc(sizeof(MPI_Request));
        } else {
          result.requests.requests = realloc(result.requests.requests, (result.requests.count + 1) * sizeof(MPI_Request));
        }

        size_t recv_data_size;
        MPI_Irecv(&recv_data_size, 1, MPI_UNSIGNED_LONG, neighbour_rank, 0, process.communicator, result.requests.requests + result.requests.count);

        if(result.halo_data == NULL) {
          result.halo_data = malloc(sizeof(float *));
          result.halo_indices = malloc(sizeof(size_t *));
          result.halo_sizes = malloc(sizeof(size_t));
        } else {
          result.halo_data = realloc(result.halo_data, (result.halo_count + 1) * sizeof(float *));
          result.halo_indices = realloc(result.halo_indices, (result.halo_count + 1) * sizeof(size_t *));
          result.halo_sizes = realloc(result.halo_sizes, (result.halo_count + 1) * sizeof(size_t));
        }

        result.halo_data[result.halo_count] = malloc(recv_data_size * sizeof(float));
        result.halo_indices[result.halo_count] = malloc(recv_data_size * sizeof(size_t));
        result.halo_sizes[result.halo_count] = recv_data_size;

        MPI_Irecv(result.halo_indices[result.halo_count], recv_data_size, MPI_UNSIGNED_LONG, neighbour_rank, 0, process.communicator, result.requests.requests + result.requests.count);
        MPI_Irecv(result.halo_data[result.halo_count], recv_data_size, MPI_FLOAT, neighbour_rank, 0, process.communicator, result.requests.requests + result.requests.count);

        result.requests.count++;
        result.halo_count++;
      }
    }
  }
  return result;
}

void dc_worker_process(dc_process_t process) {
  worker_requests_t send_requests = dc_send_halo_to_neighbours(process, process.data);
  worker_halos_t input_halos = dc_receive_halos(process);
  worker_halos_t output_halos;
  dc_concatenate_worker_requests(&send_requests, &input_halos.requests);
  MPI_Waitall(send_requests.count, send_requests.requests, MPI_STATUSES_IGNORE);
  dc_free_worker_requests(&send_requests);
  float *output_data = malloc(process.count * sizeof(float));
  memccpy(output_data, process.data, sizeof(float), process.count);
  for(unsigned int i = 0; i < process.iterations; i++) {
    output_halos = dc_receive_halos(process);
    // TODO: Process the halos
    dc_send_halo_to_neighbours(process, output_data);
    // TODO: Process my data and write to output_data
    MPI_Waitall(output_halos.requests.count, output_halos.requests.requests, MPI_STATUSES_IGNORE);
  }
  free(output_data);
  dc_free_worker_halos(&input_halos);
}

void dc_free_worker_requests(worker_requests_t *requests) {
  if (requests->requests != NULL) {
    free(requests->requests);
    requests->requests = NULL;
  }
  requests->count = 0;
}

void dc_free_worker_halos(worker_halos_t *halos) {
  for (size_t i = 0; i < halos->halo_count; i++) {
    free(halos->halo_indices[i]);
    free(halos->halo_data[i]);
  }
  free(halos->halo_indices);
  free(halos->halo_data);
  free(halos->halo_sizes);
  dc_free_worker_requests(&halos->requests);
}

void dc_worker_free(dc_process_t process) {
  free(process.indices);
  free(process.data);
}

void dc_concatenate_worker_requests(worker_requests_t *target, worker_requests_t *source) {
  target->requests = realloc(target->requests, (target->count + source->count) * sizeof(MPI_Request));
  memcpy(target->requests + target->count, source->requests, source->count);
}
