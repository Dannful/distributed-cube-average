#include <mpi.h>
#include <stdlib.h>

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

void dc_worker_process(dc_process_t process) {
  size_t cell_size = (process.stencil_size - 1) / 2;
  if (cell_size == 0) return;

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;

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
              send_data[data_index++] = process.data[index];
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

        size_t recv_data_size;
        MPI_Recv(&recv_data_size, 1, MPI_UNSIGNED_LONG, neighbour_rank, 0, process.communicator, &status);

        float *recv_data = malloc(recv_data_size * sizeof(float));
        size_t *recv_indices = malloc(recv_data_size * sizeof(size_t));

        MPI_Recv(recv_indices, recv_data_size, MPI_UNSIGNED_LONG, neighbour_rank, 0, process.communicator, &status);
        MPI_Recv(recv_data, recv_data_size, MPI_FLOAT, neighbour_rank, 0, process.communicator, &status);

        free(recv_indices);
        free(recv_data);
      }
    }
  }
}

void dc_worker_free(dc_process_t process) {
  free(process.indices);
  free(process.data);
}
