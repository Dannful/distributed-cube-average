#include <mpi.h>
#include <string.h>

#include "coordinator.h"
#include "log.h"
#include "setup.h"
#include "worker.h"

problem_data_t dc_initialize_problem(MPI_Comm comm,
                                     unsigned int topology[DIMENSIONS],
                                     unsigned int workers,
                                     unsigned int iterations,
                                     unsigned int stencil_size, size_t size_x,
                                     size_t size_y, size_t size_z) {
  problem_data_t result;
  result.communicator = comm;
  result.num_workers = workers;
  result.iterations = iterations;
  result.stencil_size = stencil_size;
  result.size_x = size_x;
  result.size_y = size_y;
  result.size_z = size_z;
  result.cube = (float *)malloc(size_x * size_y * size_z * sizeof(float));
  memcpy(result.topology, topology, DIMENSIONS * sizeof(unsigned int));
  if (result.cube == NULL) {
    dc_log_error(0, "Failed to allocate memory for cube data.");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < size_x * size_y * size_z; i++) {
    result.cube[i] = (float)(i + 1);
  }
  result.workers = calloc(workers, sizeof(float *));
  if (result.workers == NULL) {
    dc_log_error(0, "Failed to allocate memory for worker data.");
    free(result.cube);
    exit(EXIT_FAILURE);
  }
  result.worker_sizes = calloc(workers, sizeof(size_t *));
  if (result.worker_sizes == NULL) {
    dc_log_error(0, "Failed to allocate memory for worker sizes.");
    free(result.cube);
    free(result.workers);
    exit(EXIT_FAILURE);
  }
  return result;
}

void dc_partition_cube(problem_data_t problem_data) {
  size_t partition_size_x = problem_data.size_x / problem_data.topology[0];
  size_t partition_size_y = problem_data.size_y / problem_data.topology[1];
  size_t partition_size_z = problem_data.size_z / problem_data.topology[2];
  size_t remainder_x = problem_data.size_x % problem_data.topology[0];
  size_t remainder_y = problem_data.size_y % problem_data.topology[1];
  size_t remainder_z = problem_data.size_z % problem_data.topology[2];
  for (size_t worker = 0; worker < problem_data.num_workers; worker++) {
    int process_coordinates[DIMENSIONS];
    MPI_Cart_coords(problem_data.communicator, worker, DIMENSIONS,
                    process_coordinates);
    size_t worker_size_x =
        process_coordinates[0] == problem_data.topology[0] - 1
            ? partition_size_x + remainder_x
            : partition_size_x;
    size_t worker_size_y =
        process_coordinates[1] == problem_data.topology[1] - 1
            ? partition_size_y + remainder_y
            : partition_size_y;
    size_t worker_size_z =
        process_coordinates[2] == problem_data.topology[2] - 1
            ? partition_size_z + remainder_z
            : partition_size_z;
    problem_data.worker_sizes[worker] = malloc(sizeof(size_t) * DIMENSIONS);
    problem_data.worker_sizes[worker][0] = worker_size_x;
    problem_data.worker_sizes[worker][1] = worker_size_y;
    problem_data.worker_sizes[worker][2] = worker_size_z;
    problem_data.workers[worker] =
        malloc(sizeof(float) * worker_size_x * worker_size_y * worker_size_z);
    size_t count = 0;
    dc_log_info(0, "Sizes: %zu %zu %zu %zu %zu %zu", process_coordinates[0],
                process_coordinates[1], process_coordinates[2], worker_size_x,
                worker_size_y, worker_size_z);
    for (size_t z = process_coordinates[2] * partition_size_z;
         z < process_coordinates[2] * partition_size_z + worker_size_z; z++) {
      for (size_t y = process_coordinates[1] * partition_size_y;
           y < process_coordinates[1] * partition_size_y + worker_size_y; y++) {
        for (size_t x = process_coordinates[0] * partition_size_x;
             x < process_coordinates[0] * partition_size_x + worker_size_x;
             x++) {
          if (process_coordinates[0] == 1 && process_coordinates[1] == 1 &&
              process_coordinates[2] == 0) {
            dc_log_info(0, "brocio %zu %zu %zu", x, y, z);
          }
          size_t index = dc_get_index_for_coordinates(
              x, y, z, problem_data.size_x, problem_data.size_y,
              problem_data.size_z);
          problem_data.workers[worker][count++] = problem_data.cube[index];
        }
      }
    }
  }
}

void dc_send_data_to_workers(problem_data_t problem_data) {
  for (size_t worker = 0; worker < problem_data.num_workers; worker++) {
    if (worker == COORDINATOR)
      continue;
    MPI_Send(&problem_data.stencil_size, 1, MPI_UINT32_T, worker, 0,
             problem_data.communicator);
    MPI_Send(&problem_data.iterations, 1, MPI_UINT32_T, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.worker_sizes[worker], DIMENSIONS, MPI_UNSIGNED_LONG,
             worker, 0, problem_data.communicator);
    MPI_Send(problem_data.workers[worker],
             dc_compute_count_from_sizes(problem_data.worker_sizes[worker]),
             MPI_FLOAT, worker, 0, problem_data.communicator);
  }
}

float *dc_receive_data_from_workers(dc_process_t coordinator_process,
                                    size_t cube_size_x, size_t cube_size_y,
                                    size_t cube_size_z) {
  float *cube =
      (float *)malloc(sizeof(float) * cube_size_x * cube_size_y * cube_size_z);
  size_t size_x = cube_size_x / coordinator_process.topology[0];
  size_t size_y = cube_size_y / coordinator_process.topology[1];
  size_t size_z = cube_size_z / coordinator_process.topology[2];
  for (int worker_x = 0; worker_x < coordinator_process.topology[0];
       worker_x++) {
    for (int worker_y = 0; worker_y < coordinator_process.topology[1];
         worker_y++) {
      for (int worker_z = 0; worker_z < coordinator_process.topology[2];
           worker_z++) {
        int worker_coordinates[DIMENSIONS] = {worker_x, worker_y, worker_z};
        int worker_rank;
        MPI_Cart_rank(coordinator_process.communicator, worker_coordinates,
                      &worker_rank);
        size_t *worker_sizes = NULL;
        float *worker_data = NULL;
        size_t worker_count = 0;
        if (worker_rank == COORDINATOR) {
          worker_sizes = coordinator_process.sizes;
          worker_data = coordinator_process.data;
          worker_count = dc_compute_count_from_sizes(worker_sizes);
        } else {
          worker_sizes = malloc(DIMENSIONS * sizeof(size_t));
          MPI_Recv(worker_sizes, DIMENSIONS, MPI_UNSIGNED_LONG, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
          worker_count = dc_compute_count_from_sizes(worker_sizes);
          worker_data = malloc(sizeof(float) * worker_count);
          MPI_Recv(worker_data, worker_count, MPI_FLOAT, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
        }
        for (size_t x = 0; x < worker_sizes[0]; x++) {
          for (size_t y = 0; y < worker_sizes[1]; y++) {
            for (size_t z = 0; z < worker_sizes[2]; z++) {
              size_t worker_index = dc_get_index_for_coordinates(
                  x, y, z, worker_sizes[0], worker_sizes[1], worker_sizes[2]);
              size_t global_index = dc_get_index_for_coordinates(
                  worker_x * size_x + x, worker_y * size_y + y,
                  worker_z * size_z + z, cube_size_x, cube_size_y, cube_size_z);
              cube[global_index] = worker_data[worker_index];
            }
          }
        }
        if (worker_rank != COORDINATOR) {
          free(worker_sizes);
          free(worker_data);
        }
      }
    }
  }
  return cube;
}

void dc_free_problem_data_mem(problem_data_t problem_data) {
  free(problem_data.cube);
  for (size_t i = 0; i < problem_data.num_workers; i++) {
    free(problem_data.workers[i]);
    free(problem_data.worker_sizes[i]);
  }
  free(problem_data.workers);
  free(problem_data.worker_sizes);
}
