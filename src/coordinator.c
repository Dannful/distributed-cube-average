#include "../include/coordinator.h"
#include "../include/utils.h"
#include "../include/log.h"
#include "../include/setup.h"
#include "../include/worker.h"
#include "mpi.h"
#include <string.h>

problem_data_t init_problem_data(MPI_Comm comm, unsigned int topology[DIMENSIONS], unsigned int workers, unsigned int iterations, unsigned int stencil_size, size_t size_x, size_t size_y, size_t size_z) {
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
  if(result.cube == NULL) {
    log_error(0, "Failed to allocate memory for cube data.");
    exit(EXIT_FAILURE);
  }
  for(size_t i = 0; i < size_x * size_y * size_z; i++)
    result.cube[i] = (float)i;
  result.worker_count = calloc(workers, sizeof(size_t));
  if(result.worker_count == NULL) {
    log_error(0, "Failed to allocate memory for worker counts.");
    free(result.cube);
    exit(EXIT_FAILURE);
  }
  result.worker_indices = calloc(workers, sizeof(size_t *));
  if(result.worker_indices == NULL) {
    log_error(0, "Failed to allocate memory for worker indices.");
    free(result.cube);
    free(result.worker_count);
    exit(EXIT_FAILURE);
  }
  result.workers = calloc(workers, sizeof(float *));
  if(result.workers == NULL) {
    log_error(0, "Failed to allocate memory for worker data.");
    free(result.cube);
    free(result.worker_count);
    free(result.worker_indices);
    exit(EXIT_FAILURE);
  }
  result.worker_sizes = calloc(workers, sizeof(size_t*));
  if(result.worker_sizes == NULL) {
    log_error(0, "Failed to allocate memory for worker sizes.");
    free(result.cube);
    free(result.worker_count);
    free(result.worker_indices);
    free(result.workers);
    exit(EXIT_FAILURE);
  }
  return result;
}

void partition_cube(problem_data_t problem_data) {
  float ***unflattened_cube = unflatten_cube(problem_data.cube, problem_data.size_x, problem_data.size_y, problem_data.size_z);
  size_t partition_size_x = problem_data.size_x / problem_data.topology[0];
  size_t partition_size_y = problem_data.size_y / problem_data.topology[1];
  size_t partition_size_z = problem_data.size_z / problem_data.topology[2];
  size_t cell_size = (problem_data.stencil_size - 1) / 2;
  size_t remainder_x = problem_data.size_x % problem_data.topology[0];
  size_t remainder_y = problem_data.size_y % problem_data.topology[1];
  size_t remainder_z = problem_data.size_z % problem_data.topology[2];
  for(size_t worker = 0; worker < problem_data.num_workers; worker++) {
    int process_coordinates[DIMENSIONS];
    MPI_Cart_coords(problem_data.communicator, worker, DIMENSIONS, process_coordinates);
    size_t worker_size_x = process_coordinates[0] == problem_data.topology[0] - 1 ? partition_size_x + remainder_x : partition_size_x;
    size_t worker_size_y = process_coordinates[1] == problem_data.topology[1] - 1 ? partition_size_y + remainder_y : partition_size_y;
    size_t worker_size_z = process_coordinates[2] == problem_data.topology[2] - 1 ? partition_size_z + remainder_z : partition_size_z;
    problem_data.worker_sizes[worker] = malloc(sizeof(size_t) * 3);
    problem_data.worker_sizes[worker][0] = worker_size_x;
    problem_data.worker_sizes[worker][1] = worker_size_y;
    problem_data.worker_sizes[worker][2] = worker_size_z;
    problem_data.workers[worker] = malloc(sizeof(float) * worker_size_x * worker_size_y * worker_size_z);
    problem_data.worker_indices[worker] = malloc(sizeof(size_t) * worker_size_x * worker_size_y * worker_size_z);
    for(size_t x =  process_coordinates[0] * partition_size_x; x < process_coordinates[0] * partition_size_x + worker_size_x; x++) {
      for(size_t y = process_coordinates[1] * partition_size_y; y < process_coordinates[1] * partition_size_y + worker_size_y; y++) {
        for(size_t z = process_coordinates[2] * partition_size_z; z < process_coordinates[2] * partition_size_z + worker_size_z; z++) {
          problem_data.workers[worker][problem_data.worker_count[worker]] = unflattened_cube[x][y][z];
          problem_data.worker_indices[worker][problem_data.worker_count[worker]++] = get_index_for_coordinates(x, y, z, problem_data.size_x, problem_data.size_y, problem_data.size_z);
        }
      }
    }
  }
  free_unflattened_cube(unflattened_cube, problem_data.size_x, problem_data.size_y, problem_data.size_z);
}

void send_data_to_workers(problem_data_t problem_data) {
  for(size_t worker = 0; worker < problem_data.num_workers; worker++) {
    if(worker == COORDINATOR) continue;
    MPI_Send(&problem_data.stencil_size, 1, MPI_UINT32_T, worker, 0, problem_data.communicator);
    MPI_Send(&problem_data.iterations, 1, MPI_UINT32_T, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.worker_count + worker, 1, MPI_UNSIGNED_LONG, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.worker_sizes[worker], DIMENSIONS, MPI_UNSIGNED_LONG, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.worker_indices[worker], problem_data.worker_count[worker], MPI_UNSIGNED_LONG, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.workers[worker], problem_data.worker_count[worker], MPI_FLOAT, worker, 0, problem_data.communicator);
  }
}

void free_problem_data(problem_data_t problem_data) {
  free(problem_data.cube);
  for(size_t i = 0; i < problem_data.num_workers; i++) {
    free(problem_data.workers[i]);
    free(problem_data.worker_indices[i]);
    free(problem_data.worker_sizes[i]);
  }
  free(problem_data.workers);
  free(problem_data.worker_indices);
  free(problem_data.worker_count);
  free(problem_data.worker_sizes);
}
