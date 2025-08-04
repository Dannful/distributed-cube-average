#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/coordinator.h"
#include "../include/setup.h"
#include "../include/log.h"
#include "../include/worker.h"

int main(int argc, char **argv) {
  MPI_Comm communicator;
  int topology[DIMENSIONS] = {0};
  int rank, size;
  char hostname[256];
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, DIMENSIONS, topology);
  mpi_world_init(&communicator, topology);
  MPI_Comm_rank(communicator, &rank);

  mpi_process_t mpi_process = mpi_process_init(communicator, rank, topology);

  const size_t size_x = 100;
  const size_t size_y = 100;
  const size_t size_z = 100;
  const unsigned int iterations = 6, stencil_size = 3;

  if(rank == COORDINATOR) {
    log_info(rank, "Initializing problem data...");
    problem_data_t problem_data = init_problem_data(communicator, (unsigned int*) topology, size, iterations, stencil_size, size_x, size_y, size_z);
    log_info(rank, "Partitioning cube...");
    partition_cube(problem_data);
    log_info(rank, "Partition completed. Sending data to workers...");
    send_data_to_workers(problem_data);
    mpi_process.data = malloc(sizeof(float) * problem_data.worker_count[0]);
    mpi_process.count = problem_data.worker_count[0];
    mpi_process.indices = malloc(sizeof(size_t) * problem_data.worker_count[0]);
    memmove(mpi_process.data, problem_data.workers[0], sizeof(float) * problem_data.worker_count[0]);
    memmove(mpi_process.indices, problem_data.worker_indices[0], sizeof(size_t) * problem_data.worker_count[0]);
    memmove(mpi_process.sizes, problem_data.worker_sizes[0], sizeof(size_t) * DIMENSIONS);
    mpi_process.stencil_size = problem_data.stencil_size;
    mpi_process.iterations = problem_data.iterations;
    mpi_process.rank = rank;
    free_problem_data(problem_data);
  } else {
    log_info(rank, "Awaiting data from coordinator...");
    receive_worker_data(&mpi_process);
    log_info(rank, "Data received from coordinator.");
  }
  log_info(rank, "Starting worker process...");
  worker_process(mpi_process);
  worker_free(mpi_process);
  MPI_Finalize();
  return 0;
}
