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
    free_problem_data(problem_data);
  } else {
    mpi_process_t process = receive_worker_data(communicator, rank, topology);
    worker_free(process);
  }
  MPI_Finalize();
  return 0;
}
