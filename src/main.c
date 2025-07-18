#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/coordinator.h"
#include "../include/setup.h"
#include "../include/log.h"

int main(int argc, char **argv) {
  MPI_Comm communicator;
  const int topology[DIMENSIONS] = {3, 3, 3};
  int rank, size;
  char hostname[256];
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  mpi_world_init(&communicator, topology);
  MPI_Comm_rank(communicator, &rank);

  if(size != topology[0] * topology[1] * topology[2]) {
    fprintf(stderr, "Error: The number of MPI processes (%d) does not match the expected topology (%d x %d x %d = %d).\n", size, topology[0], topology[1], topology[2], topology[0] * topology[1] * topology[2]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  mpi_process_t mpi_process = mpi_process_init(communicator, rank, (unsigned int*) topology);

  const size_t size_x = 100;
  const size_t size_y = 100;
  const size_t size_z = 100;
  const unsigned int iterations = 6, stencil_size = 3;

  if(rank == 0) {
    log_info(rank, "Initializing problem data...");
    problem_data_t problem_data = init_problem_data(communicator, (unsigned int*) topology, size, iterations, stencil_size, size_x, size_y, size_z);
    log_info(rank, "Partitioning cube...");
    partition_cube(problem_data);
    free_problem_data(problem_data);
  }
  MPI_Finalize();
  return 0;
}
