#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/setup.h"

int main(int argc, char **argv) {
  MPI_Comm communicator;
  const int topology[DIMENSIONS] = {3, 3, 3};
  int rank, size;
  char hostname[256];
  MPI_Init(&argc, &argv);
  mpi_world_init(&communicator, topology);
  MPI_Comm_size(communicator, &size);

  if(size != topology[0] * topology[1] * topology[2]) {
    fprintf(stderr, "Error: The number of MPI processes (%d) does not match the expected topology (%d x %d x %d = %d).\n", size, topology[0], topology[1], topology[2], topology[0] * topology[1] * topology[2]);
    MPI_Abort(communicator, EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  mpi_process_t mpi_process = mpi_process_init(communicator, topology);

  const size_t size_x = 100;
  const size_t size_y = 100;
  const size_t size_z = 100;
  const unsigned int iterations = 6, stencil_size = 3;

  if(rank == 0) {
    problem_data_t problem_data = init_problem_data(iterations, stencil_size, size_x, size_y, size_z);
    free(problem_data.cube);
  }
  MPI_Finalize();
  return 0;
}
