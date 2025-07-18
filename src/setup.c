#include "../include/setup.h"
#include "mpi.h"
#include <stdlib.h>
#include "../include/log.h"

void mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]) {
  const int period[DIMENSIONS] = {0, 0, 0};
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, topology, period, reorder, communicator);
}

mpi_process_t mpi_process_init(MPI_Comm communicator, const int topology[DIMENSIONS]) {
  mpi_process_t process;
  process.communicator = communicator;
  MPI_Comm_rank(communicator, &process.rank);
  MPI_Cart_coords(communicator, process.rank, DIMENSIONS, process.coordinates);

  for (int i = 0; i < DIMENSIONS; i++) {
    MPI_Cart_shift(communicator, i, 1, process.neighbours[i], process.neighbours[i] + 1);
  }

  return process;
}

problem_data_t init_problem_data(unsigned int iterations, unsigned int stencil_size, size_t size_x, size_t size_y, size_t size_z) {
  problem_data_t result;
  result.iterations = iterations;
  result.stencil_size = stencil_size;
  result.size_x = size_x;
  result.size_y = size_y;
  result.size_z = size_z;
  result.cube = (float *)malloc(size_x * size_y * size_z * sizeof(float));
  if(result.cube == NULL) {
    log_error(0, "Failed to allocate memory for cube data.");
    exit(EXIT_FAILURE);
  }
  for(size_t i = 0; i < size_x * size_y * size_z; i++)
    result.cube[i] = (float)i;
  return result;
}
