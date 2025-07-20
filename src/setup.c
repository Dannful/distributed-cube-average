#include "../include/setup.h"
#include "mpi.h"

void mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]) {
  const int periods[DIMENSIONS] = {0, 0, 0};
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, topology, periods, reorder, communicator);
}

mpi_process_t mpi_process_init(MPI_Comm communicator, int rank, unsigned int topology[DIMENSIONS]) {
  mpi_process_t process;
  process.communicator = communicator;
  MPI_Cart_coords(communicator, rank, DIMENSIONS, process.coordinates);

  for (int i = 0; i < DIMENSIONS; i++) {
    MPI_Cart_shift(communicator, i, 1, process.neighbours + 2 * i, process.neighbours + 2 * i + 1);
  }

  return process;
}
