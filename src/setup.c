#include <mpi.h>
#include <string.h>

#include "setup.h"

void dc_mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]) {
  const int periods[DIMENSIONS] = {0, 0, 0};
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, topology, periods, reorder, communicator);
}

dc_process_t dc_process_init(MPI_Comm communicator, int rank, int topology[DIMENSIONS]) {
  dc_process_t process;
  process.communicator = communicator;
  process.rank = rank;
  memccpy(process.topology, topology, sizeof(int), DIMENSIONS);
  MPI_Cart_coords(communicator, rank, DIMENSIONS, process.coordinates);

  MPI_Cart_shift(communicator, 0, 1, process.neighbours + LEFT, process.neighbours + RIGHT);
  MPI_Cart_shift(communicator, 1, 1, process.neighbours + UP, process.neighbours + DOWN);
  MPI_Cart_shift(communicator, 2, 1, process.neighbours + FRONT, process.neighbours + BACK);

  return process;
}
