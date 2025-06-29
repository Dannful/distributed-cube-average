#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/coordinator.h"
#include "../include/log.h"

int main(int argc, char **argv) {
  int rank;
  char hostname[256];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const size_t sizeX = 100;
  const size_t sizeY = 100;
  const size_t sizeZ = 100;

  float *cube = (float *)malloc(sizeX * sizeY * sizeZ * sizeof(float));

  if (rank == 0) {
    log_info(rank, "Initializing cube...");
    initialize_cube(cube, sizeX, sizeY, sizeZ);
    log_info(rank, "Cube initialized.");
  }

  MPI_Finalize();
  return 0;
}
