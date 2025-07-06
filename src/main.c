#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/coordinator.h"
#include "../include/log.h"

int main(int argc, char **argv) {
  int rank, size;
  char hostname[256];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const size_t size_x = 100;
  const size_t size_y = 100;
  const size_t size_z = 100;

  size_t total_size_x;
  size_t total_size_y;
  size_t total_size_z;

  float *cube = (float *)malloc(size_x * size_y * size_z * sizeof(float));
  float **partitions = (float **)calloc(size, sizeof(float *));

  if (rank == 0) {
    log_info(rank, "Initializing cube with values X = %d, Y = %d and Z = %d...", size_x, size_y, size_z);
    initialize_cube(cube, size_x, size_y, size_z);
    log_info(rank, "Cube initialized.");
    log_info(rank, "Initializing cube partitioning with %d workers.", size);
    partition_cube(partitions, &total_size_x, &total_size_y, &total_size_z, size, cube, size_x, size_y, size_z);
    log_info(rank, "Partitioning succeeded.");
    for(int i = 0; i < size; i++) {
      if(partitions[i] != NULL)
        free(partitions[i]);
    }
    free(partitions);
  }
  MPI_Finalize();
  return 0;
}
