#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/coordinator.h"
#include "../include/log.h"
#include "../include/io.h"
#include "../include/worker.h"

int main(int argc, char **argv) {
  int rank, size;
  char hostname[256];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const size_t size_x = 100;
  const size_t size_y = 100;
  const size_t size_z = 100;

  const unsigned int iterations = 6;

  size_t total_size_x;
  size_t total_size_y;
  size_t total_size_z;

  if (rank == 0) {
    float *cube = (float *)malloc(size_x * size_y * size_z * sizeof(float));
    float **partitions = (float **)calloc(size, sizeof(float *));
    size_t *worker_count = calloc(size, sizeof(size_t));
    size_t **worker_indices = calloc(size, sizeof(size_t));

    log_info(rank, "Initializing cube with values X = %d, Y = %d and Z = %d...", size_x, size_y, size_z);
    initialize_cube(cube, size_x, size_y, size_z);
    log_info(rank, "Cube initialized.");
    log_info(rank, "Initializing cube partitioning with %d workers.", size);
    partition_cube(partitions, worker_indices, worker_count, &total_size_x, &total_size_y, &total_size_z, size, cube, size_x, size_y, size_z);
    log_info(rank, "Partitioning succeeded.");
    for(int i = 1; i < size; i++) {
      log_info(rank, "Sending %lu values to worker %d...", worker_count[i], i);
      MPI_Ssend(&size_x, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Ssend(&size_y, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Ssend(&size_z, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Ssend(worker_count + i, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Ssend(worker_indices[i], worker_count[i], MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Ssend(partitions[i], worker_count[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD);
    }
    for(int i = 0; i < size; i++) {
      if(partitions[i] != NULL)
        free(partitions[i]);
      if(worker_indices[i] != NULL)
        free(worker_indices[i]);
    }
    free(partitions);
    free(worker_indices);
    free(worker_count);
  } else {
    size_t count;
    log_info(rank, "Waiting for data from coordinator...");
    MPI_Safe_Recv(&total_size_x, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
    MPI_Safe_Recv(&total_size_y, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
    MPI_Safe_Recv(&total_size_z, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
    MPI_Safe_Recv(&count, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
    size_t *indices = (size_t *)malloc(count * sizeof(size_t));
    MPI_Safe_Recv(indices, count, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
    float *partition = (float *)malloc(count * sizeof(float));
    MPI_Safe_Recv(partition, count, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
    log_info(rank, "Received %lu values from coordinator.", count);
    worker_process(iterations, partition, count, indices);
    free(partition);
    free(indices);
  }
  MPI_Finalize();
  return 0;
}
