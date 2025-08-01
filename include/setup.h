#pragma once

#include <mpi.h>

#define DIMENSIONS 3

typedef enum {
  LEFT = 0,
  UP = 1,
  FRONT = 2,
  RIGHT = 3,
  DOWN = 4,
  BACK = 5
} neighbour_direction_t;

typedef struct {
  MPI_Comm communicator;
  int rank;
  int coordinates[DIMENSIONS];
  int neighbours[DIMENSIONS * 2];
  unsigned int iterations;
  size_t count;
  unsigned int stencil_size;
  size_t *indices;
  size_t sizes[DIMENSIONS];
  float *data;
} mpi_process_t;

void mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]);
mpi_process_t mpi_process_init(MPI_Comm communicator, int rank, int topology[DIMENSIONS]);
