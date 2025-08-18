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
  int topology[DIMENSIONS];
  unsigned int iterations;
  unsigned int stencil_size;
  size_t sizes[DIMENSIONS];
  float *data;
} dc_process_t;

void dc_mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]);
dc_process_t dc_process_init(MPI_Comm communicator, int rank, int topology[DIMENSIONS]);
