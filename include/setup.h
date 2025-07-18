#pragma once

#include <mpi.h>

#define DIMENSIONS 3

typedef enum {
  LEFT = 0,
  UP = 0,
  FRONT = 0,
  RIGHT = 1,
  DOWN = 1,
  BACK = 1
} neighbour_direction_t;

typedef struct {
  MPI_Comm communicator;
  int rank;
  int coordinates[DIMENSIONS];
  int neighbours[DIMENSIONS][2];
} mpi_process_t;

typedef struct {
  unsigned int iterations;
  unsigned int stencil_size;
  size_t size_x;
  size_t size_y;
  size_t size_z;
  float *cube;
} problem_data_t;

void mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]);
mpi_process_t mpi_process_init(MPI_Comm communicator, const int topology[DIMENSIONS]);
problem_data_t init_problem_data(unsigned int iterations, unsigned int stencil_size, size_t size_x, size_t size_y, size_t size_z);
