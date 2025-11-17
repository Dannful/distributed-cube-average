#pragma once

#include <mpi.h>

#include "definitions.h"

typedef enum {
  LEFT = 0,
  UP = 1,
  FRONT = 2,
  RIGHT = 3,
  DOWN = 4,
  BACK = 5
} neighbour_direction_t;

#include "dc_process.h"

void dc_mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]);
dc_process_t dc_process_init(MPI_Comm communicator, int rank,
                             int topology[DIMENSIONS], size_t sx, size_t sy,
                             size_t sz, float dx, float dy, float dz, float dt);
