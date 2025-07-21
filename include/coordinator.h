#pragma once

#include "setup.h"
#include <mpi.h>
#include <stdlib.h>

#define COORDINATOR 0

typedef struct {
  MPI_Comm communicator;
  unsigned int iterations;
  unsigned int stencil_size;
  size_t size_x;
  size_t size_y;
  size_t size_z;
  size_t num_workers;
  float *cube;
  float **workers;
  size_t **worker_indices;
  size_t *worker_count;
  unsigned int topology[DIMENSIONS];
} problem_data_t;

problem_data_t init_problem_data(MPI_Comm comm, unsigned int topology[DIMENSIONS], unsigned int workers, unsigned int iterations, unsigned int stencil_size, size_t size_x, size_t size_y, size_t size_z);
void send_data_to_workers(problem_data_t problem_data);
void partition_cube(problem_data_t problem_data);
void free_problem_data(problem_data_t problem_data);
