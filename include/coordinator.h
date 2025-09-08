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
  float *vpz;
  float *vsv;
  float *epsilon;
  float *delta;
  float *phi;
  float *theta;
  float *cube;
  float **workers;
  size_t **worker_sizes;
  unsigned int topology[DIMENSIONS];
} problem_data_t;

typedef struct {
  size_t size_x;
  size_t size_y;
  size_t size_z;
  float dx;
  float dy;
  float dz;
  float dt;
  float time_max;
  size_t stencil_size;
  size_t absorption_size;
} dc_arguments_t;

problem_data_t dc_initialize_problem(MPI_Comm comm, unsigned int topology[DIMENSIONS], unsigned int workers, dc_arguments_t arguments);
void dc_send_data_to_workers(problem_data_t problem_data);
void dc_partition_cube(problem_data_t problem_data);
void dc_free_problem_data_mem(problem_data_t problem_data);
float *dc_receive_data_from_workers(dc_process_t coordinator_process, size_t cube_size_x, size_t cube_size_y, size_t cube_size_z);
