#pragma once

#include "precomp.h"
#include "setup.h"
#include <mpi.h>
#include <stdlib.h>

#define COORDINATOR 0

typedef struct {
  MPI_Comm communicator;
  unsigned int iterations;
  size_t size_x;
  size_t size_y;
  size_t size_z;
  size_t num_workers;
  float dt;
  float *pp, *pc, *qp, *qc;
  float **pp_workers;
  float **pc_workers;
  float **qp_workers;
  float **qc_workers;
  dc_precomp_vars **precomp_vars_workers;
  int *source_index;
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
  size_t absorption_size;
  char *output_file;
} dc_arguments_t;

typedef struct {
  float *pc, *qc;
} dc_result_t;

problem_data_t dc_initialize_problem(MPI_Comm comm,
                                     unsigned int topology[DIMENSIONS],
                                     unsigned int border, unsigned int workers,
                                     dc_arguments_t arguments);
void dc_send_data_to_workers(problem_data_t problem_data);
void dc_partition_cube(problem_data_t *problem_data,
                         dc_precomp_vars precomp_vars);
void dc_free_problem_data_mem(problem_data_t *problem_data);
void dc_determine_source(size_t size_x, size_t size_y, size_t size_z,
                         size_t *source_x, size_t *source_y, size_t *source_z);
dc_result_t dc_receive_data_from_workers(dc_process_t coordinator_process,
                                         size_t cube_size_x, size_t cube_size_y,
                                         size_t cube_size_z);
