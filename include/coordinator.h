#pragma once

#include "dc_process.h"
#include <mpi.h>
#include <stdlib.h>

#define COORDINATOR 0

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
  size_t local_sizes[DIMENSIONS];
  size_t global_sizes[DIMENSIONS];
  size_t start_coords[DIMENSIONS];
  size_t problem_sizes[DIMENSIONS];
  unsigned int iterations;
  int source_index;
  size_t absorption_size;
} dc_partition_info_t;

void dc_determine_source(size_t size_x, size_t size_y, size_t size_z,
                         size_t *source_x, size_t *source_y, size_t *source_z);

void dc_distribute_partition_info(MPI_Comm comm, unsigned int *topology,
                                  dc_arguments_t arguments, size_t num_workers);

void dc_receive_and_write_results(dc_process_t coordinator_process,
                                  MPI_Comm comm, size_t global_sx,
                                  size_t global_sy, size_t global_sz,
                                  const char *output_file);
