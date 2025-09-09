#include <argp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coordinator.h"
#include "log.h"
#include "setup.h"
#include "worker.h"
#include "precomp.h"
#include "boundary.h"

static struct argp_option options[] = {
    {"size_x", 128, "INTEGER", 0, "Grid size in X"},
    {"size_y", 129, "INTEGER", 0, "Grid size in Y"},
    {"size_z", 130, "INTEGER", 0, "Grid size in Z"},
    {"dx", 131, "FLOAT", 0, "Step size in X"},
    {"dy", 132, "FLOAT", 0, "Step size in Y"},
    {"dz", 133, "FLOAT", 0, "Step size in Z"},
    {"dt", 134, "FLOAT", 0, "Step size in time"},
    {"time_max", 't', "FLOAT", 0, "Max time"},
    {"stencil_size", 's', "INTEGER", 0, "Stencil size"},
    {"absorption", 'a', "INTEGER", 0, "Absorption zone size"}};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  dc_arguments_t *arguments = state->input;
  switch (key) {
  case 130:
    arguments->size_x = atoi(arg);
    break;
  case 131:
    arguments->size_y = atoi(arg);
    break;
  case 132:
    arguments->size_z = atoi(arg);
    break;
  case 133:
    arguments->dx = atof(arg);
    break;
  case 134:
    arguments->dy = atof(arg);
    break;
  case 135:
    arguments->dz = atof(arg);
    break;
  case 136:
    arguments->dt = atof(arg);
    break;
  case 't':
    arguments->time_max = atof(arg);
    break;
  case 's':
    arguments->stencil_size = atoi(arg);
    break;
  case 'a':
    arguments->absorption_size = atoi(arg);
    break;
  case 'o':
    strcpy(arguments->output_file, arg);
    break;
  case ARGP_KEY_END:
    if (arguments->size_x == 0 || arguments->size_y == 0 ||
        arguments->size_z == 0 || arguments->dx == 0 || arguments->dy == 0 ||
        arguments->dz == 0 || arguments->time_max == 0 || arguments->dt == 0 ||
        arguments->absorption_size == 0 || arguments->output_file == NULL) {
      argp_usage(state);
    }
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {
    options, parse_opt, NULL,
    "A program that solves Fletcher equations in a distributed setup"};

int main(int argc, char **argv) {
  dc_arguments_t arguments = {0};
  argp_parse(&argp, argc, argv, 0, 0, &arguments);
  MPI_Comm communicator;
  int topology[DIMENSIONS] = {0};
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, DIMENSIONS, topology);
  dc_mpi_world_init(&communicator, topology);
  MPI_Comm_rank(communicator, &rank);

  dc_process_t mpi_process = dc_process_init(communicator, rank, topology);

  double start_time = MPI_Wtime();

  const size_t border = 4;
  const size_t sx = arguments.size_x + 2 * arguments.absorption_size + 2 * border;
  const size_t sy = arguments.size_y + 2 * arguments.absorption_size + 2 * border;
  const size_t sz = arguments.size_z + 2 * arguments.absorption_size + 2 * border;

  dc_anisotropy_t anisotropy_vars = dc_compute_anisotropy_vars(sx, sy, sz);
  dc_precomp_vars precomp_vars = dc_compute_precomp_vars(sx, sy, sz, anisotropy_vars, border);

  randomVelocityBoundary(sx, sy, sz, arguments.size_x, arguments.size_y, arguments.size_z, border, arguments.absorption_size, anisotropy_vars.vpz, anisotropy_vars.vsv);

  if (rank == COORDINATOR) {
    dc_log_info(rank, "Initializing problem data...");
    problem_data_t problem_data = dc_initialize_problem(
        communicator, (unsigned int *)topology, size, arguments);
    dc_log_info(rank, "Partitioning cube...");
    dc_partition_cube(problem_data);
    dc_log_info(rank, "Partition completed. Sending data to workers...");
    dc_send_data_to_workers(problem_data);
    memmove(mpi_process.sizes, problem_data.worker_sizes[0],
            sizeof(size_t) * DIMENSIONS);
    mpi_process.data =
        malloc(sizeof(float) *
               dc_compute_count_from_sizes(problem_data.worker_sizes[0]));
    memmove(mpi_process.data, problem_data.workers[0],
            sizeof(float) *
                dc_compute_count_from_sizes(problem_data.worker_sizes[0]));
    mpi_process.stencil_size = problem_data.stencil_size;
    mpi_process.iterations = problem_data.iterations;
    mpi_process.rank = rank;
    dc_free_problem_data_mem(problem_data);
  } else {
    dc_log_info(rank, "Awaiting data from coordinator...");
    dc_worker_receive_data(&mpi_process);
    dc_log_info(rank, "Data received from coordinator.");
  }
  dc_log_info(rank, "Starting worker process...");
  dc_worker_process(mpi_process);
  dc_send_data_to_coordinator(mpi_process);
  if (rank == COORDINATOR) {
    size_t total_size = arguments.size_x * arguments.size_y * arguments.size_z;
    float *cube = dc_receive_data_from_workers(
        mpi_process, arguments.size_x, arguments.size_y, arguments.size_z);
    FILE *output = fopen(arguments.output_file, "wb");
    fwrite(cube, sizeof(float), total_size, output);
    fclose(output);
    free(cube);
  }
  dc_worker_free(mpi_process);

  double end_time = MPI_Wtime();

  free_anisotropy_vars(&anisotropy_vars);
  free_precomp_vars(&precomp_vars);

  dc_log_info(rank, "Elapsed time: %lf seconds", end_time - start_time);

  MPI_Finalize();
  return 0;
}
