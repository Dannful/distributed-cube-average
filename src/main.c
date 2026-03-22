#include <argp.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "boundary.h"
#include "coordinator.h"
#include "indexing.h"
#include "log.h"
#include "precomp.h"
#include "setup.h"
#include "worker.h"

static struct argp_option options[] = {
    {"size-x", 128, "INTEGER", 0, "Grid size in X"},
    {"size-y", 129, "INTEGER", 0, "Grid size in Y"},
    {"size-z", 130, "INTEGER", 0, "Grid size in Z"},
    {"dx", 131, "FLOAT", 0, "Step size in X"},
    {"dy", 132, "FLOAT", 0, "Step size in Y"},
    {"dz", 133, "FLOAT", 0, "Step size in Z"},
    {"dt", 134, "FLOAT", 0, "Step size in time"},
    {"time-max", 't', "FLOAT", 0, "Max time"},
    {"absorption", 'a', "INTEGER", 0, "Absorption zone size"},
    {"output-file", 'o', "PATH", 0,
     "Path to the file to output the results to"},
    {0},
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  dc_arguments_t *arguments = state->input;
  switch (key) {
  case 128:
    arguments->size_x = atoi(arg);
    break;
  case 129:
    arguments->size_y = atoi(arg);
    break;
  case 130:
    arguments->size_z = atoi(arg);
    break;
  case 131:
    arguments->dx = atof(arg);
    break;
  case 132:
    arguments->dy = atof(arg);
    break;
  case 133:
    arguments->dz = atof(arg);
    break;
  case 134:
    arguments->dt = atof(arg);
    break;
  case 't':
    arguments->time_max = atof(arg);
    break;
  case 'a':
    arguments->absorption_size = atoi(arg);
    break;
  case 'o':
    arguments->output_file = strdup(arg);
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
  MPI_Init(&argc, &argv);
  dc_arguments_t arguments = {0};
  argp_parse(&argp, argc, argv, 0, 0, &arguments);
  MPI_Comm communicator;
  int topology[DIMENSIONS] = {0};
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, DIMENSIONS, topology);
  dc_mpi_world_init(&communicator, topology);
  MPI_Comm_rank(communicator, &rank);

  const size_t sx =
      arguments.size_x + 2 * arguments.absorption_size + 2 * STENCIL;
  const size_t sy =
      arguments.size_y + 2 * arguments.absorption_size + 2 * STENCIL;
  const size_t sz =
      arguments.size_z + 2 * arguments.absorption_size + 2 * STENCIL;

  dc_process_t mpi_process =
      dc_process_init(communicator, rank, size, topology, sx, sy, sz,
                      arguments.dx, arguments.dy, arguments.dz, arguments.dt);

  if (rank == COORDINATOR) {
    dc_log_info(rank, "Distributing partition info to workers...");
    dc_distribute_partition_info(communicator, (unsigned int *)topology,
                                 arguments, size);

    size_t partition_size_x = (sx - 2 * STENCIL) / topology[0];
    size_t partition_size_y = (sy - 2 * STENCIL) / topology[1];
    size_t partition_size_z = (sz - 2 * STENCIL) / topology[2];
    size_t remainder_x = (sx - 2 * STENCIL) % topology[0];
    size_t remainder_y = (sy - 2 * STENCIL) % topology[1];
    size_t remainder_z = (sz - 2 * STENCIL) % topology[2];

    // Coordinator is always at position (0,0,0)
    mpi_process.sizes[0] = partition_size_x + 2 * STENCIL;
    mpi_process.sizes[1] = partition_size_y + 2 * STENCIL;
    mpi_process.sizes[2] = partition_size_z + 2 * STENCIL;

    // Handle remainder for last process in each dimension
    // Coordinator at (0,0,0) doesn't get remainder unless it's also the last
    if (topology[0] == 1)
      mpi_process.sizes[0] += remainder_x;
    if (topology[1] == 1)
      mpi_process.sizes[1] += remainder_y;
    if (topology[2] == 1)
      mpi_process.sizes[2] += remainder_z;

    size_t count = dc_compute_count_from_sizes(mpi_process.sizes);
    mpi_process.iterations = ceil(arguments.time_max / arguments.dt);

    // Check if source is in coordinator's partition
    size_t source_x, source_y, source_z;
    dc_determine_source(sx, sy, sz, &source_x, &source_y, &source_z);
    if (source_x < mpi_process.sizes[0] && source_y < mpi_process.sizes[1] &&
        source_z < mpi_process.sizes[2]) {
      mpi_process.source_index = (int)dc_get_index_for_coordinates(
          source_x, source_y, source_z, mpi_process.sizes[0],
          mpi_process.sizes[1], mpi_process.sizes[2]);
    } else {
      mpi_process.source_index = -1;
    }

    mpi_process.pp = (float *)calloc(count, sizeof(float));
    mpi_process.pc = (float *)calloc(count, sizeof(float));
    mpi_process.qp = (float *)calloc(count, sizeof(float));
    mpi_process.qc = (float *)calloc(count, sizeof(float));
    if (mpi_process.pp == NULL || mpi_process.pc == NULL ||
        mpi_process.qp == NULL || mpi_process.qc == NULL) {
      dc_log_error(rank, "OOM: could not allocate field arrays");
      MPI_Finalize();
      exit(1);
    }

    // Compute anisotropy and precomp vars for coordinator's partition
    mpi_process.anisotropy_vars = dc_compute_anisotropy_vars(
        mpi_process.sizes[0], mpi_process.sizes[1], mpi_process.sizes[2]);
    unsigned int seed = 0;
    // Coordinator is at global position (0,0,0)
    randomVelocityBoundaryPartition(
        mpi_process.sizes[0], mpi_process.sizes[1],
        mpi_process.sizes[2], // Local sizes
        sx, sy, sz,           // Global sizes
        0, 0, 0,              // Start coords (coordinator at origin)
        arguments.size_x, arguments.size_y, arguments.size_z, // Problem sizes
        STENCIL, arguments.absorption_size, mpi_process.anisotropy_vars.vpz,
        mpi_process.anisotropy_vars.vsv, &seed);
    mpi_process.precomp_vars = dc_compute_precomp_vars(
        mpi_process.sizes[0], mpi_process.sizes[1], mpi_process.sizes[2],
        mpi_process.anisotropy_vars);

    dc_log_info(
        rank, "Coordinator initialized locally with sizes %zu x %zu x %zu",
        mpi_process.sizes[0], mpi_process.sizes[1], mpi_process.sizes[2]);
  } else {
    dc_log_info(rank, "Awaiting partition info from coordinator...");
    dc_worker_init_from_partition_info(&mpi_process, communicator);
    dc_log_info(rank, "Local initialization complete.");
  }
  dc_log_info(rank, "Starting worker process...");
  double start_time = MPI_Wtime();
  double msamples_per_s = dc_worker_process(&mpi_process, communicator);
  double end_time = MPI_Wtime();
  double total_time = end_time - start_time;
  dc_send_data_to_coordinator(mpi_process, communicator);
  if (rank == COORDINATOR) {
    dc_receive_and_write_results(mpi_process, communicator, sx, sy, sz,
                                 arguments.output_file);
  }
  dc_worker_free(mpi_process);

  free(arguments.output_file);
  dc_free_anisotropy_vars(&mpi_process.anisotropy_vars);
  dc_free_precomp_vars(&mpi_process.precomp_vars);

  if (rank == COORDINATOR) {
    printf("rank,total_time,msamples_per_s\n");
  }
  MPI_Barrier(communicator);
  printf("%d,%lf,%lf\n", rank, total_time, msamples_per_s);

  if (rank == COORDINATOR) {
    size_t global_compute_x = sx - 2 * STENCIL;
    size_t global_compute_y = sy - 2 * STENCIL;
    size_t global_compute_z = sz - 2 * STENCIL;
    double global_msamples = ((double)global_compute_x * global_compute_y *
                              global_compute_z * mpi_process.iterations) /
                             1000000.0;
    double global_msamples_per_s = global_msamples / total_time;
    printf("*,%lf,%lf\n", total_time, global_msamples_per_s);
  }

  MPI_Finalize();
  return 0;
}
