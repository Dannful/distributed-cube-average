#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coordinator.h"
#include "log.h"
#include "setup.h"
#include "worker.h"

int main(int argc, char **argv) {
  MPI_Comm communicator;
  int topology[DIMENSIONS] = {0};
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, DIMENSIONS, topology);
  dc_mpi_world_init(&communicator, topology);
  MPI_Comm_rank(communicator, &rank);

  dc_process_t mpi_process = dc_process_init(communicator, rank, topology);

  const size_t size_x = 3;
  const size_t size_y = 3;
  const size_t size_z = 3;
  const unsigned int iterations = 6, stencil_size = 3;

  if (rank == COORDINATOR) {
    dc_log_info(rank, "Initializing problem data...");
    problem_data_t problem_data =
        dc_initialize_problem(communicator, (unsigned int *)topology, size,
                              iterations, stencil_size, size_x, size_y, size_z);
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
    float *cube =
        dc_receive_data_from_workers(mpi_process, size_x, size_y, size_z);
    // TODO: Export data somewhere
    free(cube);
  }
  dc_worker_free(mpi_process);
  MPI_Finalize();
  return 0;
}
