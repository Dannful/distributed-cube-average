#include <math.h>
#include <mpi.h>
#include <string.h>

#include "coordinator.h"
#include "setup.h"
#include "worker.h"

problem_data_t dc_initialize_problem(MPI_Comm comm,
                                     unsigned int topology[DIMENSIONS],
                                     unsigned int border, unsigned int workers,
                                     dc_arguments_t arguments) {
  problem_data_t result;
  result.communicator = comm;
  result.num_workers = workers;
  result.iterations = ceil(arguments.time_max / arguments.dt);
  result.size_x =
      arguments.size_x + 2 * arguments.absorption_size + 2 * STENCIL;
  result.size_y =
      arguments.size_y + 2 * arguments.absorption_size + 2 * STENCIL;
  result.size_z =
      arguments.size_z + 2 * arguments.absorption_size + 2 * STENCIL;
  result.pp = (float *)calloc(result.size_x * result.size_y * result.size_z,
                              sizeof(float));
  result.pc = (float *)calloc(result.size_x * result.size_y * result.size_z,
                              sizeof(float));
  result.qp = (float *)calloc(result.size_x * result.size_y * result.size_z,
                              sizeof(float));
  result.qc = (float *)calloc(result.size_x * result.size_y * result.size_z,
                              sizeof(float));
  result.pp_workers = (float **)calloc(workers, sizeof(float *));
  result.pc_workers = (float **)calloc(workers, sizeof(float *));
  result.qp_workers = (float **)calloc(workers, sizeof(float *));
  result.qc_workers = (float **)calloc(workers, sizeof(float *));
  result.source_index = (int *)malloc(workers * sizeof(int));
  for(int i = 0; i < workers; i++) {
    result.source_index[i] = -1;
  }
  memcpy(result.topology, topology, DIMENSIONS * sizeof(unsigned int));
  result.worker_sizes = calloc(workers, sizeof(size_t *));
  return result;
}

void dc_partition_cube(problem_data_t *problem_data) {
  size_t partition_size_x =
      (problem_data->size_x - 2 * STENCIL) / problem_data->topology[0];
  size_t partition_size_y =
      (problem_data->size_y - 2 * STENCIL) / problem_data->topology[1];
  size_t partition_size_z =
      (problem_data->size_z - 2 * STENCIL) / problem_data->topology[2];
  size_t remainder_x = problem_data->size_x % problem_data->topology[0];
  size_t remainder_y = problem_data->size_y % problem_data->topology[1];
  size_t remainder_z = problem_data->size_z % problem_data->topology[2];
  for (size_t worker_z = 0; worker_z < problem_data->topology[2]; worker_z++) {
    for (size_t worker_y = 0; worker_y < problem_data->topology[1];
         worker_y++) {
      for (size_t worker_x = 0; worker_x < problem_data->topology[0];
           worker_x++) {
        int process_coordinates[DIMENSIONS] = {(int)worker_x, (int)worker_y,
                                               (int)worker_z};
        int worker;
        MPI_Cart_rank(problem_data->communicator, process_coordinates, &worker);
        size_t worker_size_x = worker_x == problem_data->topology[0] - 1
                                   ? partition_size_x + remainder_x
                                   : partition_size_x;
        size_t worker_size_y = worker_y == problem_data->topology[1] - 1
                                   ? partition_size_y + remainder_y
                                   : partition_size_y;
        size_t worker_size_z = worker_z == problem_data->topology[2] - 1
                                   ? partition_size_z + remainder_z
                                   : partition_size_z;
        worker_size_x += 2 * STENCIL;
        worker_size_y += 2 * STENCIL;
        worker_size_z += 2 * STENCIL;
        problem_data->worker_sizes[worker] =
            malloc(sizeof(size_t) * DIMENSIONS);
        problem_data->worker_sizes[worker][0] = worker_size_x;
        problem_data->worker_sizes[worker][1] = worker_size_y;
        problem_data->worker_sizes[worker][2] = worker_size_z;
        problem_data->pp_workers[worker] = (float *)malloc(
            sizeof(float) * worker_size_x * worker_size_y * worker_size_z);
        problem_data->pc_workers[worker] = (float *)malloc(
            sizeof(float) * worker_size_x * worker_size_y * worker_size_z);
        problem_data->qp_workers[worker] = (float *)malloc(
            sizeof(float) * worker_size_x * worker_size_y * worker_size_z);
        problem_data->qc_workers[worker] = (float *)malloc(
            sizeof(float) * worker_size_x * worker_size_y * worker_size_z);
        size_t count = 0;
        for (size_t z = process_coordinates[2] * partition_size_z;
             z < process_coordinates[2] * partition_size_z + worker_size_z;
             z++) {
          for (size_t y = process_coordinates[1] * partition_size_y;
               y < process_coordinates[1] * partition_size_y + worker_size_y;
               y++) {
            for (size_t x =
                     process_coordinates[0] * partition_size_x;
                 x < process_coordinates[0] * partition_size_x + worker_size_x;
                 x++) {
              size_t index = dc_get_index_for_coordinates(
                  x, y, z, problem_data->size_x, problem_data->size_y,
                  problem_data->size_z);
              size_t local_x =
                  x - process_coordinates[0] * partition_size_x;
              size_t local_y =
                  y - process_coordinates[1] * partition_size_y;
              size_t local_z =
                  z - process_coordinates[2] * partition_size_z;
              if (x == problem_data->size_x / 2 &&
                  y == problem_data->size_y / 2 &&
                  z == problem_data->size_z / 2 && local_x >= STENCIL &&
                  local_y >= STENCIL && local_z >= STENCIL &&
                  local_x < worker_size_x - STENCIL &&
                  local_y < worker_size_y - STENCIL &&
                  local_z < worker_size_z - STENCIL) {
                problem_data->source_index[worker] = count;
              }
              problem_data->pp_workers[worker][count] = problem_data->pp[index];
              problem_data->pc_workers[worker][count] = problem_data->pc[index];
              problem_data->qp_workers[worker][count] = problem_data->qp[index];
              problem_data->qc_workers[worker][count] = problem_data->qc[index];
              count++;
            }
          }
        }
      }
    }
  }
}

void dc_send_data_to_workers(problem_data_t problem_data) {
  for (size_t worker = 0; worker < problem_data.num_workers; worker++) {
    if (worker == COORDINATOR)
      continue;
    MPI_Send(&problem_data.source_index[worker], 1, MPI_INT, worker, 0,
             problem_data.communicator);
    MPI_Send(&problem_data.iterations, 1, MPI_UINT32_T, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.worker_sizes[worker], DIMENSIONS, MPI_UNSIGNED_LONG,
             worker, 0, problem_data.communicator);
    size_t count =
        dc_compute_count_from_sizes(problem_data.worker_sizes[worker]);
    MPI_Send(problem_data.pp_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.pc_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.qp_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.qc_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
  }
}

dc_result_t dc_receive_data_from_workers(dc_process_t coordinator_process,
                                         size_t cube_size_x, size_t cube_size_y,
                                         size_t cube_size_z) {
  dc_result_t result;
  result.pc =
      (float *)malloc(sizeof(float) * cube_size_x * cube_size_y * cube_size_z);
  result.qc =
      (float *)malloc(sizeof(float) * cube_size_x * cube_size_y * cube_size_z);
  size_t size_x = cube_size_x / coordinator_process.topology[0];
  size_t size_y = cube_size_y / coordinator_process.topology[1];
  size_t size_z = cube_size_z / coordinator_process.topology[2];
  for (int worker_x = 0; worker_x < coordinator_process.topology[0];
       worker_x++) {
    for (int worker_y = 0; worker_y < coordinator_process.topology[1];
         worker_y++) {
      for (int worker_z = 0; worker_z < coordinator_process.topology[2];
           worker_z++) {
        int worker_coordinates[DIMENSIONS] = {worker_x, worker_y, worker_z};
        int worker_rank;
        MPI_Cart_rank(coordinator_process.communicator, worker_coordinates,
                      &worker_rank);
        size_t *worker_sizes = NULL;
        float *pc = NULL;
        float *qc = NULL;
        size_t worker_count = 0;
        if (worker_rank == COORDINATOR) {
          worker_sizes = coordinator_process.sizes;
          pc = coordinator_process.pc;
          qc = coordinator_process.qc;
          worker_count = dc_compute_count_from_sizes(worker_sizes);
        } else {
          worker_sizes = malloc(DIMENSIONS * sizeof(size_t));
          MPI_Recv(worker_sizes, DIMENSIONS, MPI_UNSIGNED_LONG, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
          worker_count = dc_compute_count_from_sizes(worker_sizes);
          pc = malloc(sizeof(float) * worker_count);
          qc = malloc(sizeof(float) * worker_count);
          MPI_Recv(pc, worker_count, MPI_FLOAT, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
          MPI_Recv(qc, worker_count, MPI_FLOAT, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
        }
        for (size_t z = STENCIL; z < worker_sizes[2] - STENCIL; z++) {
          for (size_t y = STENCIL; y < worker_sizes[1] - STENCIL; y++) {
            for (size_t x = STENCIL; x < worker_sizes[0] - STENCIL; x++) {
              size_t worker_index = dc_get_index_for_coordinates(
                  x, y, z, worker_sizes[0], worker_sizes[1], worker_sizes[2]);
              size_t global_index = dc_get_index_for_coordinates(
                  worker_x * (size_x - 2 * STENCIL) + x,
                  worker_y * (size_y - 2 * STENCIL) + y,
                  worker_z * (size_z - 2 * STENCIL) + z, cube_size_x,
                  cube_size_y, cube_size_z);
              result.pc[global_index] = pc[worker_index];
              result.qc[global_index] = qc[worker_index];
            }
          }
        }
        if (worker_rank != COORDINATOR) {
          free(worker_sizes);
          free(pc);
          free(qc);
        }
      }
    }
  }
  return result;
}

void dc_free_problem_data_mem(problem_data_t *problem_data) {
  for (size_t i = 0; i < problem_data->num_workers; i++) {
    free(problem_data->pp_workers[i]);
    free(problem_data->pc_workers[i]);
    free(problem_data->qp_workers[i]);
    free(problem_data->qc_workers[i]);
    free(problem_data->worker_sizes[i]);
  }
  free(problem_data->source_index);
  free(problem_data->pp_workers);
  free(problem_data->pc_workers);
  free(problem_data->qp_workers);
  free(problem_data->qc_workers);
  free(problem_data->worker_sizes);

  problem_data->pp_workers = NULL;
  problem_data->qp_workers = NULL;
  problem_data->pc_workers = NULL;
  problem_data->qc_workers = NULL;
  problem_data->worker_sizes = NULL;
}
