#include <math.h>
#include <mpi.h>
#include <string.h>

#include "coordinator.h"
#include "log.h"
#include "setup.h"
#include "worker.h"

problem_data_t dc_initialize_problem(MPI_Comm comm,
                                     unsigned int topology[DIMENSIONS],
                                     unsigned int border, unsigned int workers,
                                     dc_arguments_t arguments) {
  problem_data_t result;
  result.communicator = comm;
  result.num_workers = workers;
  result.dt = arguments.dt;
  result.iterations = ceil(arguments.time_max / arguments.dt);
  result.size_x =
      arguments.size_x + 2 * arguments.absorption_size + 2 * STENCIL;
  result.size_y =
      arguments.size_y + 2 * arguments.absorption_size + 2 * STENCIL;
  result.size_z =
      arguments.size_z + 2 * arguments.absorption_size + 2 * STENCIL;
  result.pp = (float *)malloc(result.size_x * result.size_y * result.size_z *
                              sizeof(float));
  result.pc = (float *)malloc(result.size_x * result.size_y * result.size_z *
                              sizeof(float));
  result.qp = (float *)malloc(result.size_x * result.size_y * result.size_z *
                              sizeof(float));
  result.qc = (float *)malloc(result.size_x * result.size_y * result.size_z *
                              sizeof(float));
  for (int i = 0; i < result.size_x * result.size_y * result.size_z; i++) {
    result.pp[i] = 0.0f;
    result.pc[i] = 0.0f;
    result.qp[i] = 0.0f;
    result.qc[i] = 0.0f;
  }
  result.worker_indices = (size_t **)calloc(workers, sizeof(size_t *));
  result.pp_workers = (float **)calloc(workers, sizeof(float *));
  result.pc_workers = (float **)calloc(workers, sizeof(float *));
  result.qp_workers = (float **)calloc(workers, sizeof(float *));
  result.qc_workers = (float **)calloc(workers, sizeof(float *));
  result.source_index = (int *)malloc(workers * sizeof(int));
  for (int i = 0; i < workers; i++) {
    result.source_index[i] = -1;
  }
  memcpy(result.topology, topology, DIMENSIONS * sizeof(unsigned int));
  result.worker_sizes = calloc(workers, sizeof(size_t *));
  return result;
}

void dc_determine_source(size_t size_x, size_t size_y, size_t size_z,
                         unsigned int topology[DIMENSIONS], size_t *source_x,
                         size_t *source_y, size_t *source_z) {
  size_t partition_size_x = STENCIL + (size_x - 2 * STENCIL) / topology[0];
  size_t partition_size_y = STENCIL + (size_y - 2 * STENCIL) / topology[1];
  size_t partition_size_z = STENCIL + (size_z - 2 * STENCIL) / topology[2];

  dc_log_info(COORDINATOR,
              "Initializing source search with partition sizes (%d, %d, %d)...",
              partition_size_x, partition_size_y, partition_size_z);

  *source_x = size_x / 4;
  *source_y = size_y / 4;
  *source_z = size_z / 4;

  size_t x = 0, y = 0, z = 0;
  while (1) {
    size_t new_source_x = *source_x + x;
    size_t new_source_y = *source_y + y;
    size_t new_source_z = *source_z + z;
    size_t remainder_x = new_source_x % partition_size_x;
    size_t remainder_y = new_source_y % partition_size_y;
    size_t remainder_z = new_source_z % partition_size_z;

    if (remainder_x >= STENCIL && remainder_x < partition_size_x - STENCIL &&
        remainder_y >= STENCIL && remainder_y < partition_size_y - STENCIL &&
        remainder_z >= STENCIL && remainder_z < partition_size_z - STENCIL)
      break;
    int step_x = remainder_x == 0 ? 1 : (remainder_x < STENCIL ? 1 : -1);
    int step_y = remainder_y == 0 ? 1 : (remainder_y < STENCIL ? 1 : -1);
    int step_z = remainder_z == 0 ? 1 : (remainder_z < STENCIL ? 1 : -1);
    step_x *=
        remainder_x < STENCIL || remainder_x >= partition_size_x - STENCIL;
    step_y *=
        remainder_y < STENCIL || remainder_y >= partition_size_y - STENCIL;
    step_z *=
        remainder_z < STENCIL || remainder_z >= partition_size_z - STENCIL;
    dc_log_info(COORDINATOR, "Stepping (%d, %d, %d) towards new source.",
                step_x, step_y, step_z);
    x += step_x;
    y += step_y;
    z += step_z;
  }
  *source_x += x;
  *source_y += y;
  *source_z += z;
}

void dc_partition_cube(problem_data_t *problem_data) {
  size_t partition_size_x =
      (problem_data->size_x - 2 * STENCIL) / problem_data->topology[0];
  size_t partition_size_y =
      (problem_data->size_y - 2 * STENCIL) / problem_data->topology[1];
  size_t partition_size_z =
      (problem_data->size_z - 2 * STENCIL) / problem_data->topology[2];
  size_t remainder_x =
      (problem_data->size_x - 2 * STENCIL) % problem_data->topology[0];
  size_t remainder_y =
      (problem_data->size_y - 2 * STENCIL) % problem_data->topology[1];
  size_t remainder_z =
      (problem_data->size_z - 2 * STENCIL) % problem_data->topology[2];
  dc_log_info(COORDINATOR, "Partition sizes: X: %d, Y: %d, Z: %d",
              partition_size_x, partition_size_y, partition_size_z);
  dc_log_info(COORDINATOR, "Partition remainders: X: %d, Y: %d, Z: %d",
              remainder_x, remainder_y, remainder_z);
  size_t source_x, source_y, source_z;
  dc_log_info(COORDINATOR, "Calculating source...");
  dc_determine_source(problem_data->size_x, problem_data->size_y,
                      problem_data->size_z, problem_data->topology, &source_x,
                      &source_y, &source_z);
  dc_log_info(COORDINATOR,
              "Source is located at coordinates X: %d, Y: %d and Z: %d",
              source_x, source_y, source_z);
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
        problem_data->worker_indices[worker] = (size_t *)malloc(
            sizeof(size_t) * worker_size_x * worker_size_y * worker_size_z);
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
            for (size_t x = process_coordinates[0] * partition_size_x;
                 x < process_coordinates[0] * partition_size_x + worker_size_x;
                 x++) {
              size_t index = dc_get_index_for_coordinates(
                  x, y, z, problem_data->size_x, problem_data->size_y,
                  problem_data->size_z);
              problem_data->pp_workers[worker][count] = problem_data->pp[index];
              problem_data->pc_workers[worker][count] = problem_data->pc[index];
              problem_data->qp_workers[worker][count] = problem_data->qp[index];
              problem_data->qc_workers[worker][count] = problem_data->qc[index];
              problem_data->worker_indices[worker][count] = index;
              if (x == source_x && y == source_y && z == source_z) {
                problem_data->source_index[worker] = count;
                dc_log_info(COORDINATOR,
                            "Worker %d at coordinates (%d, %d, %d) will handle "
                            "the source at global index %d and local index %d",
                            worker, process_coordinates[0],
                            process_coordinates[1], process_coordinates[2],
                            index, count);
              }
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
    MPI_Send(problem_data.worker_indices[worker], count, MPI_UNSIGNED_LONG,
             worker, 0, problem_data.communicator);
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
      (float *)malloc(cube_size_x * cube_size_y * cube_size_z * sizeof(float));
  result.qc =
      (float *)malloc(cube_size_x * cube_size_y * cube_size_z * sizeof(float));
  for (int i = 0; i < cube_size_x * cube_size_y * cube_size_z; i++) {
    result.pc[i] = 0.0f;
    result.qc[i] = 0.0f;
  }
  size_t size_x = (cube_size_x - 2 * STENCIL) / coordinator_process.topology[0];
  size_t size_y = (cube_size_y - 2 * STENCIL) / coordinator_process.topology[1];
  size_t size_z = (cube_size_z - 2 * STENCIL) / coordinator_process.topology[2];
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
        size_t worker_sizes[DIMENSIONS] = {0};
        size_t *indices = NULL;
        float *pc = NULL;
        float *qc = NULL;
        size_t worker_count = 0;
        if (worker_rank == COORDINATOR) {
          memcpy(worker_sizes, coordinator_process.sizes,
                 DIMENSIONS * sizeof(size_t));
          indices = coordinator_process.indices;
          pc = coordinator_process.pc;
          qc = coordinator_process.qc;
          worker_count = dc_compute_count_from_sizes(worker_sizes);
        } else {
          MPI_Recv(worker_sizes, DIMENSIONS, MPI_UNSIGNED_LONG, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
          worker_count = dc_compute_count_from_sizes(worker_sizes);
          indices = malloc(sizeof(size_t) * worker_count);
          pc = malloc(sizeof(float) * worker_count);
          qc = malloc(sizeof(float) * worker_count);
          MPI_Recv(indices, worker_count, MPI_UNSIGNED_LONG, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
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
              result.pc[indices[worker_index]] = pc[worker_index];
              result.qc[indices[worker_index]] = qc[worker_index];
            }
          }
        }
        // for (size_t z = STENCIL; z < worker_sizes[2] - STENCIL; z++) {
        //   for (size_t y = STENCIL; y < worker_sizes[1] - STENCIL; y++) {
        //     for (size_t x = STENCIL; x < worker_sizes[0] - STENCIL; x++) {
        //       size_t local_x = x - STENCIL;
        //       size_t local_y = y - STENCIL;
        //       size_t local_z = z - STENCIL;
        //       size_t worker_index = dc_get_index_for_coordinates(
        //           x, y, z, worker_sizes[0], worker_sizes[1],
        //           worker_sizes[2]);
        //       size_t global_sizes[DIMENSIONS] = {cube_size_x, cube_size_y,
        //                                          cube_size_z};
        //       size_t local_coordinates[DIMENSIONS] = {x, y, z};
        //       size_t global_index = dc_get_global_coordinates(
        //           worker_coordinates, worker_sizes, global_sizes,
        //           local_coordinates, coordinator_process.topology);
        //       global_index = dc_get_index_for_coordinates(
        //           STENCIL + worker_coordinates[0] * size_x + local_x,
        //           STENCIL + worker_coordinates[1] * size_y + local_y,
        //           STENCIL + worker_coordinates[2] * size_z + local_z,
        //           global_sizes[0], global_sizes[1], global_sizes[2]);
        //       result.pc[global_index] = pc[worker_index];
        //       result.qc[global_index] = qc[worker_index];
        //     }
        //   }
        // }
        if (worker_rank != COORDINATOR) {
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
    free(problem_data->worker_indices[i]);
  }
  free(problem_data->source_index);
  free(problem_data->pp_workers);
  free(problem_data->pc_workers);
  free(problem_data->qp_workers);
  free(problem_data->qc_workers);
  free(problem_data->worker_sizes);
  free(problem_data->worker_indices);

  problem_data->pp_workers = NULL;
  problem_data->qp_workers = NULL;
  problem_data->pc_workers = NULL;
  problem_data->qc_workers = NULL;
  problem_data->worker_sizes = NULL;
  problem_data->worker_indices = NULL;
}
