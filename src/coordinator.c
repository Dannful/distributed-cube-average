#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include "coordinator.h"
#include "log.h"
#include "precomp.h"
#include "setup.h"
#include "stdlib.h"
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
  if (result.pp == NULL) {
    dc_log_error(
        COORDINATOR,
        "OOM: could not allocate memory for pp in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.pc = (float *)malloc(result.size_x * result.size_y * result.size_z *
                              sizeof(float));
  if (result.pc == NULL) {
    dc_log_error(
        COORDINATOR,
        "OOM: could not allocate memory for pc in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.qp = (float *)malloc(result.size_x * result.size_y * result.size_z *
                              sizeof(float));
  if (result.qp == NULL) {
    dc_log_error(
        COORDINATOR,
        "OOM: could not allocate memory for qp in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.qc = (float *)malloc(result.size_x * result.size_y * result.size_z *
                              sizeof(float));
  if (result.qc == NULL) {
    dc_log_error(
        COORDINATOR,
        "OOM: could not allocate memory for qc in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  for (int i = 0; i < result.size_x * result.size_y * result.size_z; i++) {
    result.pp[i] = 0.0f;
    result.pc[i] = 0.0f;
    result.qp[i] = 0.0f;
    result.qc[i] = 0.0f;
  }
  result.pp_workers = (float **)calloc(workers, sizeof(float *));
  if (result.pp_workers == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for pp_workers "
                              "in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.pc_workers = (float **)calloc(workers, sizeof(float *));
  if (result.pc_workers == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for pc_workers "
                              "in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.qp_workers = (float **)calloc(workers, sizeof(float *));
  if (result.qp_workers == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for qp_workers "
                              "in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.qc_workers = (float **)calloc(workers, sizeof(float *));
  if (result.qc_workers == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for qc_workers "
                              "in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.source_index = (int *)malloc(workers * sizeof(int));
  if (result.source_index == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for source_index "
                              "in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  for (int i = 0; i < workers; i++) {
    result.source_index[i] = -1;
  }
  memcpy(result.topology, topology, DIMENSIONS * sizeof(unsigned int));
  result.worker_sizes = calloc(workers, sizeof(size_t *));
  if (result.worker_sizes == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for worker_sizes "
                              "in dc_initialize_problem");
    MPI_Finalize();
    exit(1);
  }
  result.precomp_vars_workers = calloc(workers, sizeof(dc_precomp_vars *));
  if (result.precomp_vars_workers == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for "
                              "precomp_vars_workers in dc_initialize_problem");
    exit(1);
  }
  return result;
}

void dc_determine_source(size_t size_x, size_t size_y, size_t size_z,
                         size_t *source_x, size_t *source_y, size_t *source_z) {
  *source_x = size_x / 2;
  *source_y = size_y / 2;
  *source_z = size_z / 2;
}

#define OOM_CHECK(ptr, var, worker)                                            \
  if (ptr == NULL) {                                                           \
    dc_log_error(COORDINATOR,                                                  \
                 "OOM: could not allocate memory for " var " for worker %d",   \
                 worker);                                                      \
    MPI_Finalize();                                                            \
    exit(1);                                                                   \
  }

void dc_partition_cube(problem_data_t *problem_data,
                       dc_precomp_vars precomp_vars) {
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
                      problem_data->size_z, &source_x, &source_y, &source_z);
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
        OOM_CHECK(problem_data->worker_sizes[worker], "worker_sizes", worker);
        problem_data->worker_sizes[worker][0] = worker_size_x;
        problem_data->worker_sizes[worker][1] = worker_size_y;
        problem_data->worker_sizes[worker][2] = worker_size_z;
        size_t worker_domain_size =
            worker_size_x * worker_size_y * worker_size_z;
        problem_data->pp_workers[worker] =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->pp_workers[worker], "pp_workers", worker);
        problem_data->pc_workers[worker] =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->pc_workers[worker], "pc_workers", worker);
        problem_data->qp_workers[worker] =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->qp_workers[worker], "qp_workers", worker);
        problem_data->qc_workers[worker] =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->qc_workers[worker], "qc_workers", worker);
        problem_data->precomp_vars_workers[worker] =
            (dc_precomp_vars *)malloc(sizeof(dc_precomp_vars));
        OOM_CHECK(problem_data->precomp_vars_workers[worker],
                  "precomp_vars_workers", worker);
        problem_data->precomp_vars_workers[worker]->ch1dxx =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->ch1dxx,
                  "precomp_vars_workers.ch1dxx", worker);
        problem_data->precomp_vars_workers[worker]->ch1dyy =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->ch1dyy,
                  "precomp_vars_workers.ch1dyy", worker);
        problem_data->precomp_vars_workers[worker]->ch1dzz =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->ch1dzz,
                  "precomp_vars_workers.ch1dzz", worker);
        problem_data->precomp_vars_workers[worker]->ch1dxy =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->ch1dxy,
                  "precomp_vars_workers.ch1dxy", worker);
        problem_data->precomp_vars_workers[worker]->ch1dyz =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->ch1dyz,
                  "precomp_vars_workers.ch1dyz", worker);
        problem_data->precomp_vars_workers[worker]->ch1dxz =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->ch1dxz,
                  "precomp_vars_workers.ch1dxz", worker);
        problem_data->precomp_vars_workers[worker]->v2px =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->v2px,
                  "precomp_vars_workers.v2px", worker);
        problem_data->precomp_vars_workers[worker]->v2pz =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->v2pz,
                  "precomp_vars_workers.v2pz", worker);
        problem_data->precomp_vars_workers[worker]->v2sz =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->v2sz,
                  "precomp_vars_workers.v2sz", worker);
        problem_data->precomp_vars_workers[worker]->v2pn =
            (float *)malloc(sizeof(float) * worker_domain_size);
        OOM_CHECK(problem_data->precomp_vars_workers[worker]->v2pn,
                  "precomp_vars_workers.v2pn", worker);
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
              problem_data->precomp_vars_workers[worker]->ch1dxx[count] =
                  precomp_vars.ch1dxx[index];
              problem_data->precomp_vars_workers[worker]->ch1dyy[count] =
                  precomp_vars.ch1dyy[index];
              problem_data->precomp_vars_workers[worker]->ch1dzz[count] =
                  precomp_vars.ch1dzz[index];
              problem_data->precomp_vars_workers[worker]->ch1dxy[count] =
                  precomp_vars.ch1dxy[index];
              problem_data->precomp_vars_workers[worker]->ch1dyz[count] =
                  precomp_vars.ch1dyz[index];
              problem_data->precomp_vars_workers[worker]->ch1dxz[count] =
                  precomp_vars.ch1dxz[index];
              problem_data->precomp_vars_workers[worker]->v2px[count] =
                  precomp_vars.v2px[index];
              problem_data->precomp_vars_workers[worker]->v2pz[count] =
                  precomp_vars.v2pz[index];
              problem_data->precomp_vars_workers[worker]->v2sz[count] =
                  precomp_vars.v2sz[index];
              problem_data->precomp_vars_workers[worker]->v2pn[count] =
                  precomp_vars.v2pn[index];
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
    MPI_Send(problem_data.pp_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.pc_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.qp_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.qc_workers[worker], count, MPI_FLOAT, worker, 0,
             problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->ch1dxx, count,
             MPI_FLOAT, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->ch1dyy, count,
             MPI_FLOAT, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->ch1dzz, count,
             MPI_FLOAT, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->ch1dxy, count,
             MPI_FLOAT, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->ch1dyz, count,
             MPI_FLOAT, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->ch1dxz, count,
             MPI_FLOAT, worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->v2px, count, MPI_FLOAT,
             worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->v2pz, count, MPI_FLOAT,
             worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->v2sz, count, MPI_FLOAT,
             worker, 0, problem_data.communicator);
    MPI_Send(problem_data.precomp_vars_workers[worker]->v2pn, count, MPI_FLOAT,
             worker, 0, problem_data.communicator);
  }
}

dc_result_t dc_receive_data_from_workers(dc_process_t coordinator_process,
                                         size_t cube_size_x, size_t cube_size_y,
                                         size_t cube_size_z) {
  dc_result_t result;
  result.pc =
      (float *)malloc(cube_size_x * cube_size_y * cube_size_z * sizeof(float));
  if (result.pc == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for result.pc in "
                              "dc_receive_data_from_workers");
    MPI_Finalize();
    exit(1);
  }
  result.qc =
      (float *)malloc(cube_size_x * cube_size_y * cube_size_z * sizeof(float));
  if (result.qc == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for result.qc in "
                              "dc_receive_data_from_workers");
    MPI_Finalize();
    exit(1);
  }
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
        float *pc = NULL;
        float *qc = NULL;
        size_t worker_count = 0;
        if (worker_rank == COORDINATOR) {
          memcpy(worker_sizes, coordinator_process.sizes,
                 DIMENSIONS * sizeof(size_t));
          pc = coordinator_process.pc;
          qc = coordinator_process.qc;
          worker_count = dc_compute_count_from_sizes(worker_sizes);
        } else {
          MPI_Recv(worker_sizes, DIMENSIONS, MPI_UNSIGNED_LONG, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
          worker_count = dc_compute_count_from_sizes(worker_sizes);
          pc = malloc(sizeof(float) * worker_count);
          if (pc == NULL) {
            dc_log_error(COORDINATOR,
                         "OOM: could not allocate memory for pc for worker %d "
                         "in dc_receive_data_from_workers",
                         worker_rank);
            MPI_Finalize();
            exit(1);
          }
          qc = malloc(sizeof(float) * worker_count);
          if (qc == NULL) {
            dc_log_error(COORDINATOR,
                         "OOM: could not allocate memory for qc for worker %d "
                         "in dc_receive_data_from_workers",
                         worker_rank);
            MPI_Finalize();
            exit(1);
          }
          MPI_Recv(pc, worker_count, MPI_FLOAT, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
          MPI_Recv(qc, worker_count, MPI_FLOAT, worker_rank, 0,
                   coordinator_process.communicator, MPI_STATUSES_IGNORE);
        }
        for (size_t z = STENCIL; z < worker_sizes[2] - STENCIL; z++) {
          for (size_t y = STENCIL; y < worker_sizes[1] - STENCIL; y++) {
            for (size_t x = STENCIL; x < worker_sizes[0] - STENCIL; x++) {
              size_t local_x = x - STENCIL;
              size_t local_y = y - STENCIL;
              size_t local_z = z - STENCIL;
              size_t worker_index = dc_get_index_for_coordinates(
                  x, y, z, worker_sizes[0], worker_sizes[1], worker_sizes[2]);
              size_t global_sizes[DIMENSIONS] = {cube_size_x, cube_size_y,
                                                 cube_size_z};
              size_t local_coordinates[DIMENSIONS] = {x, y, z};
              size_t global_index = dc_get_global_coordinates(
                  worker_coordinates, worker_sizes, global_sizes,
                  local_coordinates, coordinator_process.topology);
              global_index = dc_get_index_for_coordinates(
                  STENCIL + worker_coordinates[0] * size_x + local_x,
                  STENCIL + worker_coordinates[1] * size_y + local_y,
                  STENCIL + worker_coordinates[2] * size_z + local_z,
                  global_sizes[0], global_sizes[1], global_sizes[2]);
              result.pc[global_index] = pc[worker_index];
              result.qc[global_index] = qc[worker_index];
            }
          }
        }
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
    if (problem_data->precomp_vars_workers[i] != NULL) {
      free(problem_data->precomp_vars_workers[i]->ch1dxx);
      free(problem_data->precomp_vars_workers[i]->ch1dyy);
      free(problem_data->precomp_vars_workers[i]->ch1dzz);
      free(problem_data->precomp_vars_workers[i]->ch1dxy);
      free(problem_data->precomp_vars_workers[i]->ch1dyz);
      free(problem_data->precomp_vars_workers[i]->ch1dxz);
      free(problem_data->precomp_vars_workers[i]->v2px);
      free(problem_data->precomp_vars_workers[i]->v2pz);
      free(problem_data->precomp_vars_workers[i]->v2sz);
      free(problem_data->precomp_vars_workers[i]->v2pn);
      free(problem_data->precomp_vars_workers[i]);
    }
  }
  free(problem_data->pp);
  free(problem_data->qp);
  free(problem_data->pc);
  free(problem_data->qc);
  free(problem_data->source_index);
  free(problem_data->pp_workers);
  free(problem_data->pc_workers);
  free(problem_data->qp_workers);
  free(problem_data->qc_workers);
  free(problem_data->worker_sizes);
  free(problem_data->precomp_vars_workers);

  problem_data->pp_workers = NULL;
  problem_data->qp_workers = NULL;
  problem_data->pc_workers = NULL;
  problem_data->qc_workers = NULL;
  problem_data->worker_sizes = NULL;
  problem_data->precomp_vars_workers = NULL;
}
