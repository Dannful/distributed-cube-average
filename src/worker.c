#include <stdlib.h>

#include "../include/worker.h"
#include "../include/utils.h"
#include "../include/setup.h"
#include "../include/io.h"
#include "../include/coordinator.h"
#include "mpi.h"

void extract_coordinates(size_t *position_x, size_t *position_y, size_t *position_z, size_t size_x, size_t size_y, size_t size_z, int index) {
  *position_x = index % size_x;
  *position_y = (index / size_x) % size_y;
  *position_z = index / (size_x * size_y);
}

unsigned int get_index_for_coordinates(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z) {
  return position_x + position_y * size_x + position_z * size_x * size_y;
}

unsigned int get_assigned_worker(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z, size_t worker_count) {
  unsigned int partition_count_x = greatest_common_divisor(size_x, worker_count);
  unsigned int partition_count_y = greatest_common_divisor(size_y, worker_count);
  unsigned int partition_count_z = greatest_common_divisor(size_z, worker_count);
  unsigned char reversed = position_z / size_z % 2;
  unsigned int worker = reversed ? worker_count - (position_x / partition_count_x + position_y / partition_count_y) % worker_count : (position_x / partition_count_x + position_y / partition_count_y) % worker_count;
  return worker;
}

mpi_process_t receive_worker_data(MPI_Comm communicator, int rank, int topology[DIMENSIONS]) {
  mpi_process_t process = mpi_process_init(communicator, rank, topology);
  MPI_Safe_Recv(&process.stencil_size, 1, MPI_UINT32_T, COORDINATOR, MPI_ANY_TAG, communicator);
  MPI_Safe_Recv(&process.iterations, 1, MPI_UINT32_T, COORDINATOR, MPI_ANY_TAG, communicator);
  MPI_Safe_Recv(&process.count, 1, MPI_UNSIGNED_LONG, COORDINATOR, MPI_ANY_TAG, communicator);
  process.indices = malloc(process.count * sizeof(size_t));
  process.data = malloc(process.count * sizeof(float));
  MPI_Safe_Recv(process.indices, process.count, MPI_UNSIGNED_LONG, COORDINATOR, MPI_ANY_TAG, communicator);
  MPI_Safe_Recv(process.data, process.count, MPI_FLOAT, COORDINATOR, MPI_ANY_TAG, communicator);
  return process;
}

void worker_process(mpi_process_t process) {
}

void worker_free(mpi_process_t process) {
  free(process.indices);
  free(process.data);
}
