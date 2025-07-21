#pragma once

#include "setup.h"
#include <stddef.h>

void extract_coordinates(size_t *position_x, size_t *position_y, size_t *position_z, size_t size_x, size_t size_y, size_t size_z, int index);
unsigned int get_index_for_coordinates(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z);
unsigned int get_assigned_worker(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z, size_t worker_count);
mpi_process_t receive_worker_data(MPI_Comm communicator, int rank, int topology[DIMENSIONS]);
void worker_process(mpi_process_t process);
void worker_free(mpi_process_t process);
