#include <mpi.h>
#include <string.h>

#include "setup.h"

void dc_mpi_world_init(MPI_Comm *communicator, const int topology[DIMENSIONS]) {
  const int periods[DIMENSIONS] = {0, 0, 0};
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, topology, periods, reorder,
                  communicator);
}

dc_process_t dc_process_init(MPI_Comm communicator, int rank,
                             int topology[DIMENSIONS], size_t sx, size_t sy,
                             size_t sz, float dx, float dy, float dz,
                             float dt) {
  dc_process_t process;
  process.communicator = communicator;
  process.rank = rank;
  process.dx = dx;
  process.dy = dy;
  process.dz = dz;
  process.dt = dt;
  process.source_index = -1;
  process.global_sizes[0] = sx;
  process.global_sizes[1] = sy;
  process.global_sizes[2] = sz;
  memcpy(process.topology, topology, sizeof(int) * DIMENSIONS);
  MPI_Cart_coords(communicator, rank, DIMENSIONS, process.coordinates);

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        size_t face_index = 9 * (dz + 1) + 3 * (dy + 1) + dx + 1;
        if (dx == 0 && dy == 0 && dz == 0) {
          process.neighbours[face_index] = MPI_PROC_NULL;
          continue;
        }

        int displacement[DIMENSIONS] = {dx, dy, dz};
        int target_coords[DIMENSIONS];
        unsigned char is_displacement_valid = 1;
        for (unsigned int i = 0; i < DIMENSIONS; i++) {
          if (process.coordinates[i] + displacement[i] < 0 ||
              process.coordinates[i] + displacement[i] >= process.topology[i]) {
            is_displacement_valid = 0;
            break;
          }
          target_coords[i] = process.coordinates[i] + displacement[i];
        }
        if (!is_displacement_valid) {
          process.neighbours[face_index] = MPI_PROC_NULL;
          continue;
        }
        MPI_Cart_rank(process.communicator, target_coords, process.neighbours + face_index);
      }
    }
  }

  return process;
}
