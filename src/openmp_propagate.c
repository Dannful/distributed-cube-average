#include "precomp.h"
#include "propagate.h"
#include "sample_compute.h"
#include "setup.h"

void dc_propagate(const size_t start_coords[DIMENSIONS],
                  const size_t end_coords[DIMENSIONS],
                  const size_t sizes[DIMENSIONS],
                  const int process_coordinates[DIMENSIONS],
                  const int topology[DIMENSIONS], dc_device_data *data,
                  const float dx, const float dy, const float dz,
                  const float dt) {
#pragma omp parallel for collapse(3)
  for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
    for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
      for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
        sample_compute(x, y, z, sizes, process_coordinates, topology, dx, dy,
                       dz, dt, data->pc, data->qc, data->pp_copy, data->qp_copy,
                       data->pp, data->qp, &data->precomp_vars);
      }
    }
  }
}
