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
#pragma omp parallel
  {
#pragma omp for
    for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
      for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
        for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
          sample_compute(x, y, z, sizes[0], sizes[1], sizes[2],
                         process_coordinates[0], process_coordinates[1], process_coordinates[2],
                         topology[0], topology[1], topology[2], dx, dy,
                         dz, dt, data->pc, data->qc, data->pp, data->qp,
                         data->precomp_vars.ch1dxx, data->precomp_vars.ch1dyy,
                         data->precomp_vars.ch1dzz, data->precomp_vars.ch1dxy,
                         data->precomp_vars.ch1dyz, data->precomp_vars.ch1dxz,
                         data->precomp_vars.v2px, data->precomp_vars.v2pz,
                         data->precomp_vars.v2sz, data->precomp_vars.v2pn);
        }
      }
    }
  }
}
