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
  const dc_precomp_vars *pv = &data->precomp_vars;
#pragma omp parallel
  {
#pragma omp for
    for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
      // Offset the arrays to the current z plane (64-bit) so the stencil
      // indexes within the plane, which keeps its index arithmetic 32-bit even
      // when the partition has more than 2^31 points.
      const size_t plane = z * sizes[0] * sizes[1];
      for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
        for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
          sample_compute(x, y, 0, sizes[0], sizes[1], sizes[2],
                         process_coordinates[0], process_coordinates[1], process_coordinates[2],
                         topology[0], topology[1], topology[2], dx, dy,
                         dz, dt, data->pc + plane, data->qc + plane,
                         data->pp + plane, data->qp + plane,
                         pv->ch1dxx + plane, pv->ch1dyy + plane,
                         pv->ch1dzz + plane, pv->ch1dxy + plane,
                         pv->ch1dyz + plane, pv->ch1dxz + plane,
                         pv->v2px + plane, pv->v2pz + plane,
                         pv->v2sz + plane, pv->v2pn + plane);
        }
      }
    }
  }
}
