#include "precomp.h"
#include "propagate.h"
#include "sample_compute.h"
#include "setup.h"

void dc_propagate(const size_t start_coords[DIMENSIONS],
                  const size_t end_coords[DIMENSIONS],
                  const size_t sizes[DIMENSIONS],
                  const int process_coordinates[DIMENSIONS],
                  const int topology[DIMENSIONS],
                  const dc_precomp_vars *precomp_vars, const float dx,
                  const float dy, const float dz, const float dt, float *pp_out,
                  float *pc, float *qp_out, float *qc, const float *pp_in,
                  const float *qp_in) {
#pragma omp parallel for collapse(3)
  for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
    for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
      for (size_t x = 4; x < sizes[0] - 4; x++) {
        sample_compute(x, y, z, sizes, process_coordinates, topology, dx, dy,
                       dz, dt, pc, qc, pp_in, qp_in, pp_out, qp_out,
                       (dc_precomp_vars *)precomp_vars);
      }
    }
  }
}
