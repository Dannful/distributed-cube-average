#include "precomp.h"
#include "propagate.h"
#include "sample.h"
#include "setup.h"

void dc_propagate(const size_t start_coords[DIMENSIONS],
                  const size_t end_coords[DIMENSIONS],
                  const size_t sizes[DIMENSIONS],
                  const int process_coordinates[DIMENSIONS],
                  const size_t global_sizes[DIMENSIONS],
                  const int topology[DIMENSIONS],
                  const dc_precomp_vars *precomp_vars, const float dx,
                  const float dy, const float dz, const float dt, float *pp_out,
                  float *pc, float *qp_out, float *qc, const float *pp_in,
                  const float *qp_in) {
#pragma omp parallel for collapse(3)
  for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
    for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
      for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
        sample_compute(sizes, dx, dy, dz, dt, process_coordinates, global_sizes,
                       topology, pc, qc, pp_out, qp_out, *precomp_vars, pp_in,
                       qp_in, x, y, z);
      }
    }
  }
}
