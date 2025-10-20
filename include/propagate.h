#pragma once

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
                  const float *qp_in);
