#pragma once

#include "setup.h"
#include <stddef.h>

/**
 * Compute one sample iteration for wave propagation.
 *
 * @param sizes Dimensions of the local grid
 * @param dx Grid spacing in x
 * @param dy Grid spacing in y
 * @param dz Grid spacing in z
 * @param dt Timestep
 * @param coordinates Coordinates of the worker in the topology
 * @param global_sizes Dimensions of the global grid
 * @param topology Dimensions of the worker topology
 * @param pc Current p field (p_t)
 * @param qc Current q field (q_t)
 * @param pp_out Next p field (p_{t+1})
 * @param qp_out Next q field (q_{t+1})
 * @param precomp_vars Precomputed variables
 * @param pp_in Previous p field (p_{t-1})
 * @param qp_in Previous q field (q_{t-1})
 * @param ix X coordinate of the point to compute
 * @param iy Y coordinate of the point to compute
 * @param iz Z coordinate of the point to compute
 */
void sample_compute(const size_t sizes[DIMENSIONS], float dx, float dy, float dz,
                    float dt, const int coordinates[DIMENSIONS],
                    const size_t global_sizes[DIMENSIONS],
                    const int topology[DIMENSIONS], float *pc, float *qc,
                    float *pp_out, float *qp_out,
                    dc_precomp_vars precomp_vars, const float *pp_in,
                    const float *qp_in, size_t ix, size_t iy, size_t iz);
