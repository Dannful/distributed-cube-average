#pragma once

#include "setup.h"
#include <stddef.h>

/**
 * Compute one sample iteration for wave propagation.
 *
 * @param process The MPI process data structure, containing most inputs.
 * @param pp_in Previous p field (p_{t-1})
 * @param qp_in Previous q field (q_{t-1})
 * @param ix X coordinate of the point to compute
 * @param iy Y coordinate of the point to compute
 * @param iz Z coordinate of the point to compute
 */
void sample_compute(const dc_process_t *process, float *pp_in, float *qp_in,
                    size_t ix, size_t iy, size_t iz);
