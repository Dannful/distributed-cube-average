#pragma once

#include <stddef.h>
#include "precomp.h"

/**
 * Compute one sample iteration for wave propagation.
 * 
 * @param pc Current p field
 * @param qc Current q field
 * @param pp Previous p field (will be updated to next)
 * @param qp Previous q field (will be updated to next)
 * @param precomp_vars Precomputed coefficient variables
 * @param ix X coordinate
 * @param iy Y coordinate
 * @param iz Z coordinate
 * @param size_x Grid size in x dimension
 * @param size_y Grid size in y dimension
 * @param size_z Grid size in z dimension
 * @param dx Grid spacing in x dimension
 * @param dy Grid spacing in y dimension
 * @param dz Grid spacing in z dimension
 * @param dt Time step
 */
void sample_compute(float* pc, float* qc, float* pp, float* qp,
                    dc_precomp_vars precomp_vars,
                    size_t ix, size_t iy, size_t iz, 
                    size_t size_x, size_t size_y, size_t size_z,
                    float dx, float dy, float dz, float dt);