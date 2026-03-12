#include "propagate.h"
#include "sample_compute.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#define THREADS_PER_BLOCK 16

__global__ void
propagate_kernel(const size_t start_x, const size_t start_y, const size_t start_z,
                 const size_t end_x, const size_t end_y, const size_t end_z,
                 const size_t size_x, const size_t size_y, const size_t size_z,
                 const int process_coord_x, const int process_coord_y, const int process_coord_z,
                 const int topology_x, const int topology_y, const int topology_z,
                 const float dx, const float dy, const float dz, const float dt,
                 float *pp, float *pc, float *qp, float *qc,
                 const float *ch1dxx, const float *ch1dyy, const float *ch1dzz,
                 const float *ch1dxy, const float *ch1dyz, const float *ch1dxz,
                 const float *v2px, const float *v2pz, const float *v2sz,
                 const float *v2pn) {
  const size_t x = start_x + blockIdx.x * blockDim.x + threadIdx.x;
  const size_t y = start_y + blockIdx.y * blockDim.y + threadIdx.y;

  if (x >= end_x || y >= end_y) {
    return;
  }

  for (size_t z = start_z; z < end_z; z++) {
    sample_compute(x, y, z, size_x, size_y, size_z,
                   process_coord_x, process_coord_y, process_coord_z,
                   topology_x, topology_y, topology_z,
                   dx, dy, dz, dt, pc, qc, pp, qp,
                   ch1dxx, ch1dyy, ch1dzz, ch1dxy, ch1dyz, ch1dxz,
                   v2px, v2pz, v2sz, v2pn);
  }
}

extern "C" void dc_propagate(const size_t start_coords[DIMENSIONS],
                             const size_t end_coords[DIMENSIONS],
                             const size_t sizes[DIMENSIONS],
                             const int process_coordinates[DIMENSIONS],
                             const int topology[DIMENSIONS],
                             dc_device_data *data, const float dx,
                             const float dy, const float dz, const float dt) {

  const dim3 threadsPerBlock(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
  const size_t nx = end_coords[0] - start_coords[0];
  const size_t ny = end_coords[1] - start_coords[1];

  const dim3 numBlocks(nx / threadsPerBlock.x, ny / threadsPerBlock.y);

  propagate_kernel<<<numBlocks, threadsPerBlock>>>(
      start_coords[0], start_coords[1], start_coords[2],
      end_coords[0], end_coords[1], end_coords[2],
      sizes[0], sizes[1], sizes[2],
      process_coordinates[0], process_coordinates[1], process_coordinates[2],
      topology[0], topology[1], topology[2],
      dx, dy, dz, dt, data->pp, data->pc, data->qp, data->qc,
      data->precomp_vars.ch1dxx, data->precomp_vars.ch1dyy,
      data->precomp_vars.ch1dzz, data->precomp_vars.ch1dxy,
      data->precomp_vars.ch1dyz, data->precomp_vars.ch1dxz,
      data->precomp_vars.v2px, data->precomp_vars.v2pz,
      data->precomp_vars.v2sz, data->precomp_vars.v2pn);

  cudaDeviceSynchronize();
}
