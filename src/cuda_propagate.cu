#include "propagate.h"
#include "sample_compute.h"
#include <cuda_runtime.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define THREADS_PER_BLOCK 16

// The stencil uses 32-bit indices. When the partition has more than INT_MAX
// points (large_partition), the arrays are offset to the current z plane
// (64-bit) and the stencil indexes within it. That path is much slower on
// GPUs, so it is only used when needed.
template <bool large_partition>
__global__ void
propagate_kernel(const size_t start_x, const size_t start_y, const size_t start_z,
                 const size_t end_x, const size_t end_y, const size_t end_z,
                 const size_t size_x, const size_t size_y, const size_t size_z,
                 const int process_coord_x, const int process_coord_y, const int process_coord_z,
                 const int topology_x, const int topology_y, const int topology_z,
                 const float dx, const float dy, const float dz, const float dt,
                 float *pp, float *pc, float *qp, float *qc,
                 const float *vpz, const float *vsv) {
  const size_t x = start_x + blockIdx.x * blockDim.x + threadIdx.x;
  const size_t y = start_y + blockIdx.y * blockDim.y + threadIdx.y;

  if (x >= end_x || y >= end_y) {
    return;
  }

  for (size_t z = start_z; z < end_z; z++) {
    if (large_partition) {
      const size_t plane = z * size_x * size_y;
      sample_compute_optimized(x, y, 0, size_x, size_y, size_z, process_coord_x,
                               process_coord_y, process_coord_z, topology_x,
                               topology_y, topology_z, dx, dy, dz, dt,
                               pc + plane, qc + plane, pp + plane, qp + plane,
                               vpz + plane, vsv + plane);
    } else {
      sample_compute_optimized(x, y, z, size_x, size_y, size_z, process_coord_x,
                               process_coord_y, process_coord_z, topology_x,
                               topology_y, topology_z, dx, dy, dz, dt, pc, qc,
                               pp, qp, vpz, vsv);
    }
  }
}

extern "C" void dc_propagate(const size_t start_coords[DIMENSIONS],
                             const size_t end_coords[DIMENSIONS],
                             const size_t sizes[DIMENSIONS],
                             const int process_coordinates[DIMENSIONS],
                             const int topology[DIMENSIONS],
                             dc_device_data *data, const float dx,
                             const float dy, const float dz, const float dt) {

  if (end_coords[0] <= start_coords[0] ||
      end_coords[1] <= start_coords[1] ||
      end_coords[2] <= start_coords[2]) {
    return;
  }

  const size_t nx = end_coords[0] - start_coords[0];
  const size_t ny = end_coords[1] - start_coords[1];
  const size_t nz = end_coords[2] - start_coords[2];

  const dim3 threadsPerBlock(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
  const dim3 numBlocks((nx + threadsPerBlock.x - 1) / threadsPerBlock.x,
                       (ny + threadsPerBlock.y - 1) / threadsPerBlock.y);

  const bool large_partition = sizes[0] * sizes[1] * sizes[2] > INT_MAX;
  auto kernel = large_partition ? propagate_kernel<true> : propagate_kernel<false>;
  kernel<<<numBlocks, threadsPerBlock>>>(
      start_coords[0], start_coords[1], start_coords[2],
      end_coords[0], end_coords[1], end_coords[2],
      sizes[0], sizes[1], sizes[2],
      process_coordinates[0], process_coordinates[1], process_coordinates[2],
      topology[0], topology[1], topology[2],
      dx, dy, dz, dt, data->pp, data->pc, data->qp, data->qc,
      data->vpz, data->vsv);

  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "CUDA kernel error in dc_propagate: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}
