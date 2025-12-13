#include "propagate.h"
#include "sample_compute.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#define THREADS_PER_BLOCK 16

__global__ void
propagate_kernel(const size_t *start_coords, const size_t *end_coords,
                 const size_t *sizes, const int *process_coordinates,
                 const int *topology, const dc_precomp_vars *precomp_vars,
                 const float dx, const float dy, const float dz, const float dt,
                 float *pp_out, float *pc, float *qp_out, float *qc,
                 const float *pp_in, const float *qp_in) {
  const size_t x = start_coords[0] + blockIdx.x * blockDim.x + threadIdx.x;
  const size_t y = start_coords[1] + blockIdx.y * blockDim.y + threadIdx.y;

  if (x >= end_coords[0] || y >= end_coords[1]) {
    return;
  }

  for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
    sample_compute(x, y, z, sizes, process_coordinates, topology, dx, dy, dz,
                   dt, pc, qc, pp_in, qp_in, pp_out, qp_out,
                   (dc_precomp_vars *)precomp_vars);
  }
}

extern "C" void dc_propagate(const size_t start_coords[DIMENSIONS],
                             const size_t end_coords[DIMENSIONS],
                             const size_t sizes[DIMENSIONS],
                             const int process_coordinates[DIMENSIONS],
                             const int topology[DIMENSIONS],
                             dc_device_data *data, const float dx,
                             const float dy, const float dz, const float dt) {

  size_t *d_start_coords, *d_end_coords, *d_sizes;
  int *d_process_coordinates, *d_topology;

  cudaMalloc(&d_start_coords, sizeof(size_t) * DIMENSIONS);
  cudaMalloc(&d_end_coords, sizeof(size_t) * DIMENSIONS);
  cudaMalloc(&d_sizes, sizeof(size_t) * DIMENSIONS);
  cudaMalloc(&d_process_coordinates, sizeof(int) * DIMENSIONS);
  cudaMalloc(&d_topology, sizeof(int) * DIMENSIONS);

  cudaMemcpy(d_start_coords, start_coords, sizeof(size_t) * DIMENSIONS,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_end_coords, end_coords, sizeof(size_t) * DIMENSIONS,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_sizes, sizes, sizeof(size_t) * DIMENSIONS,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_process_coordinates, process_coordinates,
             sizeof(int) * DIMENSIONS, cudaMemcpyHostToDevice);
  cudaMemcpy(d_topology, topology, sizeof(int) * DIMENSIONS,
             cudaMemcpyHostToDevice);

  const dim3 threadsPerBlock(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
  const size_t nx = end_coords[0] - start_coords[0];
  const size_t ny = end_coords[1] - start_coords[1];

  const dim3 numBlocks(nx / threadsPerBlock.x, ny / threadsPerBlock.y);

  propagate_kernel<<<numBlocks, threadsPerBlock>>>(
      d_start_coords, d_end_coords, d_sizes, d_process_coordinates, d_topology,
      data->d_precomp_vars, dx, dy, dz, dt, data->pp, data->pc, data->qp,
      data->qc, data->pp_copy, data->qp_copy);

  cudaFree(d_start_coords);
  cudaFree(d_end_coords);
  cudaFree(d_sizes);
  cudaFree(d_process_coordinates);
  cudaFree(d_topology);

  cudaDeviceSynchronize();
}
