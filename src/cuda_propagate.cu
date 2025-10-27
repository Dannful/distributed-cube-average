#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "derivatives.h"
#include "propagate.h"
#include "worker.h"

// Helper macro for CUDA error checking
#define CUDA_CHECK(call)                                                       \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA Error in %s at line %d: %s\n", __FILE__, __LINE__, \
              cudaGetErrorString(err));                                        \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

__global__ void propagate_kernel(
    const size_t *start_coords, const size_t *end_coords, const size_t *sizes,
    const int *process_coordinates, const size_t *global_sizes,
    const int *topology, const dc_precomp_vars *precomp_vars, const float dx,
    const float dy, const float dz, const float dt, float *pp_out, float *pc,
    float *qp_out, float *qc, const float *pp_in, const float *qp_in) {
  const size_t x = start_coords[0] + blockIdx.x * blockDim.x + threadIdx.x;
  const size_t y = start_coords[1] + blockIdx.y * blockDim.y + threadIdx.y;
  const size_t z = start_coords[2] + blockIdx.z * blockDim.z + threadIdx.z;

  if (x >= end_coords[0] || y >= end_coords[1] || z >= end_coords[2]) {
    return;
  }

#include "sample_compute.h"
}

extern "C" void dc_propagate(
    const size_t start_coords[DIMENSIONS], const size_t end_coords[DIMENSIONS],
    const size_t sizes[DIMENSIONS], const int process_coordinates[DIMENSIONS],
    const size_t global_sizes[DIMENSIONS], const int topology[DIMENSIONS],
    const dc_precomp_vars *h_precomp_vars, const float dx, const float dy,
    const float dz, const float dt, float *h_pp_out, float *h_pc,
    float *h_qp_out, float *h_qc, const float *h_pp_in, const float *h_qp_in) {

  size_t wave_size = sizes[0] * sizes[1] * sizes[2] * sizeof(float);
  size_t precomp_size =
      global_sizes[0] * global_sizes[1] * global_sizes[2] * sizeof(float);

  // Device pointers for wavefields
  float *d_pp_out, *d_pc, *d_qp_out, *d_qc, *d_pp_in, *d_qp_in;
  CUDA_CHECK(cudaMalloc(&d_pp_out, wave_size));
  CUDA_CHECK(cudaMalloc(&d_pc, wave_size));
  CUDA_CHECK(cudaMalloc(&d_qp_out, wave_size));
  CUDA_CHECK(cudaMalloc(&d_qc, wave_size));
  CUDA_CHECK(cudaMalloc(&d_pp_in, wave_size));
  CUDA_CHECK(cudaMalloc(&d_qp_in, wave_size));

  // Copy input wavefields to device
  CUDA_CHECK(cudaMemcpy(d_pc, h_pc, wave_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_qc, h_qc, wave_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_pp_in, h_pp_in, wave_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_qp_in, h_qp_in, wave_size, cudaMemcpyHostToDevice));

  // Device pointers for precomp_vars
  dc_precomp_vars d_precomp_vars_st;
  dc_precomp_vars *d_precomp_vars;

  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.ch1dxx, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.ch1dyy, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.ch1dzz, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.ch1dxy, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.ch1dyz, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.ch1dxz, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.v2px, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.v2pz, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.v2sz, precomp_size));
  CUDA_CHECK(cudaMalloc(&d_precomp_vars_st.v2pn, precomp_size));

  // Copy precomp_vars data to device
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.ch1dxx, h_precomp_vars->ch1dxx,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.ch1dyy, h_precomp_vars->ch1dyy,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.ch1dzz, h_precomp_vars->ch1dzz,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.ch1dxy, h_precomp_vars->ch1dxy,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.ch1dyz, h_precomp_vars->ch1dyz,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.ch1dxz, h_precomp_vars->ch1dxz,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.v2px, h_precomp_vars->v2px,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.v2pz, h_precomp_vars->v2pz,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.v2sz, h_precomp_vars->v2sz,
                        precomp_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars_st.v2pn, h_precomp_vars->v2pn,
                        precomp_size, cudaMemcpyHostToDevice));

  // Copy the struct of device pointers to the device
  CUDA_CHECK(cudaMalloc(&d_precomp_vars, sizeof(dc_precomp_vars)));
  CUDA_CHECK(cudaMemcpy(d_precomp_vars, &d_precomp_vars_st,
                        sizeof(dc_precomp_vars), cudaMemcpyHostToDevice));

  // Allocate and copy small configuration arrays
  size_t *d_start_coords, *d_end_coords, *d_sizes, *d_global_sizes;
  int *d_process_coordinates, *d_topology;

  CUDA_CHECK(cudaMalloc(&d_start_coords, sizeof(size_t) * DIMENSIONS));
  CUDA_CHECK(cudaMalloc(&d_end_coords, sizeof(size_t) * DIMENSIONS));
  CUDA_CHECK(cudaMalloc(&d_sizes, sizeof(size_t) * DIMENSIONS));
  CUDA_CHECK(cudaMalloc(&d_global_sizes, sizeof(size_t) * DIMENSIONS));
  CUDA_CHECK(cudaMalloc(&d_process_coordinates, sizeof(int) * DIMENSIONS));
  CUDA_CHECK(cudaMalloc(&d_topology, sizeof(int) * DIMENSIONS));

  CUDA_CHECK(cudaMemcpy(d_start_coords, start_coords,
                        sizeof(size_t) * DIMENSIONS, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_end_coords, end_coords, sizeof(size_t) * DIMENSIONS,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_sizes, sizes, sizeof(size_t) * DIMENSIONS,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_global_sizes, global_sizes,
                        sizeof(size_t) * DIMENSIONS, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_process_coordinates, process_coordinates,
                        sizeof(int) * DIMENSIONS, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_topology, topology, sizeof(int) * DIMENSIONS,
                        cudaMemcpyHostToDevice));

  const dim3 threadsPerBlock(8, 8, 8);
  const size_t nx = end_coords[0] - start_coords[0];
  const size_t ny = end_coords[1] - start_coords[1];
  const size_t nz = end_coords[2] - start_coords[2];

  const dim3 numBlocks((nx + threadsPerBlock.x - 1) / threadsPerBlock.x,
                       (ny + threadsPerBlock.y - 1) / threadsPerBlock.y,
                       (nz + threadsPerBlock.z - 1) / threadsPerBlock.z);

  propagate_kernel<<<numBlocks, threadsPerBlock>>>(
      d_start_coords, d_end_coords, d_sizes, d_process_coordinates,
      d_global_sizes, d_topology, d_precomp_vars, dx, dy, dz, dt, d_pp_out,
      d_pc, d_qp_out, d_qc, d_pp_in, d_qp_in);

  CUDA_CHECK(cudaDeviceSynchronize());

  // Copy results back to host
  CUDA_CHECK(cudaMemcpy(h_pp_out, d_pp_out, wave_size, cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(h_qp_out, d_qp_out, wave_size, cudaMemcpyDeviceToHost));

  // Free device memory
  CUDA_CHECK(cudaFree(d_pp_out));
  CUDA_CHECK(cudaFree(d_pc));
  CUDA_CHECK(cudaFree(d_qp_out));
  CUDA_CHECK(cudaFree(d_qc));
  CUDA_CHECK(cudaFree(d_pp_in));
  CUDA_CHECK(cudaFree(d_qp_in));

  CUDA_CHECK(cudaFree(d_precomp_vars_st.ch1dxx));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.ch1dyy));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.ch1dzz));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.ch1dxy));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.ch1dyz));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.ch1dxz));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.v2px));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.v2pz));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.v2sz));
  CUDA_CHECK(cudaFree(d_precomp_vars_st.v2pn));
  CUDA_CHECK(cudaFree(d_precomp_vars));

  CUDA_CHECK(cudaFree(d_start_coords));
  CUDA_CHECK(cudaFree(d_end_coords));
  CUDA_CHECK(cudaFree(d_sizes));
  CUDA_CHECK(cudaFree(d_global_sizes));
  CUDA_CHECK(cudaFree(d_process_coordinates));
  CUDA_CHECK(cudaFree(d_topology));
}
