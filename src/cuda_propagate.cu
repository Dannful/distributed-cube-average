#include "propagate.h"
#include "sample_compute.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

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

  sample_compute(x, y, z, sizes, global_sizes, process_coordinates, topology,
                 dx, dy, dz, dt, pc, qc, pp_in, qp_in, pp_out, qp_out,
                 (dc_precomp_vars *)precomp_vars);
}

static void
allocate_device_memory(float **d_pp_out, float **d_pc, float **d_qp_out,
                       float **d_qc, float **d_pp_in, float **d_qp_in,
                       size_t total_size, dc_precomp_vars *d_precomp_vars_st,
                       size_t precomp_size, dc_precomp_vars **d_precomp_vars,
                       size_t **d_start_coords, size_t **d_end_coords,
                       size_t **d_sizes, size_t **d_global_sizes,
                       int **d_process_coordinates, int **d_topology) {
  cudaMalloc(d_pp_out, total_size);
  cudaMalloc(d_pc, total_size);
  cudaMalloc(d_qp_out, total_size);
  cudaMalloc(d_qc, total_size);
  cudaMalloc(d_pp_in, total_size);
  cudaMalloc(d_qp_in, total_size);

  cudaMalloc(&d_precomp_vars_st->ch1dxx, precomp_size);
  cudaMalloc(&d_precomp_vars_st->ch1dyy, precomp_size);
  cudaMalloc(&d_precomp_vars_st->ch1dzz, precomp_size);
  cudaMalloc(&d_precomp_vars_st->ch1dxy, precomp_size);
  cudaMalloc(&d_precomp_vars_st->ch1dyz, precomp_size);
  cudaMalloc(&d_precomp_vars_st->ch1dxz, precomp_size);
  cudaMalloc(&d_precomp_vars_st->v2px, precomp_size);
  cudaMalloc(&d_precomp_vars_st->v2pz, precomp_size);
  cudaMalloc(&d_precomp_vars_st->v2sz, precomp_size);
  cudaMalloc(&d_precomp_vars_st->v2pn, precomp_size);

  cudaMalloc(d_precomp_vars, sizeof(dc_precomp_vars));

  cudaMalloc(d_start_coords, sizeof(size_t) * DIMENSIONS);
  cudaMalloc(d_end_coords, sizeof(size_t) * DIMENSIONS);
  cudaMalloc(d_sizes, sizeof(size_t) * DIMENSIONS);
  cudaMalloc(d_global_sizes, sizeof(size_t) * DIMENSIONS);
  cudaMalloc(d_process_coordinates, sizeof(int) * DIMENSIONS);
  cudaMalloc(d_topology, sizeof(int) * DIMENSIONS);
}

static void copy_data_to_device(
    float *d_pc, const float *h_pc, float *d_qc, const float *h_qc,
    float *d_pp_in, const float *h_pp_in, float *d_qp_in, const float *h_qp_in,
    size_t total_size, const dc_precomp_vars *h_precomp_vars,
    dc_precomp_vars *d_precomp_vars_st, size_t precomp_size,
    dc_precomp_vars *d_precomp_vars, const size_t *start_coords,
    size_t *d_start_coords, const size_t *end_coords, size_t *d_end_coords,
    const size_t *sizes, size_t *d_sizes, const size_t *global_sizes,
    size_t *d_global_sizes, const int *process_coordinates,
    int *d_process_coordinates, const int *topology, int *d_topology) {
  cudaMemcpy(d_pc, h_pc, total_size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_qc, h_qc, total_size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_pp_in, h_pp_in, total_size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_qp_in, h_qp_in, total_size, cudaMemcpyHostToDevice);

  cudaMemcpy(d_precomp_vars_st->ch1dxx, h_precomp_vars->ch1dxx, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->ch1dyy, h_precomp_vars->ch1dyy, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->ch1dzz, h_precomp_vars->ch1dzz, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->ch1dxy, h_precomp_vars->ch1dxy, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->ch1dyz, h_precomp_vars->ch1dyz, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->ch1dxz, h_precomp_vars->ch1dxz, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->v2px, h_precomp_vars->v2px, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->v2pz, h_precomp_vars->v2pz, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->v2sz, h_precomp_vars->v2sz, precomp_size,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_precomp_vars_st->v2pn, h_precomp_vars->v2pn, precomp_size,
             cudaMemcpyHostToDevice);

  cudaMemcpy(d_precomp_vars, d_precomp_vars_st, sizeof(dc_precomp_vars),
             cudaMemcpyHostToDevice);

  cudaMemcpy(d_start_coords, start_coords, sizeof(size_t) * DIMENSIONS,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_end_coords, end_coords, sizeof(size_t) * DIMENSIONS,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_sizes, sizes, sizeof(size_t) * DIMENSIONS,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_global_sizes, global_sizes, sizeof(size_t) * DIMENSIONS,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_process_coordinates, process_coordinates,
             sizeof(int) * DIMENSIONS, cudaMemcpyHostToDevice);
  cudaMemcpy(d_topology, topology, sizeof(int) * DIMENSIONS,
             cudaMemcpyHostToDevice);
}

static void copy_data_from_device(float *h_pp_out, const float *d_pp_out,
                                  float *h_qp_out, const float *d_qp_out,
                                  size_t total_size) {
  cudaMemcpy(h_pp_out, d_pp_out, total_size, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_qp_out, d_qp_out, total_size, cudaMemcpyDeviceToHost);
}

static void free_device_memory(float *d_pp_out, float *d_pc, float *d_qp_out,
                               float *d_qc, float *d_pp_in, float *d_qp_in,
                               dc_precomp_vars *d_precomp_vars_st,
                               dc_precomp_vars *d_precomp_vars,
                               size_t *d_start_coords, size_t *d_end_coords,
                               size_t *d_sizes, size_t *d_global_sizes,
                               int *d_process_coordinates, int *d_topology) {
  cudaFree(d_pp_out);
  cudaFree(d_pc);
  cudaFree(d_qp_out);
  cudaFree(d_qc);
  cudaFree(d_pp_in);
  cudaFree(d_qp_in);

  cudaFree(d_precomp_vars_st->ch1dxx);
  cudaFree(d_precomp_vars_st->ch1dyy);
  cudaFree(d_precomp_vars_st->ch1dzz);
  cudaFree(d_precomp_vars_st->ch1dxy);
  cudaFree(d_precomp_vars_st->ch1dyz);
  cudaFree(d_precomp_vars_st->ch1dxz);
  cudaFree(d_precomp_vars_st->v2px);
  cudaFree(d_precomp_vars_st->v2pz);
  cudaFree(d_precomp_vars_st->v2sz);
  cudaFree(d_precomp_vars_st->v2pn);
  cudaFree(d_precomp_vars);

  cudaFree(d_start_coords);
  cudaFree(d_end_coords);
  cudaFree(d_sizes);
  cudaFree(d_global_sizes);
  cudaFree(d_process_coordinates);
  cudaFree(d_topology);
}

extern "C" void dc_propagate(
    const size_t start_coords[DIMENSIONS], const size_t end_coords[DIMENSIONS],
    const size_t sizes[DIMENSIONS], const int process_coordinates[DIMENSIONS],
    const size_t global_sizes[DIMENSIONS], const int topology[DIMENSIONS],
    const dc_precomp_vars *h_precomp_vars, const float dx, const float dy,
    const float dz, const float dt, float *h_pp_out, float *h_pc,
    float *h_qp_out, float *h_qc, const float *h_pp_in, const float *h_qp_in) {

  size_t total_size = sizes[0] * sizes[1] * sizes[2] * sizeof(float);
  size_t precomp_size =
      global_sizes[0] * global_sizes[1] * global_sizes[2] * sizeof(float);

  float *d_pp_out, *d_pc, *d_qp_out, *d_qc, *d_pp_in, *d_qp_in;
  dc_precomp_vars d_precomp_vars_st;
  dc_precomp_vars *d_precomp_vars;
  size_t *d_start_coords, *d_end_coords, *d_sizes, *d_global_sizes;
  int *d_process_coordinates, *d_topology;

  allocate_device_memory(&d_pp_out, &d_pc, &d_qp_out, &d_qc, &d_pp_in, &d_qp_in,
                         total_size, &d_precomp_vars_st, precomp_size,
                         &d_precomp_vars, &d_start_coords, &d_end_coords,
                         &d_sizes, &d_global_sizes, &d_process_coordinates,
                         &d_topology);

  copy_data_to_device(d_pc, h_pc, d_qc, h_qc, d_pp_in, h_pp_in, d_qp_in,
                      h_qp_in, total_size, h_precomp_vars, &d_precomp_vars_st,
                      precomp_size, d_precomp_vars, start_coords,
                      d_start_coords, end_coords, d_end_coords, sizes, d_sizes,
                      global_sizes, d_global_sizes, process_coordinates,
                      d_process_coordinates, topology, d_topology);

  const dim3 threadsPerBlock(16, 16, 16);
  const size_t nx = end_coords[0] - start_coords[0];
  const size_t ny = end_coords[1] - start_coords[1];
  const size_t nz = end_coords[2] - start_coords[2];

  const dim3 numBlocks(nx / threadsPerBlock.x, ny / threadsPerBlock.y,
                       nz / threadsPerBlock.z);

  propagate_kernel<<<numBlocks, threadsPerBlock>>>(
      d_start_coords, d_end_coords, d_sizes, d_process_coordinates,
      d_global_sizes, d_topology, d_precomp_vars, dx, dy, dz, dt, d_pp_out,
      d_pc, d_qp_out, d_qc, d_pp_in, d_qp_in);

  cudaDeviceSynchronize();

  copy_data_from_device(h_pp_out, d_pp_out, h_qp_out, d_qp_out, total_size);

  free_device_memory(d_pp_out, d_pc, d_qp_out, d_qc, d_pp_in, d_qp_in,
                     &d_precomp_vars_st, d_precomp_vars, d_start_coords,
                     d_end_coords, d_sizes, d_global_sizes,
                     d_process_coordinates, d_topology);
}
