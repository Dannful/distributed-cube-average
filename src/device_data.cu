#include "device_data.h"
#include "log.h"
#include "worker.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

static void check_cuda_error(cudaError_t err, int rank, const char *msg) {
  if (err != cudaSuccess) {
    dc_log_error(rank, "CUDA Error: %s - %s", msg, cudaGetErrorString(err));
    MPI_Finalize();
    exit(1);
  }
}

dc_device_data *dc_device_data_init(dc_process_t *process) {
  dc_device_data *data = (dc_device_data *)malloc(sizeof(dc_device_data));
  if (data == NULL) {
    dc_log_error(process->rank,
                 "OOM: could not allocate memory for device_data in "
                 "dc_device_data_init");
    MPI_Finalize();
    exit(1);
  }

  size_t total_size = dc_compute_count_from_sizes(process->sizes);
  size_t total_size_bytes = total_size * sizeof(float);

  check_cuda_error(cudaMalloc(&data->pp, total_size_bytes), process->rank,
                   "cudaMalloc pp");
  check_cuda_error(cudaMalloc(&data->pc, total_size_bytes), process->rank,
                   "cudaMalloc pc");
  check_cuda_error(cudaMalloc(&data->qp, total_size_bytes), process->rank,
                   "cudaMalloc qp");
  check_cuda_error(cudaMalloc(&data->qc, total_size_bytes), process->rank,
                   "cudaMalloc qc");
  check_cuda_error(cudaMalloc(&data->pp_copy, total_size_bytes), process->rank,
                   "cudaMalloc pp_copy");
  check_cuda_error(cudaMalloc(&data->qp_copy, total_size_bytes), process->rank,
                   "cudaMalloc qp_copy");

  check_cuda_error(cudaMemcpy(data->pp, process->pp, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy pp");
  check_cuda_error(cudaMemcpy(data->pc, process->pc, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy pc");
  check_cuda_error(cudaMemcpy(data->qp, process->qp, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy qp");
  check_cuda_error(cudaMemcpy(data->qc, process->qc, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy qc");

  // Allocate and copy precomp_vars
  check_cuda_error(cudaMalloc(&data->precomp_vars.ch1dxx, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.ch1dxx");
  check_cuda_error(cudaMalloc(&data->precomp_vars.ch1dyy, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.ch1dyy");
  check_cuda_error(cudaMalloc(&data->precomp_vars.ch1dzz, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.ch1dzz");
  check_cuda_error(cudaMalloc(&data->precomp_vars.ch1dxy, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.ch1dxy");
  check_cuda_error(cudaMalloc(&data->precomp_vars.ch1dyz, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.ch1dyz");
  check_cuda_error(cudaMalloc(&data->precomp_vars.ch1dxz, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.ch1dxz");
  check_cuda_error(cudaMalloc(&data->precomp_vars.v2px, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.v2px");
  check_cuda_error(cudaMalloc(&data->precomp_vars.v2pz, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.v2pz");
  check_cuda_error(cudaMalloc(&data->precomp_vars.v2sz, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.v2sz");
  check_cuda_error(cudaMalloc(&data->precomp_vars.v2pn, total_size_bytes),
                   process->rank, "cudaMalloc precomp_vars.v2pn");

  check_cuda_error(cudaMemcpy(data->precomp_vars.ch1dxx,
                              process->precomp_vars.ch1dxx, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.ch1dxx");
  check_cuda_error(cudaMemcpy(data->precomp_vars.ch1dyy,
                              process->precomp_vars.ch1dyy, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.ch1dyy");
  check_cuda_error(cudaMemcpy(data->precomp_vars.ch1dzz,
                              process->precomp_vars.ch1dzz, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.ch1dzz");
  check_cuda_error(cudaMemcpy(data->precomp_vars.ch1dxy,
                              process->precomp_vars.ch1dxy, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.ch1dxy");
  check_cuda_error(cudaMemcpy(data->precomp_vars.ch1dyz,
                              process->precomp_vars.ch1dyz, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.ch1dyz");
  check_cuda_error(cudaMemcpy(data->precomp_vars.ch1dxz,
                              process->precomp_vars.ch1dxz, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.ch1dxz");
  check_cuda_error(cudaMemcpy(data->precomp_vars.v2px,
                              process->precomp_vars.v2px, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.v2px");
  check_cuda_error(cudaMemcpy(data->precomp_vars.v2pz,
                              process->precomp_vars.v2pz, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.v2pz");
  check_cuda_error(cudaMemcpy(data->precomp_vars.v2sz,
                              process->precomp_vars.v2sz, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.v2sz");
  check_cuda_error(cudaMemcpy(data->precomp_vars.v2pn,
                              process->precomp_vars.v2pn, total_size_bytes,
                              cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy precomp_vars.v2pn");

  check_cuda_error(cudaMalloc(&data->d_precomp_vars, sizeof(dc_precomp_vars)),
                   process->rank, "cudaMalloc d_precomp_vars");
  check_cuda_error(cudaMemcpy(data->d_precomp_vars, &data->precomp_vars,
                              sizeof(dc_precomp_vars), cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy d_precomp_vars");

  return data;
}

void dc_device_data_free(dc_device_data *data) {
  cudaFree(data->pp);
  cudaFree(data->pc);
  cudaFree(data->qp);
  cudaFree(data->qc);
  cudaFree(data->pp_copy);
  cudaFree(data->qp_copy);

  cudaFree(data->precomp_vars.ch1dxx);
  cudaFree(data->precomp_vars.ch1dyy);
  cudaFree(data->precomp_vars.ch1dzz);
  cudaFree(data->precomp_vars.ch1dxy);
  cudaFree(data->precomp_vars.ch1dyz);
  cudaFree(data->precomp_vars.ch1dxz);
  cudaFree(data->precomp_vars.v2px);
  cudaFree(data->precomp_vars.v2pz);
  cudaFree(data->precomp_vars.v2sz);
  cudaFree(data->precomp_vars.v2pn);
  cudaFree(data->d_precomp_vars);

  free(data);
}

void dc_device_data_get_results(dc_process_t *process, dc_device_data *data) {
  size_t total_size = dc_compute_count_from_sizes(process->sizes);
  size_t total_size_bytes = total_size * sizeof(float);

  check_cuda_error(cudaMemcpy(process->pc, data->pc, total_size_bytes,
                              cudaMemcpyDeviceToHost),
                   process->rank, "cudaMemcpy pc to host");
  check_cuda_error(cudaMemcpy(process->qc, data->qc, total_size_bytes,
                              cudaMemcpyDeviceToHost),
                   process->rank, "cudaMemcpy qc to host");
}

void dc_device_swap_arrays(dc_device_data *data) {
  float *temp;

  temp = data->pp;
  data->pp = data->pc;
  data->pc = temp;

  temp = data->qp;
  data->qp = data->qc;
  data->qc = temp;
}

__global__ void add_source_kernel(float *pc, float *qc, size_t index,
                                  float source) {
  pc[index] += source;
  qc[index] += source;
}

void dc_device_add_source(dc_device_data *data, size_t index, float source) {
  add_source_kernel<<<1, 1>>>(data->pc, data->qc, index, source);
}

void dc_device_extract_halo_face(dc_device_data *data, float *buffer,
                                 const size_t start_coords[DIMENSIONS],
                                 const size_t end_coords[DIMENSIONS],
                                 const size_t sizes[DIMENSIONS],
                                 const float *from_array) {
  size_t width = end_coords[0] - start_coords[0];
  size_t height = end_coords[1] - start_coords[1];
  size_t depth = end_coords[2] - start_coords[2];

  cudaMemcpy3DParms params = {0};
  params.srcPos = make_cudaPos(start_coords[0] * sizeof(float), start_coords[1],
                               start_coords[2]);
  params.dstPos = make_cudaPos(0, 0, 0);
  params.srcPtr = make_cudaPitchedPtr(
      (void *)from_array, sizes[0] * sizeof(float), sizes[0], sizes[1]);
  params.dstPtr =
      make_cudaPitchedPtr(buffer, width * sizeof(float), width, height);
  params.extent = make_cudaExtent(width * sizeof(float), height, depth);
  params.kind = cudaMemcpyDeviceToHost;
  cudaMemcpy3D(&params);
}

void dc_device_insert_halo_face(dc_device_data *data, const float *buffer,
                                const size_t start_coords[DIMENSIONS],
                                const size_t end_coords[DIMENSIONS],
                                const size_t sizes[DIMENSIONS],
                                float *to_array) {
  size_t width = end_coords[0] - start_coords[0];
  size_t height = end_coords[1] - start_coords[1];
  size_t depth = end_coords[2] - start_coords[2];

  cudaMemcpy3DParms params = {0};
  params.srcPos = make_cudaPos(0, 0, 0);
  params.dstPos = make_cudaPos(start_coords[0] * sizeof(float), start_coords[1],
                               start_coords[2]);
  params.srcPtr =
      make_cudaPitchedPtr((void *)buffer, width * sizeof(float), width, height);
  params.dstPtr = make_cudaPitchedPtr(to_array, sizes[0] * sizeof(float),
                                      sizes[0], sizes[1]);
  params.extent = make_cudaExtent(width * sizeof(float), height, depth);
  params.kind = cudaMemcpyHostToDevice;
  cudaMemcpy3D(&params);
}

void dc_device_data_copy_to_device_copies(dc_device_data *data,
                                          const size_t sizes[DIMENSIONS]) {
  size_t total_size = dc_compute_count_from_sizes((size_t *)sizes);
  size_t total_size_bytes = total_size * sizeof(float);
  cudaMemcpy(data->pp_copy, data->pp, total_size_bytes,
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(data->qp_copy, data->qp, total_size_bytes,
             cudaMemcpyDeviceToDevice);
}
