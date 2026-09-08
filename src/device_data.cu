#include "dc_process.h"
#include "device_data.h"
#include "indexing.h"
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <stdio.h>
#include <stdlib.h>

static void check_cuda_error(cudaError_t err, int rank, const char *msg) {
  if (err != cudaSuccess) {
    fprintf(stderr, "[%d] CUDA Error: %s - %s\n", rank, msg,
            cudaGetErrorString(err));
    exit(1);
  }
}

int select_device(dc_process_t *process) {
  int device_count;
  check_cuda_error(cudaGetDeviceCount(&device_count), process->rank,
                   "cudaGetDeviceCount");
  const int max_hostname_length = 256;
  char *current_hostname =
      process->hostnames + max_hostname_length * process->rank;
  size_t current_index = 0;
  for (int i = 0; i < process->num_workers && i != process->rank; i++) {
    char *hostname = process->hostnames + max_hostname_length * i;
    if (strcmp(current_hostname, hostname) == 0) {
      current_index++;
    }
  }
  return current_index % device_count;
}

dc_device_data *dc_device_data_init(dc_process_t *process) {
  dc_device_data *data = (dc_device_data *)malloc(sizeof(dc_device_data));
  if (data == NULL) {
    fprintf(stderr,
            "[%d] OOM: could not allocate memory for device_data in "
            "dc_device_data_init\n",
            process->rank);
    exit(1);
  }

  const int device = select_device(process);
  cudaDeviceProp device_prop;
  check_cuda_error(cudaGetDeviceProperties(&device_prop, device), process->rank,
                   "cudaGetDeviceProperties");
  check_cuda_error(cudaSetDevice(device), process->rank, "cudaSetDevice");
  printf("CUDA source using device (%d) %s with compute capability %d.%d\n",
         device, device_prop.name, device_prop.major, device_prop.minor);

  size_t total_size = dc_compute_count_from_sizes(process->sizes);
  size_t total_size_bytes = total_size * sizeof(float);

  // Allocate field arrays (4 arrays)
  check_cuda_error(cudaMalloc(&data->pp, total_size_bytes), process->rank,
                   "cudaMalloc pp");
  check_cuda_error(cudaMalloc(&data->pc, total_size_bytes), process->rank,
                   "cudaMalloc pc");
  check_cuda_error(cudaMalloc(&data->qp, total_size_bytes), process->rank,
                   "cudaMalloc qp");
  check_cuda_error(cudaMalloc(&data->qc, total_size_bytes), process->rank,
                   "cudaMalloc qc");

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

  // Allocate only vpz and vsv (2 arrays instead of 10 precomp arrays)
  // Other values (ch1d*, v2*) are computed on-the-fly in the kernel
  check_cuda_error(cudaMalloc(&data->vpz, total_size_bytes), process->rank,
                   "cudaMalloc vpz");
  check_cuda_error(cudaMalloc(&data->vsv, total_size_bytes), process->rank,
                   "cudaMalloc vsv");

  check_cuda_error(cudaMemcpy(data->vpz, process->anisotropy_vars.vpz,
                              total_size_bytes, cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy vpz");
  check_cuda_error(cudaMemcpy(data->vsv, process->anisotropy_vars.vsv,
                              total_size_bytes, cudaMemcpyHostToDevice),
                   process->rank, "cudaMemcpy vsv");

  return data;
}

void dc_device_data_free(dc_device_data *data) {
  cudaFree(data->pp);
  cudaFree(data->pc);
  cudaFree(data->qp);
  cudaFree(data->qc);
  cudaFree(data->vpz);
  cudaFree(data->vsv);

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
