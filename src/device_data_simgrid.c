#include "device_data.h"
#include <stdlib.h>
#include <stdio.h>

// SIMGRID mode: No CUDA/GPU operations, everything is simulated via SMPI_SAMPLE_FLOPS

dc_device_data *dc_device_data_init(dc_process_t *process) {
  dc_device_data *data = (dc_device_data *)malloc(sizeof(dc_device_data));
  if (data == NULL) {
    fprintf(stderr, "[%d] OOM: could not allocate device_data\n", process->rank);
    exit(1);
  }
  data->pp = NULL;
  data->pc = NULL;
  data->qp = NULL;
  data->qc = NULL;
  data->precomp_vars.ch1dxx = NULL;
  data->precomp_vars.ch1dyy = NULL;
  data->precomp_vars.ch1dzz = NULL;
  data->precomp_vars.ch1dxy = NULL;
  data->precomp_vars.ch1dyz = NULL;
  data->precomp_vars.ch1dxz = NULL;
  data->precomp_vars.v2px = NULL;
  data->precomp_vars.v2pz = NULL;
  data->precomp_vars.v2sz = NULL;
  data->precomp_vars.v2pn = NULL;
  return data;
}

void dc_device_data_free(dc_device_data *data) {
  free(data);
}

void dc_device_data_get_results(dc_process_t *process, dc_device_data *data) {
  (void)process;
  (void)data;
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

void dc_device_add_source(dc_device_data *data, size_t index, float source) {
  (void)data;
  (void)index;
  (void)source;
}

void dc_device_extract_halo_face(dc_device_data *data, float *buffer,
                                 const size_t start_coords[DIMENSIONS],
                                 const size_t end_coords[DIMENSIONS],
                                 const size_t sizes[DIMENSIONS],
                                 const float *from_array) {
  (void)data;
  (void)buffer;
  (void)start_coords;
  (void)end_coords;
  (void)sizes;
  (void)from_array;
}

void dc_device_insert_halo_face(dc_device_data *data, const float *buffer,
                                const size_t start_coords[DIMENSIONS],
                                const size_t end_coords[DIMENSIONS],
                                const size_t sizes[DIMENSIONS],
                                float *to_array) {
  (void)data;
  (void)buffer;
  (void)start_coords;
  (void)end_coords;
  (void)sizes;
  (void)to_array;
}
