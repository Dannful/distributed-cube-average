#ifndef DC_DEVICE_DATA_H
#define DC_DEVICE_DATA_H

#include "dc_process.h"
#include "definitions.h"
#include "precomp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  float *pp, *pc, *qp, *qc;
  float *pp_copy, *qp_copy;
  dc_precomp_vars precomp_vars;

  // CUDA-specific fields for kernel launch parameters
#ifdef __CUDACC__
  size_t *d_start_coords, *d_end_coords, *d_sizes;
  int *d_process_coordinates, *d_topology;
  dc_precomp_vars *d_precomp_vars;
#endif
} dc_device_data;

void dc_device_data_copy_to_device_copies(dc_device_data *data,
                                          const size_t sizes[DIMENSIONS]);

dc_device_data *dc_device_data_init(dc_process_t *process);

void dc_device_data_free(dc_device_data *data);

void dc_device_data_get_results(dc_process_t *process, dc_device_data *data);

void dc_device_swap_arrays(dc_device_data *data);

void dc_device_add_source(dc_device_data *data, size_t index, float source);

void dc_device_extract_halo_face(dc_device_data *data, float *buffer,
                                 const size_t start_coords[DIMENSIONS],
                                 const size_t end_coords[DIMENSIONS],
                                 const size_t sizes[DIMENSIONS],
                                 const float *from_array);

void dc_device_insert_halo_face(dc_device_data *data, const float *buffer,
                                const size_t start_coords[DIMENSIONS],
                                const size_t end_coords[DIMENSIONS],
                                const size_t sizes[DIMENSIONS],
                                float *to_array);

#ifdef __cplusplus
}
#endif

#endif // DC_DEVICE_DATA_H
