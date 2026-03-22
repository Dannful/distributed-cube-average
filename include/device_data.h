#pragma once

#include "dc_process.h"
#include "definitions.h"
#include "precomp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  float *pp, *pc, *qp, *qc;
  float *vpz, *vsv;
  dc_precomp_vars precomp_vars;
} dc_device_data;

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
