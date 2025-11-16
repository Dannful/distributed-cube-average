#include "device_data.h"
#include "indexing.h"
#include "log.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

  data->pp = process->pp;
  data->pc = process->pc;
  data->qp = process->qp;
  data->qc = process->qc;
  data->precomp_vars = process->precomp_vars;

  data->pp_copy = (float *)malloc(total_size_bytes);
  if (data->pp_copy == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for pp_copy in "
                                "dc_device_data_init");
    MPI_Finalize();
    exit(1);
  }
  data->qp_copy = (float *)malloc(total_size_bytes);
  if (data->qp_copy == NULL) {
    dc_log_error(process->rank, "OOM: could not allocate memory for qp_copy in "
                                "dc_device_data_init");
    MPI_Finalize();
    exit(1);
  }

  return data;
}

void dc_device_data_free(dc_device_data *data) {
  free(data->pp_copy);
  free(data->qp_copy);
  free(data);
}

void dc_device_data_get_results(dc_process_t *process, dc_device_data *data) {
  // No-op for OpenMP, as data is already on the host
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
  data->pc[index] += source;
  data->qc[index] += source;
}

void dc_device_extract_halo_face(dc_device_data *data, float *buffer,
                                 const size_t start_coords[DIMENSIONS],
                                 const size_t end_coords[DIMENSIONS],
                                 const size_t sizes[DIMENSIONS],
                                 const float *from_array) {
  size_t data_index = 0;
  for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
    for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
      for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
        size_t from_idx =
            dc_get_index_for_coordinates(x, y, z, sizes[0], sizes[1], sizes[2]);
        buffer[data_index++] = from_array[from_idx];
      }
    }
  }
}

void dc_device_insert_halo_face(dc_device_data *data, const float *buffer,
                                const size_t start_coords[DIMENSIONS],
                                const size_t end_coords[DIMENSIONS],
                                const size_t sizes[DIMENSIONS],
                                float *to_array) {
  size_t data_index = 0;
  for (size_t z = start_coords[2]; z < end_coords[2]; z++) {
    for (size_t y = start_coords[1]; y < end_coords[1]; y++) {
      for (size_t x = start_coords[0]; x < end_coords[0]; x++) {
        to_array[dc_get_index_for_coordinates(x, y, z, sizes[0], sizes[1],
                                              sizes[2])] = buffer[data_index++];
      }
    }
  }
}

void dc_device_data_copy_to_device_copies(dc_device_data *data,
                                          const size_t sizes[DIMENSIONS]) {
  size_t total_size = dc_compute_count_from_sizes((size_t *)sizes);
  size_t total_size_bytes = total_size * sizeof(float);
  memcpy(data->pp_copy, data->pp, total_size_bytes);
  memcpy(data->qp_copy, data->qp, total_size_bytes);
}
