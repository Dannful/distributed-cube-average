#pragma once

#include "setup.h"
#include <stddef.h>

void dc_extract_coordinates(size_t *position_x, size_t *position_y, size_t *position_z, size_t size_x, size_t size_y, size_t size_z, int index);
unsigned int dc_get_index_for_coordinates(size_t position_x, size_t position_y, size_t position_z, size_t size_x, size_t size_y, size_t size_z);
void dc_worker_receive_data(dc_process_t *process);
void dc_worker_process(dc_process_t process);
void dc_worker_free(dc_process_t process);
