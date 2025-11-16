#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DC_PROPAGATE_H
#define DC_PROPAGATE_H

#include "device_data.h"
#include "definitions.h"

void dc_propagate(const size_t start_coords[DIMENSIONS],
                  const size_t end_coords[DIMENSIONS],
                  const size_t sizes[DIMENSIONS],
                  const int process_coordinates[DIMENSIONS],
                  const int topology[DIMENSIONS], dc_device_data *data,
                  const float dx, const float dy, const float dz,
                  const float dt);

#endif // DC_PROPAGATE_H

#ifdef __cplusplus
}
#endif
