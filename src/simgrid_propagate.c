#include "propagate.h"

// SIMGRID mode: Computation is modeled via SMPI_SAMPLE_FLOPS in worker.c

void dc_propagate(const size_t start_coords[DIMENSIONS],
                  const size_t end_coords[DIMENSIONS],
                  const size_t sizes[DIMENSIONS],
                  const int process_coordinates[DIMENSIONS],
                  const int topology[DIMENSIONS],
                  dc_device_data *data, const float dx,
                  const float dy, const float dz, const float dt) {
  (void)start_coords;
  (void)end_coords;
  (void)sizes;
  (void)process_coordinates;
  (void)topology;
  (void)data;
  (void)dx;
  (void)dy;
  (void)dz;
  (void)dt;
}
