#pragma once

#include "definitions.h"
#include "precomp.h"
#include <stddef.h>

typedef struct {
  int rank;
  int coordinates[DIMENSIONS];
  int neighbours[NEIGHBOURHOOD];
  int topology[DIMENSIONS];
  unsigned int iterations;
  int source_index;
  float dx, dy, dz, dt;
  size_t sizes[DIMENSIONS];
  dc_anisotropy_t anisotropy_vars;
  dc_precomp_vars precomp_vars;
  float *pp, *pc, *qp, *qc;
  char *hostnames;
  size_t num_workers;
} dc_process_t;
