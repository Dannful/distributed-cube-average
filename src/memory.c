#include "memory.h"
#ifdef SIMGRID
#include <smpi/smpi.h>
#endif
#include <stdlib.h>
#include <string.h>

void *shared_malloc(size_t bytes) {
#ifdef SIMGRID
  return SMPI_SHARED_MALLOC(bytes);
#else
  return malloc(bytes);
#endif
}

void *shared_calloc(size_t n, size_t bytes) {
#ifdef SIMGRID
  void *ptr = SMPI_SHARED_MALLOC(n * bytes);
  memset(ptr, 0, n * bytes);
  return ptr;
#else
  return calloc(n, bytes);
#endif
}
