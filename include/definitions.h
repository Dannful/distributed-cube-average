#pragma once

#define DIMENSIONS 3
#define STENCIL 4
#define NEIGHBOURHOOD 27

#if defined(__CUDACC__)
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

#ifdef USE_SIMGRID
#define dc_malloc(bytes) SMPI_SHARED_MALLOC(bytes)
#define dc_calloc(n, bytes) ({ \
  size_t _n = (n); \
  size_t _s = (s); \
  size_t _sz = _n * _s; \
  void *_p = SMPI_SHARED_MALLOC(_sz); \
  if(_p) memset(_p, 0, _sz); \
  _p; \
})
#define dc_free(ptr) SMPI_SHARED_FREE(ptr)
#else
#define dc_malloc(bytes) malloc(bytes)
#define dc_calloc(n, bytes) calloc(n, bytes)
#define dc_free(ptr) free(ptr)
#endif
