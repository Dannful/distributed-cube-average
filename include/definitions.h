#pragma once

#define DIMENSIONS 3
#define STENCIL 4
#define NEIGHBOURHOOD 27

#if defined(__CUDACC__)
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
