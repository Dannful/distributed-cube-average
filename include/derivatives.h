#pragma once

#if defined(__cplusplus)
#include <cstddef>
#else
#include <stddef.h>
#endif

#if defined(__CUDACC__)
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

/* Eight order finite differences coefficients of the first derivative */
static const float L1 = 0.8f;                    // 4/5
static const float L2 = -0.2f;                   // -1/5
static const float L3 = 0.0380952380952381f;     // 4/105
static const float L4 = -0.0035714285714285713f; // -1/280

/* Eight order finite differences coefficients of the cross second derivative */
static const float L11 = 0.64f;  // L1*L1
static const float L12 = -0.16f; // L1*L2
static const float L13 = 0.03047619047619047618f;
static const float L14 = -0.00285714285714285713f;
static const float L22 = 0.04f;
static const float L23 = -0.00761904761904761904f;
static const float L24 = 0.00071428571428571428f;
static const float L33 = 0.00145124716553287981f;
static const float L34 = -0.00013605442176870748f;
static const float L44 = 0.00001275510204081632f;

/* Eight order finite differences coefficients of the second derivative */
static const float K0 = -2.84722222222222222222f; // -205/72
static const float K1 = 1.6f;                     // 8/5
static const float K2 = -0.2f;                    // -1/5
static const float K3 = 0.02539682539682539682f;  // 8/315
static const float K4 = -0.00178571428571428571f; // -1/560

static inline HOST_DEVICE float der1(const float *p, int i, int s, float dinv) {
  return (L1 * (p[i + s] - p[i - s]) + L2 * (p[i + 2 * s] - p[i - 2 * s]) +
          L3 * (p[i + 3 * s] - p[i - 3 * s]) +
          L4 * (p[i + 4 * s] - p[i - 4 * s])) *
         dinv;
}

static inline HOST_DEVICE float der2(const float *p, int i, int s,
                                     float d2inv) {
  return (K0 * p[i] + K1 * (p[i + s] + p[i - s]) +
          K2 * (p[i + 2 * s] + p[i - 2 * s]) +
          K3 * (p[i + 3 * s] + p[i - 3 * s]) +
          K4 * (p[i + 4 * s] + p[i - 4 * s])) *
         d2inv;
}

static inline HOST_DEVICE float derCross(const float *p, int i, int s11,
                                         int s21, float dinv) {
  float term1 = L11 * (p[i + s21 + s11] - p[i + s21 - s11] - p[i - s21 + s11] +
                       p[i - s21 - s11]);
  float term2 = L12 * (p[i + s21 + (2 * s11)] - p[i + s21 - (2 * s11)] -
                       p[i - s21 + (2 * s11)] + p[i - s21 - (2 * s11)] +
                       p[i + (2 * s21) + s11] - p[i + (2 * s21) - s11] -
                       p[i - (2 * s21) + s11] + p[i - (2 * s21) - s11]);
  float term3 = L13 * (p[i + s21 + (3 * s11)] - p[i + s21 - (3 * s11)] -
                       p[i - s21 + (3 * s11)] + p[i - s21 - (3 * s11)] +
                       p[i + (3 * s21) + s11] - p[i + (3 * s21) - s11] -
                       p[i - (3 * s21) + s11] + p[i - (3 * s21) - s11]);
  float term4 = L14 * (p[i + s21 + (4 * s11)] - p[i + s21 - (4 * s11)] -
                       p[i - s21 + (4 * s11)] + p[i - s21 - (4 * s11)] +
                       p[i + (4 * s21) + s11] - p[i + (4 * s21) - s11] -
                       p[i - (4 * s21) + s11] + p[i - (4 * s21) - s11]);
  float term5 =
      L22 * (p[i + (2 * s21) + (2 * s11)] - p[i + (2 * s21) - (2 * s11)] -
             p[i - (2 * s21) + (2 * s11)] + p[i - (2 * s21) - (2 * s11)]);
  float term6 =
      L23 * (p[i + (2 * s21) + (3 * s11)] - p[i + (2 * s21) - (3 * s11)] -
             p[i - (2 * s21) + (3 * s11)] + p[i - (2 * s21) - (3 * s11)] +
             p[i + (3 * s21) + (2 * s11)] - p[i + (3 * s21) - (2 * s11)] -
             p[i - (3 * s21) + (2 * s11)] + p[i - (3 * s21) - (2 * s11)]);
  float term7 =
      L24 * (p[i + (2 * s21) + (4 * s11)] - p[i + (2 * s21) - (4 * s11)] -
             p[i - (2 * s21) + (4 * s11)] + p[i - (2 * s21) - (4 * s11)] +
             p[i + (4 * s21) + (2 * s11)] - p[i + (4 * s21) - (2 * s11)] -
             p[i - (4 * s21) + (2 * s11)] + p[i - (4 * s21) - (2 * s11)]);
  float term8 =
      L33 * (p[i + (3 * s21) + (3 * s11)] - p[i + (3 * s21) - (3 * s11)] -
             p[i - (3 * s21) + (3 * s11)] + p[i - (3 * s21) - (3 * s11)]);
  float term9 =
      L34 * (p[i + (3 * s21) + (4 * s11)] - p[i + (3 * s21) - (4 * s11)] -
             p[i - (3 * s21) + (4 * s11)] + p[i - (3 * s21) - (4 * s11)] +
             p[i + (4 * s21) + (3 * s11)] - p[i + (4 * s21) - (3 * s11)] -
             p[i - (4 * s21) + (3 * s11)] + p[i - (4 * s21) - (3 * s11)]);
  float term10 =
      L44 * (p[i + (4 * s21) + (4 * s11)] - p[i + (4 * s21) - (4 * s11)] -
             p[i - (4 * s21) + (4 * s11)] + p[i - (4 * s21) - (4 * s11)]);
  return (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 +
          term9 + term10) *
         dinv;
}
