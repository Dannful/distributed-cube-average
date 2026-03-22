#pragma once

#include "derivatives.h"
#include "indexing.h"
#include "precomp.h"

// Precomputed constants for theta = atan(1) = pi/4, phi = 1.0
// These are constant across the entire domain
#define CH1DXX_CONST 0.1459632909f  // sin²(π/4) * cos²(1)
#define CH1DYY_CONST 0.3540367091f  // sin²(π/4) * sin²(1)
#define CH1DZZ_CONST 0.5f           // cos²(π/4)
#define CH1DXY_CONST 0.4546487134f  // sin²(π/4) * sin(2)
#define CH1DYZ_CONST 0.8414709848f  // sin(π/2) * sin(1)
#define CH1DXZ_CONST 0.5403023059f  // sin(π/2) * cos(1)

// Constant multipliers for v2px and v2pn (from epsilon=0.24, delta=0.1)
#define V2PX_FACTOR 1.48f  // 1 + 2*epsilon
#define V2PN_FACTOR 1.2f   // 1 + 2*delta

// Optimized version: takes vpz/vsv and computes v2* on-the-fly
static inline HOST_DEVICE void
sample_compute_optimized(size_t x, size_t y, size_t z, size_t size_x, size_t size_y, size_t size_z,
               int process_coord_x, int process_coord_y, int process_coord_z,
               int topology_x, int topology_y, int topology_z, float dx,
               float dy, float dz, float dt, const float *pc, const float *qc,
               float *pp, float *qp, const float *vpz, const float *vsv) {

  // Calculate strides for each dimension
  const int strideX =
      dc_get_index_for_coordinates(1, 0, 0, size_x, size_y, size_z) -
      dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);
  const int strideY =
      dc_get_index_for_coordinates(0, 1, 0, size_x, size_y, size_z) -
      dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);
  const int strideZ =
      dc_get_index_for_coordinates(0, 0, 1, size_x, size_y, size_z) -
      dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);

  // Calculate inverse values for derivatives
  const float dxxinv = 1.0f / (dx * dx);
  const float dyyinv = 1.0f / (dy * dy);
  const float dzzinv = 1.0f / (dz * dz);
  const float dxyinv = 1.0f / (dx * dy);
  const float dxzinv = 1.0f / (dx * dz);
  const float dyzinv = 1.0f / (dy * dz);

  // Calculate index for current position
  const int i = dc_get_index_for_coordinates(x, y, z, size_x, size_y, size_z);

  // Compute v2* values on-the-fly from vpz/vsv
  const float v2pz = vpz[i] * vpz[i];
  const float v2sz = vsv[i] * vsv[i];
  const float v2px = v2pz * V2PX_FACTOR;
  const float v2pn = v2pz * V2PN_FACTOR;

  // p derivatives, H1(p) and H2(p)
  const float pxx = der2(pc, i, strideX, dxxinv);
  const float pyy = der2(pc, i, strideY, dyyinv);
  const float pzz = der2(pc, i, strideZ, dzzinv);
  const float pxy = derCross(pc, i, strideX, strideY, dxyinv);
  const float pyz = derCross(pc, i, strideY, strideZ, dyzinv);
  const float pxz = derCross(pc, i, strideX, strideZ, dxzinv);

  // Use constant ch1d* values
  const float cpxx = CH1DXX_CONST * pxx;
  const float cpyy = CH1DYY_CONST * pyy;
  const float cpzz = CH1DZZ_CONST * pzz;
  const float cpxy = CH1DXY_CONST * pxy;
  const float cpxz = CH1DXZ_CONST * pxz;
  const float cpyz = CH1DYZ_CONST * pyz;
  const float h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
  const float h2p = pxx + pyy + pzz - h1p;

  // q derivatives, H1(q) and H2(q)
  const float qxx = der2(qc, i, strideX, dxxinv);
  const float qyy = der2(qc, i, strideY, dyyinv);
  const float qzz = der2(qc, i, strideZ, dzzinv);
  const float qxy = derCross(qc, i, strideX, strideY, dxyinv);
  const float qyz = derCross(qc, i, strideY, strideZ, dyzinv);
  const float qxz = derCross(qc, i, strideX, strideZ, dxzinv);

  const float cqxx = CH1DXX_CONST * qxx;
  const float cqyy = CH1DYY_CONST * qyy;
  const float cqzz = CH1DZZ_CONST * qzz;
  const float cqxy = CH1DXY_CONST * qxy;
  const float cqxz = CH1DXZ_CONST * qxz;
  const float cqyz = CH1DYZ_CONST * qyz;
  const float h1q = cqxx + cqyy + cqzz + cqxy + cqxz + cqyz;
  const float h2q = qxx + qyy + qzz - h1q;

  // p-q derivatives, H1(p-q)
  const float h1pmq = h1p - h1q;
  const float h2pmq = h2p - h2q;

  // rhs of p and q equations
  float rhsp = v2px * h2p + v2pz * h1q + v2sz * h1pmq;
  float rhsq = v2pn * h2p + v2pz * h1q - v2sz * h2pmq;

  // new p and q
  pp[i] = 2.0f * pc[i] - pp[i] + rhsp * dt * dt;
  qp[i] = 2.0f * qc[i] - qp[i] + rhsq * dt * dt;
}

// Legacy version for backward compatibility (OpenMP uses this)
static inline HOST_DEVICE void
sample_compute(size_t x, size_t y, size_t z, size_t size_x, size_t size_y, size_t size_z,
               int process_coord_x, int process_coord_y, int process_coord_z,
               int topology_x, int topology_y, int topology_z, float dx,
               float dy, float dz, float dt, const float *pc, const float *qc,
               float *pp, float *qp, const float *ch1dxx, const float *ch1dyy,
               const float *ch1dzz, const float *ch1dxy, const float *ch1dyz,
               const float *ch1dxz, const float *v2px, const float *v2pz,
               const float *v2sz, const float *v2pn) {


  // Calculate strides for each dimension
  const int strideX =
      dc_get_index_for_coordinates(1, 0, 0, size_x, size_y, size_z) -
      dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);
  const int strideY =
      dc_get_index_for_coordinates(0, 1, 0, size_x, size_y, size_z) -
      dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);
  const int strideZ =
      dc_get_index_for_coordinates(0, 0, 1, size_x, size_y, size_z) -
      dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);

  // Calculate inverse values for derivatives
  const float dxxinv = 1.0f / (dx * dx);
  const float dyyinv = 1.0f / (dy * dy);
  const float dzzinv = 1.0f / (dz * dz);
  const float dxyinv = 1.0f / (dx * dy);
  const float dxzinv = 1.0f / (dx * dz);
  const float dyzinv = 1.0f / (dy * dz);

  // Calculate index for current position
  const int i = dc_get_index_for_coordinates(x, y, z, size_x, size_y, size_z);

  // p derivatives, H1(p) and H2(p)
  const float pxx = der2(pc, i, strideX, dxxinv);
  const float pyy = der2(pc, i, strideY, dyyinv);
  const float pzz = der2(pc, i, strideZ, dzzinv);
  const float pxy = derCross(pc, i, strideX, strideY, dxyinv);
  const float pyz = derCross(pc, i, strideY, strideZ, dyzinv);
  const float pxz = derCross(pc, i, strideX, strideZ, dxzinv);

  const float cpxx = ch1dxx[i] * pxx;
  const float cpyy = ch1dyy[i] * pyy;
  const float cpzz = ch1dzz[i] * pzz;
  const float cpxy = ch1dxy[i] * pxy;
  const float cpxz = ch1dxz[i] * pxz;
  const float cpyz = ch1dyz[i] * pyz;
  const float h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
  const float h2p = pxx + pyy + pzz - h1p;

  // q derivatives, H1(q) and H2(q)
  const float qxx = der2(qc, i, strideX, dxxinv);
  const float qyy = der2(qc, i, strideY, dyyinv);
  const float qzz = der2(qc, i, strideZ, dzzinv);
  const float qxy = derCross(qc, i, strideX, strideY, dxyinv);
  const float qyz = derCross(qc, i, strideY, strideZ, dyzinv);
  const float qxz = derCross(qc, i, strideX, strideZ, dxzinv);

  const float cqxx = ch1dxx[i] * qxx;
  const float cqyy = ch1dyy[i] * qyy;
  const float cqzz = ch1dzz[i] * qzz;
  const float cqxy = ch1dxy[i] * qxy;
  const float cqxz = ch1dxz[i] * qxz;
  const float cqyz = ch1dyz[i] * qyz;
  const float h1q = cqxx + cqyy + cqzz + cqxy + cqxz + cqyz;
  const float h2q = qxx + qyy + qzz - h1q;

  // p-q derivatives, H1(p-q)
  const float h1pmq = h1p - h1q;
  const float h2pmq = h2p - h2q;

  // rhs of p and q equations
  float rhsp = v2px[i] * h2p + v2pz[i] * h1q +
               v2sz[i] * h1pmq;
  float rhsq = v2pn[i] * h2p + v2pz[i] * h1q -
               v2sz[i] * h2pmq;

  // new p and q
  pp[i] = 2.0f * pc[i] - pp[i] + rhsp * dt * dt;
  qp[i] = 2.0f * qc[i] - qp[i] + rhsq * dt * dt;
}
