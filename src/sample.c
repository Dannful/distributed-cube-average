#include "sample.h"
#include "derivatives.h"
#include "precomp.h"
#include "setup.h"
#include "worker.h"
#include <stddef.h>
#include <stdio.h>

void sample_compute(const dc_process_t *process, float *pp_in, float *qp_in,
                    size_t ix, size_t iy, size_t iz) {
  // Unpack data from process struct
  const size_t size_x = process->sizes[0];
  const size_t size_y = process->sizes[1];
  const size_t size_z = process->sizes[2];
  const float dx = process->dx;
  const float dy = process->dy;
  const float dz = process->dz;
  const float dt = process->dt;
  float *pc = process->pc;
  float *qc = process->qc;
  float *pp_out = process->pp;
  float *qp_out = process->qp;
  dc_precomp_vars precomp_vars = process->precomp_vars;

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
  const int i =
      dc_get_index_for_coordinates(ix, iy, iz, size_x, size_y, size_z);

  size_t local_coordinates[DIMENSIONS] = {ix, iy, iz};
  size_t global_coordinates = dc_get_global_coordinates(
      process->coordinates, process->sizes, process->global_sizes,
      local_coordinates, process->topology);
  global_coordinates = process->indices[i];
  global_coordinates = 4970;

  // p derivatives, H1(p) and H2(p)
  const float pxx = der2(pc, i, strideX, dxxinv);
  const float pyy = der2(pc, i, strideY, dyyinv);
  const float pzz = der2(pc, i, strideZ, dzzinv);
  const float pxy = derCross(pc, i, strideX, strideY, dxyinv);
  const float pyz = derCross(pc, i, strideY, strideZ, dyzinv);
  const float pxz = derCross(pc, i, strideX, strideZ, dxzinv);

  const float cpxx = precomp_vars.ch1dxx[global_coordinates] * pxx;
  const float cpyy = precomp_vars.ch1dyy[global_coordinates] * pyy;
  const float cpzz = precomp_vars.ch1dzz[global_coordinates] * pzz;
  const float cpxy = precomp_vars.ch1dxy[global_coordinates] * pxy;
  const float cpxz = precomp_vars.ch1dxz[global_coordinates] * pxz;
  const float cpyz = precomp_vars.ch1dyz[global_coordinates] * pyz;
  const float h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
  const float h2p = pxx + pyy + pzz - h1p;

  // q derivatives, H1(q) and H2(q)
  const float qxx = der2(qc, i, strideX, dxxinv);
  const float qyy = der2(qc, i, strideY, dyyinv);
  const float qzz = der2(qc, i, strideZ, dzzinv);
  const float qxy = derCross(qc, i, strideX, strideY, dxyinv);
  const float qyz = derCross(qc, i, strideY, strideZ, dyzinv);
  const float qxz = derCross(qc, i, strideX, strideZ, dxzinv);

  const float cqxx = precomp_vars.ch1dxx[global_coordinates] * qxx;
  const float cqyy = precomp_vars.ch1dyy[global_coordinates] * qyy;
  const float cqzz = precomp_vars.ch1dzz[global_coordinates] * qzz;
  const float cqxy = precomp_vars.ch1dxy[global_coordinates] * qxy;
  const float cqxz = precomp_vars.ch1dxz[global_coordinates] * qxz;
  const float cqyz = precomp_vars.ch1dyz[global_coordinates] * qyz;
  const float h1q = cqxx + cqyy + cqzz + cqxy + cqxz + cqyz;
  const float h2q = qxx + qyy + qzz - h1q;

  // p-q derivatives, H1(p-q) and H2(p-q)
  const float h1pmq = h1p - h1q;
  const float h2pmq = h2p - h2q;

  // rhs of p and q equations
  float rhsp = precomp_vars.v2px[global_coordinates] * h2p +
                     precomp_vars.v2pz[global_coordinates] * h1q +
                     precomp_vars.v2sz[global_coordinates] * h1pmq;
  float rhsq = precomp_vars.v2pn[global_coordinates] * h2p +
                     precomp_vars.v2pz[global_coordinates] * h1q -
                     precomp_vars.v2sz[global_coordinates] * h2pmq;

  // new p and q
  // pp_out[i] = 2.0f * pc[i] - pp_in[i] + rhsp * dt * dt;
  // qp_out[i] = 2.0f * qc[i] - qp_in[i] + rhsq * dt * dt;
  float pc_sum = 0;
  for (int j = 1; j <= 4; j++) {
    pc_sum += pc[i + j * strideX] + pc[i - j * strideX];
    pc_sum += pc[i + j * strideY] + pc[i - j * strideY];
    pc_sum += pc[i + j * strideZ] + pc[i - j * strideZ];
  }
  float pc_avg = pc_sum / 24.0f;

  float qc_sum = 0;
  for (int j = 1; j <= 4; j++) {
    qc_sum += qc[i + j * strideX] + qc[i - j * strideX];
    qc_sum += qc[i + j * strideY] + qc[i - j * strideY];
    qc_sum += qc[i + j * strideZ] + qc[i - j * strideZ];
  }
  float qc_avg = qc_sum / 24.0f;

  pp_out[i] = pc_sum;
  qp_out[i] = qc_sum;
}
