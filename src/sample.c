#include "sample.h"
#include "derivatives.h"
#include "worker.h"
#include "precomp.h"
#include <stddef.h>

void sample_compute(float* pc, float* qc, float* pp, float* qp,
                   dc_precomp_vars precomp_vars,
                   size_t ix, size_t iy, size_t iz,
                   size_t size_x, size_t size_y, size_t size_z,
                   float dx, float dy, float dz, float dt) {
    // Calculate strides for each dimension
    const int strideX = dc_get_index_for_coordinates(1, 0, 0, size_x, size_y, size_z) -
                       dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);
    const int strideY = dc_get_index_for_coordinates(0, 1, 0, size_x, size_y, size_z) -
                       dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);
    const int strideZ = dc_get_index_for_coordinates(0, 0, 1, size_x, size_y, size_z) -
                       dc_get_index_for_coordinates(0, 0, 0, size_x, size_y, size_z);

    // Calculate inverse values for derivatives
    const float dxxinv = 1.0f / (dx * dx);
    const float dyyinv = 1.0f / (dy * dy);
    const float dzzinv = 1.0f / (dz * dz);
    const float dxyinv = 1.0f / (dx * dy);
    const float dxzinv = 1.0f / (dx * dz);
    const float dyzinv = 1.0f / (dy * dz);

    // Calculate index for current position
    const int i = dc_get_index_for_coordinates(ix, iy, iz, size_x, size_y, size_z);

    // p derivatives, H1(p) and H2(p)
    const float pxx = der2(pc, i, strideX, dxxinv);
    const float pyy = der2(pc, i, strideY, dyyinv);
    const float pzz = der2(pc, i, strideZ, dzzinv);
    const float pxy = derCross(pc, i, strideX, strideY, dxyinv);
    const float pyz = derCross(pc, i, strideY, strideZ, dyzinv);
    const float pxz = derCross(pc, i, strideX, strideZ, dxzinv);

    const float cpxx = precomp_vars.ch1dxx[i] * pxx;
    const float cpyy = precomp_vars.ch1dyy[i] * pyy;
    const float cpzz = precomp_vars.ch1dzz[i] * pzz;
    const float cpxy = precomp_vars.ch1dxy[i] * pxy;
    const float cpxz = precomp_vars.ch1dxz[i] * pxz;
    const float cpyz = precomp_vars.ch1dyz[i] * pyz;
    const float h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
    const float h2p = pxx + pyy + pzz - h1p;

    // q derivatives, H1(q) and H2(q)
    const float qxx = der2(qc, i, strideX, dxxinv);
    const float qyy = der2(qc, i, strideY, dyyinv);
    const float qzz = der2(qc, i, strideZ, dzzinv);
    const float qxy = derCross(qc, i, strideX, strideY, dxyinv);
    const float qyz = derCross(qc, i, strideY, strideZ, dyzinv);
    const float qxz = derCross(qc, i, strideX, strideZ, dxzinv);

    const float cqxx = precomp_vars.ch1dxx[i] * qxx;
    const float cqyy = precomp_vars.ch1dyy[i] * qyy;
    const float cqzz = precomp_vars.ch1dzz[i] * qzz;
    const float cqxy = precomp_vars.ch1dxy[i] * qxy;
    const float cqxz = precomp_vars.ch1dxz[i] * qxz;
    const float cqyz = precomp_vars.ch1dyz[i] * qyz;
    const float h1q = cqxx + cqyy + cqzz + cqxy + cqxz + cqyz;
    const float h2q = qxx + qyy + qzz - h1q;

    // p-q derivatives, H1(p-q) and H2(p-q)
    const float h1pmq = h1p - h1q;
    const float h2pmq = h2p - h2q;

    // rhs of p and q equations
    const float rhsp = precomp_vars.v2px[i] * h2p + precomp_vars.v2pz[i] * h1q + precomp_vars.v2sz[i] * h1pmq;
    const float rhsq = precomp_vars.v2pn[i] * h2p + precomp_vars.v2pz[i] * h1q - precomp_vars.v2sz[i] * h2pmq;

    // new p and q
    pp[i] = 2.0f * pc[i] - pp[i] + rhsp * dt * dt;
    qp[i] = 2.0f * qc[i] - qp[i] + rhsq * dt * dt;
}
