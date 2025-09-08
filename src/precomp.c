#include "precomp.h"
#include <math.h>
#include <stdlib.h>

dc_precomp_vars dc_compute_precomp_vars(int sx, int sy, int sz,
                                     dc_anisotropy_t anisotropy, int bord) {
  dc_precomp_vars vars = {0};
  int n = sx * sy * sz;

  vars.ch1dxx = (float *)malloc(n * sizeof(float));
  vars.ch1dyy = (float *)malloc(n * sizeof(float));
  vars.ch1dzz = (float *)malloc(n * sizeof(float));
  vars.ch1dxy = (float *)malloc(n * sizeof(float));
  vars.ch1dyz = (float *)malloc(n * sizeof(float));
  vars.ch1dxz = (float *)malloc(n * sizeof(float));
  vars.v2px = (float *)malloc(n * sizeof(float));
  vars.v2pz = (float *)malloc(n * sizeof(float));
  vars.v2sz = (float *)malloc(n * sizeof(float));
  vars.v2pn = (float *)malloc(n * sizeof(float));

  for (int i = 0; i < n; i++) {
    float sinTheta = sin(anisotropy.theta[i]);
    float cosTheta = cos(anisotropy.theta[i]);
    float sin2Theta = sin(2.0f * anisotropy.theta[i]);
    float sinPhi = sin(anisotropy.phi[i]);
    float cosPhi = cos(anisotropy.phi[i]);
    float sin2Phi = sin(2.0f * anisotropy.phi[i]);
    vars.ch1dxx[i] = sinTheta * sinTheta * cosPhi * cosPhi;
    vars.ch1dyy[i] = sinTheta * sinTheta * sinPhi * sinPhi;
    vars.ch1dzz[i] = cosTheta * cosTheta;
    vars.ch1dxy[i] = sinTheta * sinTheta * sin2Phi;
    vars.ch1dyz[i] = sin2Theta * sinPhi;
    vars.ch1dxz[i] = sin2Theta * cosPhi;
  }

  for (int i = 0; i < n; i++) {
    vars.v2sz[i] = anisotropy.vsv[i] * anisotropy.vsv[i];
    vars.v2pz[i] = anisotropy.vpz[i] * anisotropy.vpz[i];
    vars.v2px[i] = vars.v2pz[i] * (1.0f + 2.0f * anisotropy.epsilon[i]);
    vars.v2pn[i] = vars.v2pz[i] * (1.0f + 2.0f * anisotropy.delta[i]);
  }

  return vars;
}

dc_anisotropy_t dc_compute_anisotropy_vars(int sx, int sy, int sz) {
  dc_anisotropy_t anisotropy;
  int n = sx * sy * sz;
  anisotropy.vpz = (float *)malloc(sizeof(float) * n);
  anisotropy.vsv = (float *)malloc(sizeof(float) * n);
  anisotropy.epsilon = (float *)malloc(sizeof(float) * n);
  anisotropy.delta = (float *)malloc(sizeof(float) * n);
  anisotropy.phi = (float *)malloc(sizeof(float) * n);
  anisotropy.theta = (float *)malloc(sizeof(float) * n);

  for (int i = 0; i < n; i++) {
    anisotropy.vpz[i] = 3000.0;
    anisotropy.epsilon[i] = 0.24;
    anisotropy.delta[i] = 0.1;
    anisotropy.phi[i] = 1.0;
    anisotropy.theta[i] = atanf(1.0);
    if (SIGMA > MAX_SIGMA) {
      anisotropy.vsv[i] = 0.0;
    } else {
      anisotropy.vsv[i] =
          anisotropy.vpz[i] *
          sqrtf(fabsf(anisotropy.epsilon[i] - anisotropy.delta[i]) / SIGMA);
    }
  }
  return anisotropy;
}

void free_anisotropy_vars(dc_anisotropy_t *anisotropy) {
  free(anisotropy->theta);
  free(anisotropy->phi);
  free(anisotropy->vsv);
  free(anisotropy->vpz);
  free(anisotropy->epsilon);
  free(anisotropy->delta);

  anisotropy->theta = NULL;
  anisotropy->phi = NULL;
  anisotropy->vsv = NULL;
  anisotropy->vpz = NULL;
  anisotropy->epsilon = NULL;
  anisotropy->delta = NULL;
}
