#include "precomp.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "coordinator.h"
#include "log.h"

dc_precomp_vars dc_compute_precomp_vars(int sx, int sy, int sz,
                                        dc_anisotropy_t anisotropy) {
  dc_precomp_vars vars = {0};
  int n = sx * sy * sz;

  vars.ch1dxx = (float *)malloc(n * sizeof(float));
  if (vars.ch1dxx == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for ch1dxx in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.ch1dyy = (float *)malloc(n * sizeof(float));
  if (vars.ch1dyy == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for ch1dyy in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.ch1dzz = (float *)malloc(n * sizeof(float));
  if (vars.ch1dzz == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for ch1dzz in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.ch1dxy = (float *)malloc(n * sizeof(float));
  if (vars.ch1dxy == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for ch1dxy in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.ch1dyz = (float *)malloc(n * sizeof(float));
  if (vars.ch1dyz == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for ch1dyz in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.ch1dxz = (float *)malloc(n * sizeof(float));
  if (vars.ch1dxz == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for ch1dxz in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.v2px = (float *)malloc(n * sizeof(float));
  if (vars.v2px == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for v2px in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.v2pz = (float *)malloc(n * sizeof(float));
  if (vars.v2pz == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for v2pz in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.v2sz = (float *)malloc(n * sizeof(float));
  if (vars.v2sz == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for v2sz in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }
  vars.v2pn = (float *)malloc(n * sizeof(float));
  if (vars.v2pn == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for v2pn in "
                              "dc_compute_precomp_vars");
    MPI_Finalize();
    exit(1);
  }

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
    vars.v2px[i] = vars.v2pz[i] * (1.0 + 2.0 * anisotropy.epsilon[i]);
    vars.v2pn[i] = vars.v2pz[i] * (1.0 + 2.0 * anisotropy.delta[i]);
  }

  return vars;
}

dc_anisotropy_t dc_compute_anisotropy_vars(int sx, int sy, int sz) {
  dc_anisotropy_t anisotropy;
  int n = sx * sy * sz;
  anisotropy.vpz = (float *)malloc(sizeof(float) * n);
  if (anisotropy.vpz == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for vpz in "
                              "dc_compute_anisotropy_vars");
    MPI_Finalize();
    exit(1);
  }
  anisotropy.vsv = (float *)malloc(sizeof(float) * n);
  if (anisotropy.vsv == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for vsv in "
                              "dc_compute_anisotropy_vars");
    MPI_Finalize();
    exit(1);
  }
  anisotropy.epsilon = (float *)malloc(sizeof(float) * n);
  if (anisotropy.epsilon == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for epsilon in "
                              "dc_compute_anisotropy_vars");
    MPI_Finalize();
    exit(1);
  }
  anisotropy.delta = (float *)malloc(sizeof(float) * n);
  if (anisotropy.delta == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for delta in "
                              "dc_compute_anisotropy_vars");
    MPI_Finalize();
    exit(1);
  }
  anisotropy.phi = (float *)malloc(sizeof(float) * n);
  if (anisotropy.phi == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for phi in "
                              "dc_compute_anisotropy_vars");
    MPI_Finalize();
    exit(1);
  }
  anisotropy.theta = (float *)malloc(sizeof(float) * n);
  if (anisotropy.theta == NULL) {
    dc_log_error(COORDINATOR, "OOM: could not allocate memory for theta in "
                              "dc_compute_anisotropy_vars");
    MPI_Finalize();
    exit(1);
  }

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

void free_precomp_vars(dc_precomp_vars *precomp) {
  free(precomp->ch1dxx);
  free(precomp->ch1dyy);
  free(precomp->ch1dzz);
  free(precomp->ch1dxy);
  free(precomp->ch1dyz);
  free(precomp->ch1dxz);
  free(precomp->v2px);
  free(precomp->v2pz);
  free(precomp->v2sz);
  free(precomp->v2pn);

  precomp->ch1dxx = NULL;
  precomp->ch1dyy = NULL;
  precomp->ch1dzz = NULL;
  precomp->ch1dxy = NULL;
  precomp->ch1dyz = NULL;
  precomp->ch1dxz = NULL;
  precomp->v2px = NULL;
  precomp->v2pz = NULL;
  precomp->v2sz = NULL;
  precomp->v2pn = NULL;
}
