#pragma once

#define SIGMA 0.75
#define MAX_SIGMA 10.0

typedef struct {
  float *ch1dxx;
  float *ch1dyy;
  float *ch1dzz;
  float *ch1dxy;
  float *ch1dyz;
  float *ch1dxz;
  float *v2px;
  float *v2pz;
  float *v2sz;
  float *v2pn;
} dc_precomp_vars;

typedef struct {
  float *theta;
  float *phi;
  float *vsv;
  float *vpz;
  float *epsilon;
  float *delta;
} dc_anisotropy_t;

dc_precomp_vars dc_compute_precomp_vars(int sx, int sy, int sz,
                                     dc_anisotropy_t anisotropy,
                                     int bord);
dc_anisotropy_t dc_compute_anisotropy_vars(int sx, int sy, int sz);
void free_precomp_vars(dc_precomp_vars *vars);
void free_anisotropy_vars(dc_anisotropy_t *anisotropy);
