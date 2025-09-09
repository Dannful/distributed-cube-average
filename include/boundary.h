#pragma once

#define FRACABS 0.03125

void randomVelocityBoundary(int sx, int sy, int sz, int nx, int ny, int nz,
                            int bord, int absorb, float *vpz, float *vsv);
