#pragma once

#define FRACABS 0.03125

void randomVelocityBoundaryPartition(
    int local_sx, int local_sy, int local_sz,
    int global_sx, int global_sy, int global_sz,
    int start_x, int start_y, int start_z,
    int nx, int ny, int nz,
    int bord, int absorb,
    float *vpz, float *vsv,
    unsigned int *seed);
