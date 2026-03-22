#include <math.h>
#include <stdlib.h>

#include "boundary.h"
#include "indexing.h"

void randomVelocityBoundaryPartition(int local_sx, int local_sy, int local_sz,
                                     int global_sx, int global_sy,
                                     int global_sz, int start_x, int start_y,
                                     int start_z, int nx, int ny, int nz,
                                     int bord, int absorb, float *vpz,
                                     float *vsv, unsigned int *seed) {

  int bordLen = bord + absorb - 1;
  int firstIn = bordLen + 1;
  float frac = 1.0f / (float)(absorb);

  float maxP = 3000.0f;
  float maxS = maxP * sqrtf(fabsf(0.24f - 0.1f) / 0.75f);

  int end_x = start_x + local_sx;
  int end_y = start_y + local_sy;
  int end_z = start_z + local_sz;

  for (int gz = 0; gz < global_sz; gz++) {
    for (int gy = 0; gy < global_sy; gy++) {
      for (int gx = 0; gx < global_sx; gx++) {
        if ((gz >= firstIn && gz <= bordLen + nz) &&
            (gy >= firstIn && gy <= bordLen + ny) &&
            (gx >= firstIn && gx <= bordLen + nx)) {
          continue;
        } else if ((gz >= bord && gz <= 2 * bordLen + nz) &&
                   (gy >= bord && gy <= 2 * bordLen + ny) &&
                   (gx >= bord && gx <= 2 * bordLen + nx)) {
          float rfac = (float)rand_r(seed) / (float)RAND_MAX;

          if (gx >= start_x && gx < end_x && gy >= start_y && gy < end_y &&
              gz >= start_z && gz < end_z) {
            int local_ix = gx - start_x;
            int local_iy = gy - start_y;
            int local_iz = gz - start_z;
            int i = dc_get_index_for_coordinates(local_ix, local_iy, local_iz,
                                                 local_sx, local_sy, local_sz);

            int distz, disty, distx;
            int ivelz, ively, ivelx;

            if (gz > bordLen + nz) {
              distz = gz - bordLen - nz;
              ivelz = bordLen + nz;
            } else if (gz < firstIn) {
              distz = firstIn - gz;
              ivelz = firstIn;
            } else {
              distz = 0;
              ivelz = gz;
            }

            if (gy > bordLen + ny) {
              disty = gy - bordLen - ny;
              ively = bordLen + ny;
            } else if (gy < firstIn) {
              disty = firstIn - gy;
              ively = firstIn;
            } else {
              disty = 0;
              ively = gy;
            }

            if (gx > bordLen + nx) {
              distx = gx - bordLen - nx;
              ivelx = bordLen + nx;
            } else if (gx < firstIn) {
              distx = firstIn - gx;
              ivelx = firstIn;
            } else {
              distx = 0;
              ivelx = gx;
            }

            int dist = (disty > distz) ? disty : distz;
            dist = (dist > distx) ? dist : distx;
            float bordDist = (float)(dist)*frac;

            // Reference velocity: convert global ref coords to local
            int ref_local_x = ivelx - start_x;
            int ref_local_y = ively - start_y;
            int ref_local_z = ivelz - start_z;

            float ref_vpz, ref_vsv;
            if (ref_local_x >= 0 && ref_local_x < local_sx &&
                ref_local_y >= 0 && ref_local_y < local_sy &&
                ref_local_z >= 0 && ref_local_z < local_sz) {
              int ref_i = dc_get_index_for_coordinates(ref_local_x, ref_local_y,
                                                       ref_local_z, local_sx,
                                                       local_sy, local_sz);
              ref_vpz = vpz[ref_i];
              ref_vsv = vsv[ref_i];
            } else {
              ref_vpz = 3000.0f;
              ref_vsv = maxS;
            }

            vpz[i] = ref_vpz * (1.0f - bordDist) + maxP * rfac * bordDist;
            vsv[i] = ref_vsv * (1.0f - bordDist) + maxS * rfac * bordDist;
          }
        } else if (gx >= start_x && gx < end_x && gy >= start_y && gy < end_y &&
                   gz >= start_z && gz < end_z) {
          int local_ix = gx - start_x;
          int local_iy = gy - start_y;
          int local_iz = gz - start_z;
          int i = dc_get_index_for_coordinates(local_ix, local_iy, local_iz,
                                               local_sx, local_sy, local_sz);
          vpz[i] = 0.0f;
          vsv[i] = 0.0f;
        }
      }
    }
  }
}
