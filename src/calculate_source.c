#include <math.h>

#include "calculate_source.h"

#define FCUT 40.0
#define PICUBE 31.00627668029982017537
#define TWOSQRTPI 3.54490770181103205458
#define THREESQRTPI 5.31736155271654808184

float dc_calculate_source(float dt, int it) {
  float tf, fc, fct, expo;
  tf = TWOSQRTPI / FCUT;
  fc = FCUT / THREESQRTPI;
  fct = fc * (((float)it) * dt - tf);
  expo = PICUBE * fct * fct;
  return ((1.0f - 2.0f * expo) * expf(-expo));
}
