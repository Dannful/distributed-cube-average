#include "../include/utils.h"
#include <stdlib.h>

unsigned int greatest_common_divisor(unsigned int a, unsigned int b) {
  unsigned int d = 0;
    while ((a & 1) == 0 && (b & 1) == 0) {
      a >>= 1;
      b >>= 1;
      d++;
    }
    while ((a & 1) == 0) {
      a >>= 1;
    }
    while ((b & 1) == 0) {
      b >>= 1;
    }
    while(a != b) {
      if (a > b) {
        a -= b;
        while((a & 1) == 0) {
          a >>= 1;
        }
      } else if(b > a) {
        b -= a;
        while((b & 1) == 0) {
          b >>= 1;
        }
      } else {
        exit(EXIT_FAILURE);
      }
    }
    return a << d;
}
