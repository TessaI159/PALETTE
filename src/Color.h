#ifndef COLOR_H
#define COLOR_H

#include <stdint.h>

struct Color {
  float *alpha;
  float *beta;
  float *gamma;
  uint32_t num;
};

#endif
