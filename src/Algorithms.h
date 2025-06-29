#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <stdint.h>
typedef struct Color Color;

struct kmeans {
  uint8_t num_centroids;
  float (*diff_func)(Color *sam, Color *ref);
};

#endif
