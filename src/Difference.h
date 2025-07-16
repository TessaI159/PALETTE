#ifndef DIFFERENCE_H
#define DIFFERENCE_H

#include <stdint.h>
typedef struct Color Color;

float *delta_diff(const struct Color colors, const struct Color centroids);

#endif
