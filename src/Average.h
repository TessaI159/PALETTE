#ifndef AVERAGE_H
#define AVERAGE_H

#include <stdint.h>

typedef struct Color Color;

void color_avg(const struct Color colors, struct Color cents,
	       uint8_t which_cent);

#endif
