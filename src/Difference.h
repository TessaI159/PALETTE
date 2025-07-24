#ifndef DIFFERENCE_H
#define DIFFERENCE_H

#include <stdint.h>
typedef struct Color Color;

void delta_diff(const struct Color *restrict colors,
		const struct Color *restrict cents, float *restrict diffs);

#endif
