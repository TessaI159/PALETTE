#include <stdint.h>
#include <stdio.h>

#include "Color.h"

struct Color lab_avg_fallback(struct Color **colors, uint16_t num_col) {
	puts("Fallback average\n");
	return Color_create(0, 0, 0);
}
