#include <stdio.h>
#ifdef __AVX__
#include <immintrin.h>
#endif

#include "Color.h"

struct Color cielab_avg_avx(const struct Color *colors, uint16_t num_col) {
	puts("AVX average\n");
	return Color_create(0, 0, 0);
}
